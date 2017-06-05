tic

dist_plotting = 0;

%Airplane, wing, and flight parameters
W = 3200;
b = 34;
c = 5;                  % chord length in feet
S = 170;                %ft^2 (FIXME)
AR = b^2 / S;
e = .8;
Lambda = 0;
h = [0 10000 25000];
rho = 0.00238;

% Conditions
PHAA = [147*1.4667 4.4];    % [speed (ft/s)     loading factor]
PLAA = [270*1.4667 4.4];    % [speed (ft/s)     loading factor]
NHAA = [114*1.4667 -1.76];  % [speed (ft/s)     loading factor]
NLAA = [230*1.4667 -1.76];  % [speed (ft/s)     loading factor]
conditions = [PHAA' PLAA' NHAA' NLAA'];
condition_names = {'PHAA', 'PLAA', 'NHAA', 'NLAA'};

%% Loop through each condition
y = 0:1:12*(b/2);
n = length(y);
for cond_iter = 1:length(conditions)
    % This loop computes everything for a single condition at a time.
    % Might need to change that later.  Anything we need outside of this
    % loop can be saved using eval and condition_names(cond_iter).
    
    %% Lift Distribution
    % So L is in terms of semi span! So everything below (basically) is
    % too!!
    L = W*conditions(2, cond_iter)/2;
    V = conditions(1, cond_iter);
    
    gamma = 4*L/(rho*pi*V*b/2*12);
    L_dist_ellp = fliplr(rho*V*gamma*sqrt(1-(2*y/(12*b)).^2));
    L_dist_rect = fliplr(L/(12*b/2)*ones(1, length(y)));
    L_dist = (L_dist_ellp + L_dist_rect)/2;
    
    
    %% Plotting Lift Distributions
    if dist_plotting
        figure(1)
        plot(y, fliplr(L_dist))
        hold on
        
        if (cond_iter == length(conditions))
            title('Lift Distributions at Sea Level for Four Critical Conditions', 'FontSize', 15)
            xlabel('Spanwise Location [inches]')
            ylabel('Lift [pounds]')
            legend(condition_names)
        end
    end
    
    
   %% Xfoil Data Import
    xfoil_data = csvread(cell2mat(strcat(['xfoil_data/NACA_2412_' condition_names(cond_iter)])));
    xf_alpha = xfoil_data(:,1);
    xf_cl = xfoil_data(:,2);
    xf_cdp = xfoil_data(:,4);
    xf_cm = xfoil_data(:,5);
    
    
    %% Determing Drag
    cL = 2*L/(rho*V^2*S);
    cl = cL;                %FIXME
    alpha = interp1(xf_cl, xf_alpha, cl);
    if alpha >= 0
        indices = xf_alpha>=0;
        xf_alpha_temp = xf_alpha(indices);
        [xf_cdp_temp, indices] = unique(xf_cdp(indices));
        cdp = -1*interp1(xf_alpha_temp(indices), (-1*xf_cdp_temp), alpha);
    else
        indicies = xf_alpha < 0;
        xf_alpha_temp = xf_alpha(indices);
        [xf_cdp_temp, indices] = unique(xf_cdp(indices));
        cdp = interp1(xf_alpha_temp(indices), xf_cdp_temp, alpha);
    end
    cDp = cdp;              %FIXME?
    cDi = cL^2 / (pi * AR * e);
    cD = (cDp + cDi);
    D = 1/2*(rho*V^2*cD*S);
    D_dist_rect = ones(1,length(y));
    D_dist_tip = D_dist_rect * 1.25;
    D_dist = fliplr(horzcat(D_dist_rect(1:uint8(.8*length(D_dist_rect))), D_dist_tip(uint8(.8*length(D_dist_tip))+1:length(D_dist_tip))));
    D_dist = (D/(12*b/2))*D_dist./mean(D_dist);

    if dist_plotting
        figure(2);
        plot(y,fliplr(D_dist))
        hold on
        if (cond_iter == length(conditions))
            title('Drag Distributions at Sea Level for Four Critical Conditions', 'FontSize', 15)
            xlabel('Spanwise Location [inches]')
            ylabel('Drag Force [pounds/in]')
            legend(condition_names)
        end
        if (cond_iter == 1)
            figure(10) 
            plot(y, fliplr(L_dist_ellp))
            hold on
            plot(y, fliplr(L_dist_rect))
            plot(y, fliplr(L_dist))
            title('Elliptical, Rectangular, and Average Lift Distributions', 'FontSize', 15)
            xlabel('Spanwise Location [inches]')
            ylabel('Drag Force [pounds/in]')
            legend('Elliptical','Rectangular','Average')
        end
    end
    
    
    %% Wing coordinate transformation
    alpha = alpha*pi/180;
    XLoad = -L_dist * sin(alpha) + D_dist * cos(alpha);
    ZLoad = L_dist * cos(alpha) + D_dist * sin(alpha);
    
    %% Shear Flow
    Vx = numInt(XLoad);
    Vz = numInt(ZLoad);
    eval(cell2mat(strcat([condition_names(cond_iter) '_VX = Vx;'])));
    eval(cell2mat(strcat([condition_names(cond_iter) '_VZ = Vz;'])));
    
    %% Moment Dis
    Mx = numInt(Vx);
    Mz = numInt(Vz);
    eval(cell2mat(strcat([condition_names(cond_iter) '_MX = Mx;'])));
    eval(cell2mat(strcat([condition_names(cond_iter) '_MZ = Mz;'])));
    
    
    %% Plotting
    if dist_plotting
    figure(3)
    subplot(3,1,1)
    hold on
    plot(y, fliplr(XLoad))
    subplot(3,1,2)
    hold on
    plot(y, fliplr(Vx))
    subplot(3,1,3)
    hold on
    plot(y, fliplr(Mx))
    
    figure(4)
    subplot(3,1,1)
    hold on
    plot(y, fliplr(ZLoad))
    subplot(3,1,2)
    hold on
    plot(y, fliplr(Vz))
    subplot(3,1,3)
    hold on
    plot(y, fliplr(Mz))
    
    if cond_iter == length(conditions)
        figure(3)
        subplot(3,1,1)
        title('X-Direction Load Distribution', 'FontSize', 12)
        xlabel('Spanwise Location [inches]')
        ylabel('Pounds/Inch')
        legend(condition_names)
        
        subplot(3,1,2)
        title('X-Direction Shear Flow Distribution', 'FontSize', 12)
        xlabel('Spanwise Location [inches]')
        ylabel('Pounds/Inch^2')
        legend(condition_names)
        
        subplot(3,1,3)
        title('X-Direction Moment Distribution', 'FontSize', 12)
        xlabel('Spanwise Location [inches]')
        ylabel('Pounds/Inch^3')
        legend(condition_names)
        
        
        figure(4)
        subplot(3,1,1)
        title('N-Direction Load Distribution', 'FontSize', 12)
        xlabel('Spanwise Location [inches]')
        ylabel('Pounds/Inch')
        legend(condition_names)
        
        subplot(3,1,2)
        title('N-Direction Shear Flow Distribution', 'FontSize', 12)
        xlabel('Spanwise Location [inches]')
        ylabel('Pounds/Inch^2')
        legend(condition_names)
        
        subplot(3,1,3)
        title('N-Direction Moment Distribution', 'FontSize', 12)
        xlabel('Spanwise Location [inches]')
        ylabel('Pounds/Inch^3')
        legend(condition_names)
        
    end
       
    end 
end

%% Determining cL
% L = 1/2 rho v^2 S cL
% cL = 2*L/(rho*v^2*s);


%% DRAG DISTRIBUTION
% Make this rectangular with 100% from 0-80, 125% from 80-100
% figure(2)

