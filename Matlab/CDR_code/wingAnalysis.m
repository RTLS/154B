%wing shear flow
clear all;
close all;

distributions
plotting = 0;


%% Define our geometry -- will optimize on this later
numTopStringers = 4;
numBottomStringers = 4;
numNoseTopStringers = 2;
numNoseBottomStringers = 2;

% Thicknesses in inches
t_upper = 0.02;
t_lower = 0.02;
t_upper_front = 0.02;
t_lower_front = 0.02;
t_frontSpar = 0.04;
t_rearSpar = 0.04;

c = 5*12;

frontSpar = 0.2 * c;
backSpar = 0.7 * c;

numRibs = 18;
ribL = b*12 / (2*(numRibs-1));

%% Solve for our conditions in unit cases, for superposition

for loading_condition = 1:3
   
% What are the units here?
if loading_condition == 1
    Vx = 1; Vz = 0; My = 0;
elseif loading_condition == 2
    Vx = 0; Vz = 1; My = 0;
elseif loading_condition == 3
    Vx = 0; Vz = 0; My = 1;
end

sparCaps(1).posX = frontSpar;
sparCaps(2).posX = frontSpar;
sparCaps(3).posX = backSpar;
sparCaps(4).posX = backSpar;

sparCaps(1).posZ = get_z(frontSpar,1);
sparCaps(2).posZ = get_z(frontSpar,0);
sparCaps(3).posZ = get_z(backSpar,1);
sparCaps(4).posZ = get_z(backSpar,0);

% in in^2
sparCaps(1).area = .1;
sparCaps(2).area = .1;
sparCaps(3).area = .1;
sparCaps(4).area = .1;

upperStringerGap = (sparCaps(3).posX - sparCaps(1).posX)/(numTopStringers + 1);
lowerStringerGap = (sparCaps(3).posX - sparCaps(1).posX)/(numBottomStringers + 1);
upperNoseStringerGap = (sparCaps(1).posX - 0)/(numNoseTopStringers + 1);
lowerNoseStringerGap = (sparCaps(1).posX - 0)/(numNoseBottomStringers + 1);


%assume stringers spaced evenly along X axis betwen Spars
%top Stringers
for i=1:numTopStringers
    topStringers(i).posX = sparCaps(1).posX + upperStringerGap*i;
    topStringers(i).posZ = get_z(topStringers(i).posX,1);
    topStringers(i).area = .1;
end

%bottom Stringers
for i=1:numBottomStringers
    bottomStringers(i).posX = sparCaps(4).posX - lowerStringerGap*i;
    bottomStringers(i).posZ = get_z(bottomStringers(i).posX,0);
    bottomStringers(i).area = .1;

end

%nose bottom Stringers
for i=1:numNoseBottomStringers
    noseBottomStringers(i).posX = sparCaps(2).posX - lowerNoseStringerGap*i;
    noseBottomStringers(i).posZ = get_z(noseBottomStringers(i).posX,0);
    noseBottomStringers(i).area = .1;
end

%nose top Stringers
for i=1:numNoseTopStringers
    noseTopStringers(i).posX = upperNoseStringerGap*i;
    noseTopStringers(i).posZ = get_z(noseTopStringers(i).posX,1);
    noseTopStringers(i).area = .1;
end


centroid.posX = sum([sparCaps.posX].*[sparCaps.area]) + ...
    sum([topStringers.posX].*[topStringers.area]) + ...
    sum([bottomStringers.posX].*[bottomStringers.area]) + ...
    sum([noseTopStringers.posX].*[noseTopStringers.area]) + ...
    sum([noseBottomStringers.posX].*[noseBottomStringers.area]);

centroid.posX = centroid.posX / ( sum([sparCaps.area]) + sum([topStringers.area]) + ...
    sum([bottomStringers.area]) + sum([noseTopStringers.area]) + sum([noseBottomStringers.area]));

centroid.posZ = sum([sparCaps.posZ].*[sparCaps.area]) + ...
    sum([topStringers.posZ].*[topStringers.area]) + ...
    sum([bottomStringers.posZ].*[bottomStringers.area]) + ...
    sum([noseTopStringers.posZ].*[noseTopStringers.area]) + ...
    sum([noseBottomStringers.posZ].*[noseBottomStringers.area]);
centroid.posZ = centroid.posZ / ( sum([sparCaps.area]) + sum([topStringers.area]) + ...
    sum([bottomStringers.area]) + sum([noseTopStringers.area]) + sum([noseBottomStringers.area]));

%summing contributions for inertia terms
Ix = 0; Iz = 0; Ixz = 0;

for i=1:4 %spar caps
    Ix = Ix + sparCaps(i).area*(sparCaps(i).posZ-centroid.posZ)^2;
    Iz = Iz + sparCaps(i).area*(sparCaps(i).posX-centroid.posX)^2;
    Ixz = Ixz + sparCaps(i).area*(sparCaps(i).posX-centroid.posX)*(sparCaps(i).posZ-centroid.posZ);
end


for i=1:numTopStringers %top stringers
    Ix = Ix + topStringers(i).area*(topStringers(i).posZ-centroid.posZ)^2;
    Iz = Iz + topStringers(i).area*(topStringers(i).posX-centroid.posX)^2;
    Ixz = Ixz + topStringers(i).area*(topStringers(i).posX-centroid.posX)*(topStringers(i).posZ-centroid.posZ);
end
for i=1:numBottomStringers %bottom stringers
    Ix = Ix + bottomStringers(i).area*(bottomStringers(i).posZ-centroid.posZ)^2;
    Iz = Iz + bottomStringers(i).area*(bottomStringers(i).posX-centroid.posX)^2;
    Ixz = Ixz + bottomStringers(i).area*(bottomStringers(i).posX-centroid.posX)*(bottomStringers(i).posZ-centroid.posZ);
end
for i=1:numNoseTopStringers %nose top stringers
    Ix = Ix + noseTopStringers(i).area*(noseTopStringers(i).posZ-centroid.posZ)^2;
    Iz = Iz + noseTopStringers(i).area*(noseTopStringers(i).posX-centroid.posX)^2;
    Ixz = Ixz + noseTopStringers(i).area*(noseTopStringers(i).posX-centroid.posX)*(noseTopStringers(i).posZ-centroid.posZ);
end
for i=1:numNoseBottomStringers %nose bottom stringers
    Ix = Ix + noseBottomStringers(i).area*(noseBottomStringers(i).posZ-centroid.posZ)^2;
    Iz = Iz + noseBottomStringers(i).area*(noseBottomStringers(i).posX-centroid.posX)^2;
    Ixz = Ixz + noseBottomStringers(i).area*(noseBottomStringers(i).posX-centroid.posX)*(noseBottomStringers(i).posZ-centroid.posZ);
end

%Ixz = -Ixz;

%define webs

%% web cell 1

%upper webs
numWebs = numTopStringers;
stringerGap = upperStringerGap;
webThickness = t_upper;
tempStringers = topStringers;

for i=1:(numWebs+1)
    web(i).xStart = sparCaps(1).posX + stringerGap*(i-1);
    web(i).xEnd = sparCaps(1).posX + stringerGap*(i);
    web(i).thickness = webThickness;
    web(i).zStart = get_z(web(i).xStart,1);
    web(i).zEnd = get_z(web(i).xEnd,1);
    if i==1
        web(i).dp_area = sparCaps(1).area;
        web(i).dP_X = 0;
        web(i).dP_Z = 0;
        web(i).qPrime_X = 0;
        web(i).qPrime_Z = 0;
    else
        web(i).dp_area = tempStringers(i-1).area;
        dx = web(i).xStart-centroid.posX; dz = web(i).zStart-centroid.posZ;
        web(i).dP_X = get_dp(dx,dz,Vx,0,Ix,Iz,Ixz,web(i).dp_area);  %just Vx
        web(i).dP_Z = get_dp(dx,dz,0,Vz,Ix,Iz,Ixz,web(i).dp_area);  %just Vz
        web(i).qPrime_X = web(i-1).qPrime_X - web(i).dP_X;
        web(i).qPrime_Z = web(i-1).qPrime_Z - web(i).dP_Z;
    end
    tempInt = get_int(web(i).xStart,web(i).xEnd,1);  %integral of airfoil function
    triangle1 = abs( (web(i).xStart - sparCaps(1).posX)*web(i).zStart/2);
    triangle2 = abs((web(i).xEnd - sparCaps(1).posX)*web(i).zEnd/2);
    web(i).Area = tempInt + triangle1 - triangle2;
    web(i).ds = get_ds(web(i).xStart,web(i).xEnd,1);
    web(i).dS_over_t = web(i).ds / web(i).thickness;
    
    web(i).q_dS_over_t_X = web(i).qPrime_X * web(i).dS_over_t;
    web(i).q_dS_over_t_Z = web(i).qPrime_Z * web(i).dS_over_t;
    web(i).two_A_qprime_X = 2*web(i).Area*web(i).qPrime_X;
    web(i).two_A_qprime_Z = 2*web(i).Area*web(i).qPrime_Z;
    web(i).qp_dx_X = web(i).qPrime_X *(web(i).xEnd-web(i).xStart);
    web(i).qp_dx_Z = web(i).qPrime_Z *(web(i).xEnd-web(i).xStart);
    web(i).qp_dz_X = web(i).qPrime_X *(web(i).zEnd-web(i).zStart);
    web(i).qp_dz_Z = web(i).qPrime_Z *(web(i).zEnd-web(i).zStart);
    %web(i).radCurv = ...   Example:  get_curve(web(i).xStart,web(i).xEnd,1)
end
webTop = web;
web = [];

%rear spar
i=1;
web(i).xStart = sparCaps(3).posX;
web(i).xEnd = sparCaps(4).posX;
web(i).thickness = t_rearSpar;
web(i).zStart = sparCaps(3).posZ;
web(i).zEnd = sparCaps(4).posZ;
web(i).dp_area = sparCaps(3).area;
dx = web(i).xStart-centroid.posX; dz = web(i).zStart-centroid.posZ;
web(i).dP_X = get_dp(dx,dz,Vx,0,Ix,Iz,Ixz,web(i).dp_area);
web(i).dP_Z = get_dp(dx,dz,0,Vz,Ix,Iz,Ixz,web(i).dp_area);
web(i).qPrime_X = webTop(numTopStringers+1).qPrime_X - web(i).dP_X;
web(i).qPrime_Z = webTop(numTopStringers+1).qPrime_Z - web(i).dP_Z;

web(i).Area = (sparCaps(3).posX-sparCaps(1).posX)*sparCaps(3).posZ/2 + ...
    abs((sparCaps(3).posX-sparCaps(1).posX)*sparCaps(4).posZ/2);
web(i).ds = abs(sparCaps(3).posZ - sparCaps(4).posZ);
web(i).dS_over_t = web(i).ds / web(i).thickness;

web(i).q_dS_over_t_X = web(i).qPrime_X * web(i).dS_over_t;
web(i).q_dS_over_t_Z = web(i).qPrime_Z * web(i).dS_over_t;
web(i).two_A_qprime_X = 2*web(i).Area*web(i).qPrime_X;
web(i).two_A_qprime_Z = 2*web(i).Area*web(i).qPrime_Z;
web(i).qp_dx_X = web(i).qPrime_X *(web(i).xEnd-web(i).xStart);
web(i).qp_dx_Z = web(i).qPrime_Z *(web(i).xEnd-web(i).xStart);
web(i).qp_dz_X = web(i).qPrime_X *(web(i).zEnd-web(i).zStart);
web(i).qp_dz_Z = web(i).qPrime_Z *(web(i).zEnd-web(i).zStart);

webRearSpar = web;
web = [];


%lower webs
numWebs = numBottomStringers;
stringerGap = lowerStringerGap;
webThickness = t_lower;
tempStringers = bottomStringers;

for i=1:(numWebs+1)
    web(i).xStart = sparCaps(4).posX - stringerGap*(i-1);
    web(i).xEnd = sparCaps(4).posX - stringerGap*(i);
    web(i).thickness = webThickness;
    web(i).zStart = get_z(web(i).xStart,0);
    web(i).zEnd = get_z(web(i).xEnd,0);
    dx = web(i).xStart-centroid.posX; dz = web(i).zStart-centroid.posZ;
    if i==1
        web(i).dp_area = sparCaps(4).area;
        web(i).dP_X = get_dp(dx,dz,Vx,0,Ix,Iz,Ixz,web(i).dp_area);
        web(i).dP_Z = get_dp(dx,dz,0,Vz,Ix,Iz,Ixz,web(i).dp_area);
        web(i).qPrime_X = webRearSpar.qPrime_X - web(i).dP_X;
        web(i).qPrime_Z = webRearSpar.qPrime_Z - web(i).dP_Z;
    else
        web(i).dp_area = tempStringers(i-1).area;
        web(i).dP_X = get_dp(dx,dz, Vx,0,Ix,Iz,Ixz,web(i).dp_area);
        web(i).dP_Z = get_dp(dx,dz, 0,Vz,Ix,Iz,Ixz,web(i).dp_area);
        web(i).qPrime_X = web(i-1).qPrime_X - web(i).dP_X;
        web(i).qPrime_Z = web(i-1).qPrime_Z - web(i).dP_Z;
    end
    
    tempInt = get_int(web(i).xEnd,web(i).xStart,0);  %integral of airfoil function
    triangle2 = abs((web(i).xStart - sparCaps(1).posX)*web(i).zStart/2);
    triangle1 = abs((web(i).xEnd - sparCaps(1).posX)*web(i).zEnd/2);
    web(i).Area = tempInt + triangle1 - triangle2;
    web(i).ds = get_ds(web(i).xStart,web(i).xEnd,0);
    web(i).dS_over_t = web(i).ds / web(i).thickness;

    web(i).q_dS_over_t_X = web(i).qPrime_X * web(i).dS_over_t;
    web(i).q_dS_over_t_Z = web(i).qPrime_Z * web(i).dS_over_t;
    web(i).two_A_qprime_X = 2*web(i).Area*web(i).qPrime_X;
    web(i).two_A_qprime_Z = 2*web(i).Area*web(i).qPrime_Z;
    web(i).qp_dx_X = web(i).qPrime_X*(web(i).xEnd-web(i).xStart);
    web(i).qp_dx_Z = web(i).qPrime_Z*(web(i).xEnd-web(i).xStart);
    web(i).qp_dz_X = web(i).qPrime_X*(web(i).zEnd-web(i).zStart);
    web(i).qp_dz_Z = web(i).qPrime_Z*(web(i).zEnd-web(i).zStart);

    %web(i).radCurv = ...   Example:  get_curve(web(i).xStart,web(i).xEnd,1)
end
webBottom = web;
web = [];

%front Spar
i=1;
web(i).xStart = sparCaps(2).posX;
web(i).xEnd = sparCaps(1).posX;
web(i).thickness = t_frontSpar;
web(i).zStart = sparCaps(2).posZ;
web(i).zEnd = sparCaps(1).posZ;
web(i).dp_area = sparCaps(2).area;
dx = web(i).xStart-centroid.posX; dz = web(i).zStart-centroid.posZ;
web(i).dP_X = get_dp(dx,dz,Vx,0,Ix,Iz,Ixz,web(i).dp_area);
web(i).dP_Z = get_dp(dx,dz,0,Vz,Ix,Iz,Ixz,web(i).dp_area);
web(i).qPrime_X = webBottom(numBottomStringers+1).qPrime_X - web(i).dP_X;
web(i).qPrime_Z = webBottom(numBottomStringers+1).qPrime_Z - web(i).dP_Z;
web(i).Area = 0;
web(i).ds = abs(sparCaps(2).posZ - sparCaps(1).posZ);
web(i).dS_over_t = web(i).ds / web(i).thickness;

web(i).q_dS_over_t_X = web(i).qPrime_X * web(i).dS_over_t;
web(i).q_dS_over_t_Z = web(i).qPrime_Z * web(i).dS_over_t;
web(i).two_A_qprime_X = 2*web(i).Area*web(i).qPrime_X;
web(i).two_A_qprime_Z = 2*web(i).Area*web(i).qPrime_Z;
web(i).qp_dx_X = web(i).qPrime_X *(web(i).xEnd-web(i).xStart);
web(i).qp_dx_Z = web(i).qPrime_Z *(web(i).xEnd-web(i).xStart);
web(i).qp_dz_X = web(i).qPrime_X *(web(i).zEnd-web(i).zStart);
web(i).qp_dz_Z = web(i).qPrime_Z *(web(i).zEnd-web(i).zStart);

webFrontSpar = web;
web = [];


%% web cell 2

%lower nose webs
numWebs = numNoseBottomStringers;
stringerGap = lowerNoseStringerGap;
webThickness = t_lower_front;
tempStringers = noseBottomStringers;

for i=1:(numWebs+1)
    web(i).xStart = sparCaps(2).posX - stringerGap*(i-1);
    web(i).xEnd = sparCaps(2).posX - stringerGap*(i);
    web(i).thickness = webThickness;
    web(i).zStart = get_z(web(i).xStart,0);
    web(i).zEnd = get_z(web(i).xEnd,0);
    dx = web(i).xStart-centroid.posX; dz = web(i).zStart-centroid.posZ;

    if i==1
        web(i).dp_area = sparCaps(2).area;
        web(i).dP_X = 0;
        web(i).dP_Z = 0;
        web(i).qPrime_X = 0;
        web(i).qPrime_Z = 0;
    else
        web(i).dp_area = tempStringers(i-1).area;
        web(i).dP_X = get_dp(dx,dz,Vx,0,Ix,Iz,Ixz,web(i).dp_area);
        web(i).dP_Z = get_dp(dx,dz,0,Vz,Ix,Iz,Ixz,web(i).dp_area);
        web(i).qPrime_X = web(i-1).qPrime_X - web(i).dP_X;
        web(i).qPrime_Z = web(i-1).qPrime_Z - web(i).dP_Z;
    end
    tempInt = get_int(web(i).xEnd,web(i).xStart,0);  %integral of airfoil function
    triangle1 = abs((web(i).xStart - sparCaps(2).posX)*web(i).zStart/2);
    triangle2 = abs((web(i).xEnd - sparCaps(2).posX)*web(i).zEnd/2);
    web(i).Area = tempInt + triangle1 - triangle2;
    web(i).ds = get_ds(web(i).xStart,web(i).xEnd,0);
    web(i).dS_over_t = web(i).ds / web(i).thickness;

    web(i).q_dS_over_t_X = web(i).qPrime_X * web(i).dS_over_t;
    web(i).q_dS_over_t_Z = web(i).qPrime_Z * web(i).dS_over_t;
    web(i).two_A_qprime_X = 2*web(i).Area*web(i).qPrime_X;
    web(i).two_A_qprime_Z = 2*web(i).Area*web(i).qPrime_Z;
    web(i).qp_dx_X = web(i).qPrime_X *(web(i).xEnd-web(i).xStart);
    web(i).qp_dx_Z = web(i).qPrime_Z *(web(i).xEnd-web(i).xStart);
    web(i).qp_dz_X = web(i).qPrime_X *(web(i).zEnd-web(i).zStart);
    web(i).qp_dz_Z = web(i).qPrime_Z *(web(i).zEnd-web(i).zStart);

    %web(i).radCurv = ...   Example:  get_curve(web(i).xStart,web(i).xEnd,1)
end
webLowerNose = web;
web = [];

%upper nose webs
numWebs = numNoseTopStringers;
stringerGap = upperNoseStringerGap;
webThickness = t_upper_front;
tempStringers = noseTopStringers;

for i=1:(numWebs+1)
    web(i).xStart = stringerGap*(i-1);
    web(i).xEnd = stringerGap*(i);
    web(i).thickness = webThickness;
    web(i).zStart = get_z(web(i).xStart,1);
    web(i).zEnd = get_z(web(i).xEnd,1);
    dx = web(i).xStart-centroid.posX; dz = web(i).zStart-centroid.posZ;
    if i==1
        web(i).dp_area = 0;
        web(i).dP_X = 0;
        web(i).dP_Z = 0;
        web(i).qPrime_X = webLowerNose(numNoseBottomStringers+1).qPrime_X - web(i).dP_X;
        web(i).qPrime_Z = webLowerNose(numNoseBottomStringers+1).qPrime_Z - web(i).dP_Z;
    else
        web(i).dp_area = tempStringers(i-1).area;
        web(i).dP_X = get_dp(dx,dz,Vx,0,Ix,Iz,Ixz,web(i).dp_area);
        web(i).dP_Z = get_dp(dx,dz,0,Vz,Ix,Iz,Ixz,web(i).dp_area);
        web(i).qPrime_X = web(i-1).qPrime_X - web(i).dP_X;
        web(i).qPrime_Z = web(i-1).qPrime_Z - web(i).dP_Z;
    end
    tempInt = get_int(web(i).xStart,web(i).xEnd,1);  %integral of airfoil function
    triangle2 = abs((web(i).xStart - sparCaps(2).posX)*web(i).zStart/2);
    triangle1 = abs((web(i).xEnd - sparCaps(2).posX)*web(i).zEnd/2);
    web(i).Area = tempInt + triangle1 - triangle2;
    web(i).ds = get_ds(web(i).xStart,web(i).xEnd,1);
    web(i).dS_over_t = web(i).ds / web(i).thickness;

    web(i).q_dS_over_t_X = web(i).qPrime_X * web(i).dS_over_t;
    web(i).q_dS_over_t_Z = web(i).qPrime_Z * web(i).dS_over_t;
    web(i).two_A_qprime_X = 2*web(i).Area*web(i).qPrime_X;
    web(i).two_A_qprime_Z = 2*web(i).Area*web(i).qPrime_Z;
    web(i).qp_dx_X = web(i).qPrime_X *(web(i).xEnd-web(i).xStart);
    web(i).qp_dx_Z = web(i).qPrime_Z *(web(i).xEnd-web(i).xStart);
    web(i).qp_dz_X = web(i).qPrime_X *(web(i).zEnd-web(i).zStart);
    web(i).qp_dz_Z = web(i).qPrime_Z *(web(i).zEnd-web(i).zStart);

    %web(i).radCurv = ...   Example:  get_curve(web(i).xStart,web(i).xEnd,1)
end
webUpperNose = web;
web = [];


%front Spar
i=1;
web(i).xStart = sparCaps(1).posX;
web(i).xEnd = sparCaps(2).posX;
web(i).thickness = t_frontSpar;
web(i).zStart = sparCaps(1).posZ;
web(i).zEnd = sparCaps(2).posZ;
web(i).dp_area = sparCaps(1).area;
dx = web(i).xStart-centroid.posX; dz = web(i).zStart-centroid.posZ;

web(i).dP_X = get_dp(dx,dz,Vx,0,Ix,Iz,Ixz,web(i).dp_area);
web(i).dP_Z = get_dp(dx,dz,0,Vz,Ix,Iz,Ixz,web(i).dp_area);
web(i).qPrime_X = webUpperNose(numNoseTopStringers+1).qPrime_X - web(i).dP_X;
web(i).qPrime_Z = webUpperNose(numNoseTopStringers+1).qPrime_Z - web(i).dP_Z;
web(i).Area = 0;
web(i).ds = abs(sparCaps(1).posZ - sparCaps(2).posZ);
web(i).dS_over_t = web(i).ds / web(i).thickness;
web(i).q_dS_over_t_X = web(i).qPrime_X * web(i).dS_over_t;
web(i).q_dS_over_t_Z = web(i).qPrime_Z * web(i).dS_over_t;
web(i).two_A_qprime_X = 2*web(i).Area*web(i).qPrime_X;
web(i).two_A_qprime_Z = 2*web(i).Area*web(i).qPrime_Z;
web(i).qp_dx_X = web(i).qPrime_X *(web(i).xEnd-web(i).xStart);
web(i).qp_dx_Z = web(i).qPrime_Z *(web(i).xEnd-web(i).xStart);
web(i).qp_dz_X = web(i).qPrime_X *(web(i).zEnd-web(i).zStart);
web(i).qp_dz_Z = web(i).qPrime_Z *(web(i).zEnd-web(i).zStart);

webFrontSparCell2 = web;
web = [];


%% Check that q'*dx sums up to Vx

Fx = sum([webTop.qp_dx_X])+webRearSpar.qp_dx_X+ sum([webBottom.qp_dx_X])+webFrontSpar.qp_dx_X;  %cell 1
Fx = Fx + sum([webLowerNose.qp_dx_X])+ sum([webUpperNose.qp_dx_X]);  %cell 2

Fz = sum([webTop.qp_dz_Z])+webRearSpar.qp_dz_Z+ sum([webBottom.qp_dz_Z])+webFrontSpar.qp_dz_Z;  %cell 1
Fz = Fz + sum([webLowerNose.qp_dz_Z])+ sum([webUpperNose.qp_dz_Z]);  %cell 2
%%

% sum up the ds/t and  q*ds/t to solve 2 equations, 2 unknowns

% [A]*[q1s q2s] = B

A11 = sum([webTop.dS_over_t])+webRearSpar.dS_over_t+ sum([webBottom.dS_over_t])+webFrontSpar.dS_over_t;
A22 = sum([webLowerNose.dS_over_t])+ sum([webUpperNose.dS_over_t])+webFrontSparCell2.dS_over_t;
A12 = -webFrontSpar.dS_over_t;
A21 = -webFrontSparCell2.dS_over_t;

B1_X = sum([webTop.q_dS_over_t_X])+webRearSpar.q_dS_over_t_X+ sum([webBottom.q_dS_over_t_X])+webFrontSpar.q_dS_over_t_X;
B2_X = sum([webLowerNose.q_dS_over_t_X])+ sum([webUpperNose.q_dS_over_t_X])+webFrontSparCell2.q_dS_over_t_X;
B1_Z = sum([webTop.q_dS_over_t_Z])+webRearSpar.q_dS_over_t_Z+ sum([webBottom.q_dS_over_t_Z])+webFrontSpar.q_dS_over_t_Z;
B2_Z = sum([webLowerNose.q_dS_over_t_Z])+ sum([webUpperNose.q_dS_over_t_Z])+webFrontSparCell2.q_dS_over_t_Z;

Amat = [A11 A12; A21 A22];
Bmat_X = -[B1_X;B2_X];
Bmat_Z = -[B1_Z;B2_Z];

qs_X = inv(Amat)*Bmat_X;
qs_Z = inv(Amat)*Bmat_Z;


sum_2_a_q_X = sum([webTop.two_A_qprime_X])+webRearSpar.two_A_qprime_X+ sum([webBottom.two_A_qprime_X]);  %cell 1 qprimes
sum_2_a_q_X = sum_2_a_q_X + sum([webLowerNose.two_A_qprime_X])+ sum([webUpperNose.two_A_qprime_X]);   %cell 2 qprimes
sum_2_a_q_X = sum_2_a_q_X +  qs_X(1)*(sum([webTop.Area])+webRearSpar.Area+ sum([webBottom.Area]));
sum_2_a_q_X = sum_2_a_q_X +  qs_X(2)*(sum([webLowerNose.Area])+ sum([webUpperNose.Area]));

sum_2_a_q_Z = sum([webTop.two_A_qprime_Z])+webRearSpar.two_A_qprime_Z+ sum([webBottom.two_A_qprime_Z]);  %cell 1 qprimes
sum_2_a_q_Z = sum_2_a_q_Z + sum([webLowerNose.two_A_qprime_Z])+ sum([webUpperNose.two_A_qprime_Z]);   %cell 2 qprimes
sum_2_a_q_Z = sum_2_a_q_Z +  qs_Z(1)*(sum([webTop.Area])+webRearSpar.Area+ sum([webBottom.Area]));
sum_2_a_q_Z = sum_2_a_q_Z +  qs_Z(2)*(sum([webLowerNose.Area])+ sum([webUpperNose.Area]));

%shear center
sc.posX =  sum_2_a_q_Z / Vz + frontSpar;
sc.posZ =  sum_2_a_q_X / Vx;


% now consider the torque representing shifting the load from the quarter
% chord to the SC  (need to check signs on these moments)
 
torque_Z = Vz*(sc.posX - 0.25);
torque_X = Vx*sc.posZ;
torque_Y = My;


Area1 = sum([webTop.Area]) + webRearSpar.Area + sum([webBottom.Area]);
%check area
Area1_check = get_int(frontSpar,backSpar,1) + get_int(frontSpar,backSpar,0);

Area2 = sum([webLowerNose.Area]) + sum([webUpperNose.Area]);
Area2_check = get_int(0,frontSpar,1) + get_int(0,frontSpar,0);


%for twist equation  (see excel spreadsheet example)

q1t_over_q2t = (A22/Area2 + webFrontSpar.dS_over_t/Area1)/(A11/Area1 + webFrontSpar.dS_over_t/Area2);

q2t = torque_X/(2*Area1*q1t_over_q2t + 2*Area2);
q1t = q2t*q1t_over_q2t;
qt_X = [q1t;q2t];

q2t = torque_Z/(2*Area1*q1t_over_q2t + 2*Area2);
q1t = q2t*q1t_over_q2t;
qt_Z = [q1t;q2t];

q2t = torque_Y/(2*Area1*q1t_over_q2t + 2*Area2);
q1t = q2t*q1t_over_q2t;
qt_Y = [q1t;q2t];


%% --- - add up all shear flows:  qtot = qPrime + qt + qs

webTop = sum_shear_flows(webTop, qt_X(1), qt_Y(1), qt_Z(1), qs_X(1), qs_Z(1));
webBottom = sum_shear_flows(webBottom, qt_X(1), qt_Y(1), qt_Z(1), qs_X(1), qs_Z(1));
webRearSpar = sum_shear_flows(webRearSpar, qt_X(1), qt_Y(1), qt_Z(1), qs_X(1), qs_Z(1));
webUpperNose = sum_shear_flows(webUpperNose, qt_X(2), qt_Y(2), qt_Z(2), qs_X(2), qs_Z(2));
webLowerNose = sum_shear_flows(webLowerNose, qt_X(2), qt_Y(2), qt_Z(2), qs_X(2), qs_Z(2));
webFrontSpar = sum_shear_flows(webFrontSpar, qt_X(1)-qt_X(2), qt_Y(1)-qt_Y(2), qt_Z(1)-qt_Z(2),...
                                qs_X(1)-qs_X(2), qs_Z(1)-qs_Z(2));


%% Shear flows and stresses, still need to superimpose tho
if loading_condition == 1
    unit_shear_flows_X = [webUpperNose.qtot_X webTop.qtot_X webRearSpar.qtot_X webBottom.qtot_X webFrontSpar.qtot_X webLowerNose.qtot_X];
    unit_shear_stresses_X = [[webUpperNose.qtot_X]/t_upper_front [webTop.qtot_X]/t_upper [webRearSpar.qtot_X]/t_rearSpar...
        [webBottom.qtot_X]/t_lower [webFrontSpar.qtot_X]/t_frontSpar [webLowerNose.qtot_X]/t_lower_front];
elseif loading_condition ==2
    unit_shear_flows_Z = [webUpperNose.qtot_Z webTop.qtot_Z webRearSpar.qtot_Z webBottom.qtot_Z webFrontSpar.qtot_Z webLowerNose.qtot_Z];
    unit_shear_stresses_Z = [[webUpperNose.qtot_Z]/t_upper_front [webTop.qtot_Z]/t_upper [webRearSpar.qtot_Z]/t_rearSpar...
        [webBottom.qtot_Z]/t_lower [webFrontSpar.qtot_Z]/t_frontSpar [webLowerNose.qtot_Z]/t_lower_front];
elseif loading_condition == 3
    unit_shear_flows_Y = [webUpperNose.qtot_Y webTop.qtot_Y webRearSpar.qtot_Y webBottom.qtot_Y webFrontSpar.qtot_Y webLowerNose.qtot_Y];
    unit_shear_stresses_Y = [[webUpperNose.qtot_Y]/t_upper_front [webTop.qtot_Y]/t_upper [webRearSpar.qtot_Y]/t_rearSpar...
        [webBottom.qtot_Y]/t_lower [webFrontSpar.qtot_Y]/t_frontSpar [webLowerNose.qtot_Y]/t_lower_front];
end


end

%% Shear Stress Distributions
numWebs = numTopStringers + numBottomStringers + numNoseTopStringers + numNoseBottomStringers + 6;
numCrossSections = length(PHAA_MX);
PHAA_shear_stress = zeros(numWebs, numCrossSections);
PLAA_shear_stress = zeros(numWebs, numCrossSections);
NHAA_shear_stress = zeros(numWebs, numCrossSections);
NLAA_shear_stress = zeros(numWebs, numCrossSections);

for i = 1:numCrossSections
   PHAA_shear_stress(:, i) = unit_shear_stresses_X'*PHAA_VX(i) + unit_shear_stresses_Z'*PHAA_VZ(i); 
   PLAA_shear_stress(:, i) = unit_shear_stresses_X'*PLAA_VX(i) + unit_shear_stresses_Z'*PLAA_VZ(i); 
   NHAA_shear_stress(:, i) = unit_shear_stresses_X'*NHAA_VX(i) + unit_shear_stresses_Z'*NHAA_VZ(i); 
   NLAA_shear_stress(:, i) = unit_shear_stresses_X'*NLAA_VX(i) + unit_shear_stresses_Z'*NLAA_VZ(i); 
end

%% Bending Stress Distribution
stringerLocationsX = [noseTopStringers.posX sparCaps(1).posX topStringers.posX sparCaps(3).posX sparCaps(4).posX...
                        bottomStringers.posX sparCaps(2).posX noseBottomStringers.posX];
stringerLocationsZ = [noseTopStringers.posZ sparCaps(1).posZ topStringers.posZ sparCaps(3).posZ sparCaps(4).posZ...
                        bottomStringers.posZ sparCaps(2).posZ noseBottomStringers.posZ];
                    
numStringers = length(stringerLocationsX);
PHAA_bending_stress = zeros(numStringers, numCrossSections);
PLAA_bending_stress = zeros(numStringers, numCrossSections);
NHAA_bending_stress = zeros(numStringers, numCrossSections);
NLAA_bending_stress = zeros(numStringers, numCrossSections);

for i = 1:numCrossSections
    frac1 = (Ix*PHAA_MZ(i) + Ixz*PHAA_MX(i))/(Ix*Iz - Ixz^2);
    frac2 = (Iz*PHAA_MX(i) + Ixz*PHAA_MZ(i))/(Ix*Iz - Ixz^2);
    PHAA_bending_stress(:, i) = (stringerLocationsX - centroid.posX)*frac1 - (stringerLocationsZ - centroid.posZ)*frac2;
    
    frac1 = (Ix*PLAA_MZ(i) + Ixz*PLAA_MX(i))/(Ix*Iz - Ixz^2);
    frac2 = (Iz*PLAA_MX(i) + Ixz*PLAA_MZ(i))/(Ix*Iz - Ixz^2);
    PLAA_bending_stress(:, i) = (stringerLocationsX - centroid.posX)*frac1 - (stringerLocationsZ - centroid.posZ)*frac2;
    
    frac1 = (Ix*NHAA_MZ(i) + Ixz*NHAA_MX(i))/(Ix*Iz - Ixz^2);
    frac2 = (Iz*NHAA_MX(i) + Ixz*NHAA_MZ(i))/(Ix*Iz - Ixz^2);
    NHAA_bending_stress(:, i) = (stringerLocationsX - centroid.posX)*frac1 - (stringerLocationsZ - centroid.posZ)*frac2;
    
    frac1 = (Ix*NLAA_MZ(i) + Ixz*NLAA_MX(i))/(Ix*Iz - Ixz^2);
    frac2 = (Iz*NLAA_MX(i) + Ixz*NLAA_MZ(i))/(Ix*Iz - Ixz^2);
    NLAA_bending_stress(:, i) = (stringerLocationsX - centroid.posX)*frac1 - (stringerLocationsZ - centroid.posZ)*frac2;
end

%% Getting F critical for bending
E = 1.044E+7;        % PSI (68.9 GPa)
k=1;
F_crit = (pi^2 * E * 0.0035) / (k * ribL)^2;
Sigma_crit = F_crit./[[noseTopStringers.area] [sparCaps(1).area] [topStringers.area] [sparCaps(3).area]...
                      [sparCaps(4).area] [bottomStringers.area] [sparCaps(2).area] [noseBottomStringers.area]];

% loading into table, FOS
stringerNames = cell(numStringers, 1);
for i = 1:numStringers
   if i <= numNoseTopStringers
       stringerNames(i) = {['topNoseStringer' int2str(i)]};
   elseif i == numNoseTopStringers + 1
       stringerNames(i) = {['frontSparTopCap']};
   elseif i <= numTopStringers + numNoseTopStringers + 1
       stringerNames(i) = {['topStringer' int2str(i - numNoseTopStringers - 1)]};
   elseif i == numTopStringers + numNoseTopStringers + 2
       stringerNames(i) = {['backSparTopCap']};
   elseif i == numTopStringers + numNoseTopStringers + 3
       stringerNames(i) = {['backSparBottomCap']};
   elseif i <= numTopStringers + numNoseTopStringers + numBottomStringers + 3
       stringerNames(i) = {['bottomStringer' int2str(i - numNoseTopStringers - numTopStringers - 3)]};
   elseif i == numTopStringers + numNoseTopStringers + numBottomStringers + 4
       stringerNames(i) = {['frontSparBottomCap']};
   else
       stringerNames(i) = {['bottomNoseStringer' int2str(i - numNoseTopStringers - numTopStringers - numBottomStringers - 4)]};
   end
end

%% Data Visualization
path = '../../Reports/CDR/plots/';

% FIXME (dividing by 10^6)
FOS = table(PHAA_bending_stress(:,205),abs(Sigma_crit'./PHAA_bending_stress(:,205)));
FOS.Properties.VariableNames = {'BendingStress', 'FOS'};
FOS.Properties.RowNames = stringerNames;
display(FOS)

%plotting airfoil cross-section

xChord = 0:.01:1;
upperSurface = zeros(1,length(xChord));
lowerSurface = zeros(1,length(xChord));

for i=1:length(xChord)
    upperSurface(i) = get_z(xChord(i),1);
    lowerSurface(i) = get_z(xChord(i),0);
end

if plotting

    figure(1); hold on; axis equal; grid on;
    %plot(xChord,z_camber,'-')
    plot(xChord,upperSurface,'-')
    plot(xChord,lowerSurface,'-')

    plot([sparCaps(1).posX sparCaps(2).posX],[sparCaps(1).posZ sparCaps(2).posZ],'-')
    plot([sparCaps(3).posX sparCaps(4).posX],[sparCaps(3).posZ sparCaps(4).posZ],'-')
    plot([sparCaps.posX],[sparCaps.posZ],'o')
    plot([topStringers.posX],[topStringers.posZ],'or')
    plot([bottomStringers.posX],[bottomStringers.posZ],'or')
    plot([noseTopStringers.posX],[noseTopStringers.posZ],'or')
    plot([noseBottomStringers.posX],[noseBottomStringers.posZ],'or')
    plot(centroid.posX/c,centroid.posZ/c,'rx')
    plot(sc.posX,sc.posZ,'gx')
    title('Cross Section Layout')
    xlabel('Chord Percent')
    ylabel('Chord Percent')
    print(strcat([path 'cross_section']), '-dpng');

    figure(2); hold on;
    scatter3(stringerLocationsX/c, stringerLocationsZ/c, PHAA_bending_stress(:,200),'filled')
    scatter3(xChord,upperSurface,zeros(1,length(upperSurface)));
    scatter3(xChord,lowerSurface,zeros(1,length(lowerSurface)));
    title('Bending Stress in each Stringer')
    xlabel('Chord Percent')
    ylabel('Chord Percent')
    zlabel('Bending Stress [psi]')
    view(-30,30)
    set(gcf, 'Position', [200, 200, 800, 400])
    print(strcat([path 'bending_stress_3d']), '-dpng');

    figure(3); hold on;
    for i = 1:length(PHAA_shear_stress(:,1));
        plot(PHAA_shear_stress(i,:))
    end
    title('Shear Stress Across the Wing')
    xlabel('X Distance [in]')
    ylabel('Shear Stress [psi]')
    print(strcat([path 'shear_stress_stringers']), '-dpng');

    figure(4); hold on;
    for i = 1:length(PHAA_bending_stress(:,1));
        plot(PHAA_bending_stress(i,:))
    end
    title('Bending Stress Across the Wing')
    xlabel('X Distance [in]')
    ylabel('Shear Stress [psi]')
    print(strcat([path 'bending_stress_stringers']), '-dpng');
end