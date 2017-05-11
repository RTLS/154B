% Define the Variables %
S=174;          %ft^2%
AR=7.44;
l=.687;
nmax=3.4;
nmin=-1.2;
CLpos=1.6;
CLneg=-.65;
p=.002377;      %slugs/ft^3%
W=2300;      %slug*ft/s^2%
 
% Define Equations %
V1=((nmax*2*W)/(CLpos*p*S))^(1/2)
V2=((nmin*2*W)/(CLneg*p*S))^(1/2)
Vpos=[0:V1];
Vneg=[0:V2];
q=50.8;
n1=(CLpos*p*S*Vpos.^2)/(2*W);
n2=(CLneg*p*S*Vneg.^2)/(2*W);
Vmax=((2*q)/p)^(1/2)
 
 
% Plotting Vn diagram %
plot(Vpos, n1)
hold on
plot(Vneg, n2)
plot([V1:.5:Vmax],nmax*ones(length([V1:.5:Vmax])))
plot([V2:.5:Vmax],nmin*ones(length([V2:.5:Vmax])))
plot(Vmax*ones(length([nmin:.01:nmax])),[nmin:.01:nmax])
 
 
% Labeling Vn Diagram %
% gtext ('Positive Stall Limit')
% gtext ('Negative Stall Limit')
% gtext ('Positive Structural Limit')
% gtext ('Negative Structural Limit')
xlabel ('Calibrated Airspeed, ft/sec')
ylabel ('Load Factor, n')
