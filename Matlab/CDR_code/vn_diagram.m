%Script for constructing V-n diagram.  Modify with better estimates for
%aero terms

close all;

W = 3200;
S = 5*17*2;
rho = 0.00238;
Vc = 230;
Vd = 270;

Cz_max = 1.5;   %rough estimate for now
Cz_neg = -1.0;  %rough estimate for now
kg = 0.7;       %rough estimate for now
Cza = 3.7;        %rough estimate for now


v = 0:270; %mph
vfps = v*5280/3600; %fps
n_lift = (0.5*rho*S*Cz_max/W).*vfps.^2;
n_neg =  (0.5*rho*S*Cz_neg/W).*vfps.^2;

lim_load_plus = 4.4;
lim_load_minus = -1.76;

g_50_plus = 1+kg*rho*50*Cza*S/(2*W).*vfps;
g_30_plus = 1+kg*rho*30*Cza*S/(2*W).*vfps;
g_50_minus = 1-kg*rho*50*Cza*S/(2*W).*vfps;
g_30_minus = 1-kg*rho*30*Cza*S/(2*W).*vfps;


figure; grid on; hold on;set(gcf,'color',[1 1 1]);

tempInd = find(n_lift<lim_load_plus); plot(v(tempInd),n_lift(tempInd),'linewidth',2); 
tempInd = find(n_neg>lim_load_minus); plot(v(tempInd),n_neg(tempInd),'linewidth',2);

tempInd = find(n_lift>lim_load_plus,1); plot([v(tempInd) v(end)],[lim_load_plus lim_load_plus],'linewidth',2);

tempInd = find(n_neg<lim_load_minus,1); plot([v(tempInd) Vc],[lim_load_minus lim_load_minus],'linewidth',2);

plot([Vd Vd],[-1 lim_load_plus],'linewidth',2)
plot([Vc Vd],[lim_load_minus -1],'linewidth',2)
plot(v(1:230),g_50_plus(1:230),'--g','linewidth',2)
plot(v(1:270),g_30_plus(1:270),'--g','linewidth',2)
plot(v(1:230),g_50_minus(1:230),'--g','linewidth',2)
plot(v(1:270),g_30_minus(1:270),'--g','linewidth',2)
plot([Vc Vd],[g_50_plus(230) g_30_plus(270)],'--g','linewidth',2)
plot([Vc Vd],[g_50_minus(230) g_30_minus(270)],'--g','linewidth',2)



plot([230 230],[lim_load_minus lim_load_plus],'--r','linewidth',1)
xlabel('V (mph)','fontsize',16,'fontweight','bold');ylabel('n','fontsize',16,'fontweight','bold')
set(gca,'FontSize',16,'fontweight','bold');
ylim([-3 6])
title('V-N Diagram')