clc
clear all
to = 0.3;
fo = 10;
c_gauss = (5-3*1i);
xi_arm = 0.25;

%% Timeseries prop
Fmax = 8*(fo+5*(1+abs(c_gauss)));
dfmin = 0.1;
dt = 1/Fmax;
NFFT = 2^ceil(log2(1/(dfmin*dt)));
df = 1/(NFFT*dt);
NUP = NFFT/2+1;
f = 1/(2*dt)*linspace(0,1,NUP).';
t = linspace(0,(NUP-1)*dt,NUP);

%% Gaussian pulse
[VFgauss] = get_gaussian_noise(f,to,fo,c_gauss);
[VTgauss,~] = Get_TS(VFgauss,f);
AFgauss = VFgauss*2*pi.*f;
[ATgauss,~] = Get_TS(AFgauss,f);
FSgauss = 1/max(abs(VTgauss));
VFgauss = VFgauss*FSgauss;
VTgauss = VTgauss*FSgauss;
AFgauss = AFgauss*FSgauss;
ATgauss = ATgauss*FSgauss;

%% Damped armonic pulse
Vpulse = @(t,to,fo,xi) sin(2*pi*fo*(t-to)).*exp(-2*pi*fo*xi*(t-to)).*(t>=to);
VTarm = Vpulse(t,to,fo,xi_arm).';
VTarm(isnan(VTarm)) = 0;
[VFarm,farm] = Get_FS(VTarm,t);
AFarm = VFarm*2*pi.*farm;
[ATarm,~] = Get_TS(AFarm,farm);
FSarm = 1/max(abs(VTarm));
VFarm = VFarm*FSarm;
VTarm = VTarm*FSarm;
AFarm = AFarm*FSarm;
ATarm = ATarm*FSarm;

%% Plots

idx1 = find(cumsum(VTgauss.^2)/(dot(VTgauss,VTgauss))>0.999,1);
idx2 = find(cumsum(VTarm.^2)/(dot(VTarm,VTarm))>0.999,1);
idx = max(idx1,idx2);

close all
hfig = figure(1);
set(hfig,'Color',[1 1 1],'Position',[30,30,1000,300]);
hold on
plot(f,abs(VFgauss),'-r','linewidth',2);
plot(farm,abs(VFarm),'-b','linewidth',2);
hold off
xlim([0 2*fo]);
xlabel('f');
ylabel('|VF|');
ylim([0 max(1.2*max(abs(VFgauss)),1.2*max(abs(VFarm)))])
grid on
legend({'Gaussian pulse','Armonic pulse'});
set(gca,'Position',[0.05,0.14,0.85,0.79]);

hfig = figure(2);
set(hfig,'Color',[1 1 1],'Position',[150,200,1000,300]);
hold on
plot(t(1:idx),VTgauss(1:idx),'-r','linewidth',2);
plot(t(1:idx),VTarm(1:idx),'-b','linewidth',1);
hold off
grid on
xlabel('t');
ylabel('VT');
legend({'Gaussian pulse','Armonic pulse'});
set(gca,'Position',[0.05,0.14,0.85,0.79]);

hfig = figure(3);
set(hfig,'Color',[1 1 1],'Position',[300,400,1000,300]);
hold on
plot(t(1:idx),ATgauss(1:idx),'-r','linewidth',2);
plot(t(1:idx),ATarm(1:idx),'-b','linewidth',1);
hold off
grid on
xlabel('t');
ylabel('AT');
legend({'Gaussian pulse','Armonic pulse'});
set(gca,'Position',[0.05,0.14,0.85,0.79]);