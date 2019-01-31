
clc
clear all
Ka = 0.2; % Acceleration amplification factor
f = 30; % Expected frequency
K = 2094/1000; % Parameter of attenuation law
N = 1.686;  % Parameter of attenuation law

W = @(D,ky) (D./(((ky*9.81)./(2*pi*f*Ka*K)).^(1/-N))).^2;

D = linspace(0,300,301).';

close all
hfig = figure(1);
set(hfig,'color',[1 1 1],'Position',[40 40 800 500]);
hplot = gobjects(4,1);
hold on
hplot(1) = plot(D,W(D,0.12),'-b','linewidth',2);
hplot(2) = plot(D,W(D,0.019),'-r','linewidth',2);
hplot(3) = plot(D,W(D,0.43),'-k','linewidth',2);
hplot(4) = plot(D,W(D,0.45),'-k','linewidth',2);
hold off
grid on
xlabel('Distance [m]','fontsize',12);
ylabel('Maximum charge weight [kg]','fontsize',12);
legend(hplot,{'Previous study(ky=0.12 K 95%)','Mining road(ky=0.019 K 95%)'...
             ,'Raise (ky=0.43 K 95%)','Global wedge(ky=0.45 K 95%)'}...
    ,'location','northwest','fontsize',12);
