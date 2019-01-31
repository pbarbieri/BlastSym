
% Set exitation parameters
Tmax = 1;
fo = 200;
xi = 0.25;
to = 0.3;

% Set timeseries parameters
dt = 10^(floor(log10(1/(8*min(fo(fo>0))))));
NPo = ceil(Tmax/dt)+1;
t = linspace(0,(NPo-1)*dt,NPo).';

% Single blast generatorn
Delta = 0.01; % remaining amplitud at -t1
t1 = 10*dt;
k = log(1/Delta-1)/t1;
L = @(t) 1./(1+exp(-k*t));
Vpulse = @(t,to,fo,xi) sin(2*pi*fo*(t-to)).*exp(-xi*2*pi*fo*(t-to)).*L(t-to);
Apulse = @(t,to,fo,xi) 2*pi*fo*cos(2*pi*fo*(t-to)).*exp(-xi*2*pi*fo*(t-to)).*L(t-to)...
                       - xi*2*pi*fo* sin(2*pi*fo*(t-to)).*exp(-xi*2*pi*fo*(t-to)).*L(t-to)...
                       + sin(2*pi*fo*(t-to)).*exp(-xi*2*pi*fo*(t-to)).*L(t-to).*(1-L(t-to));
                   
                   
VT = Vpulse(t,to,fo,xi);
AT_teo = Apulse(t,to,fo,xi);
AT_time = zeros(NPo,1);
for k = 2:NPo-1
    AT_time(k) = (VT(k+1)-VT(k-1))/(2*dt);
end
AT_time(NPo) = (VT(end)-VT(end-1))/dt;

[VF,f] = Get_FS(VT,t);
[AT_f,~] = Get_TS(2*pi*f.*VF,f);


close all
hfig = figure(1);
set(hfig,'Color',[1 1 1],'Position',[50,50,1000,300]);
hold on
plot(t(1:NPo),VT(1:NPo),'-b','linewidth',1);
hold off
grid on
xlabel('t');
ylabel('VT');
set(gca,'Position',[0.07,0.14,0.85,0.83]);


hfig = figure(2);
set(hfig,'Color',[1 1 1],'Position',[500,200,1000,300]);
hold on
plot(t(1:NPo),AT_teo(1:NPo),'-b','linewidth',1);
plot(t(1:NPo),AT_f(1:NPo),'-r','linewidth',1);
plot(t(1:NPo),AT_time(1:NPo),'-k','linewidth',1);
hold off
grid on
xlabel('t');
ylabel('AT');
legend({'Theorical integration','Freq. integration','Time integration'});
set(gca,'Position',[0.07,0.14,0.85,0.83]);


                   
                   