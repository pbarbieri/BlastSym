
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
Upulse = @(t,to,fo,xi) sin(2*pi*fo*(t-to)).*exp(-xi*2*pi*fo*(t-to)).*L(t-to);

Vpulse = @(t,to,fo,xi) 2*pi*fo*cos(2*pi*fo*(t-to)).*exp(-xi*2*pi*fo*(t-to)).*L(t-to)...
                       - xi*2*pi*fo*sin(2*pi*fo*(t-to)).*exp(-xi*2*pi*fo*(t-to)).*L(t-to)...
                       + sin(2*pi*fo*(t-to)).*exp(-xi*2*pi*fo*(t-to)).*k.*L(t-to).*(1-L(t-to));

Apulse = @(t,to,fo,xi) -(2*pi*fo).^2*sin(2*pi*fo*(t-to)).*exp(-xi*2*pi*fo*(t-to)).*L(t-to)...
                       -(2*pi*fo).^2*xi*cos(2*pi*fo*(t-to)).*exp(-xi*2*pi*fo*(t-to)).*L(t-to)...
                       + 2*pi*fo*cos(2*pi*fo*(t-to)).*exp(-xi*2*pi*fo*(t-to)).*k.*L(t-to).*(1-L(t-to))...
                       - xi*(2*pi*fo)^2*cos(2*pi*fo*(t-to)).*exp(-xi*2*pi*fo*(t-to)).*L(t-to)...
                       + (xi*2*pi*fo)^2*sin(2*pi*fo*(t-to)).*exp(-xi*2*pi*fo*(t-to)).*L(t-to)...
                       - xi*2*pi*fo* sin(2*pi*fo*(t-to)).*exp(-xi*2*pi*fo*(t-to)).*k.*L(t-to).*(1-L(t-to))... 
                       + 2*pi*fo*cos(2*pi*fo*(t-to)).*exp(-xi*2*pi*fo*(t-to)).*k.*L(t-to).*(1-L(t-to))...
                       - 2*pi*fo*xi*sin(2*pi*fo*(t-to)).*exp(-xi*2*pi*fo*(t-to)).*k.*L(t-to).*(1-L(t-to))...
                       + sin(2*pi*fo*(t-to)).*exp(-xi*2*pi*fo*(t-to)).*k.^2.*(L(t-to).*(1-L(t-to)).^2-L(t-to).^2.*(1-L(t-to)));                

AT = Apulse(t,to,fo,xi);
UT= Upulse(t,to,fo,xi);
VT_teo = Vpulse(t,to,fo,xi);

VT_time = zeros(NPo,1);
for k = 2:NPo-1
    VT_time(k) = (UT(k+1)-UT(k-1))/(2*dt);
end
VT_time(NPo) = (UT(end)-UT(end-1))/dt;


AT_time = zeros(NPo,1);
for k = 2:NPo-1
    AT_time(k) = (VT_time(k+1)-VT_time(k-1))/(2*dt);
end
AT_time(NPo) = (VT_time(end)-VT_time(end-1))/dt;


[UF,f] = Get_FS(UT,t);
[VT_f,~] = Get_TS(2*pi*f*1i.*UF,f);
[AT_f,~] = Get_TS((2*pi*f*1i).^2.*UF,f);

[VT_int,UT_int] = Get_VUT(AT_time,t);

close all
hfig = figure(1);
set(hfig,'Color',[1 1 1],'Position',[50,50,1000,300]);
hold on
plot(t(1:NPo),UT(1:NPo),'-b','linewidth',1);
plot(t(1:NPo),UT_int(1:NPo),'-r','linewidth',1);
hold off
grid on
xlabel('t');
ylabel('UT');
legend({'Seed','differentiatied and integrated'});
set(gca,'Position',[0.07,0.14,0.85,0.83]);

hfig = figure(2);
set(hfig,'Color',[1 1 1],'Position',[250,200,1000,300]);
hold on
plot(t(1:NPo),VT_teo(1:NPo),'-b','linewidth',1);
plot(t(1:NPo),VT_f(1:NPo),'-r','linewidth',1);
plot(t(1:NPo),VT_time(1:NPo),'-k','linewidth',1);
hold off
grid on
xlabel('t');
ylabel('VT');
legend({'Analitical differentiation','Freq. differentiation','Time differentiation'});
set(gca,'Position',[0.07,0.14,0.85,0.83]);

hfig = figure(3);
set(hfig,'Color',[1 1 1],'Position',[500,400,1000,300]);
hold on
plot(t(1:NPo),AT(1:NPo),'-b','linewidth',1);
plot(t(1:NPo),AT_f(1:NPo),'-r','linewidth',1);
plot(t(1:NPo),AT_time(1:NPo),'-k','linewidth',1);
hold off
grid on
xlabel('t');
ylabel('AT');
legend({'Analitical differentiation','Freq. differentiation','Time differentiation'});
set(gca,'Position',[0.07,0.14,0.85,0.83]);



                   
                   