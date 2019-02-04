function [] = main_r2(XlsFile,fo,xi,sigmaDelay,Vw,PPV,SiteX,SiteY,R,OutputFileName)
%% ========================================================================
% Copyright SRK/FIUBA (C) 2019
% Coded By: P. Barbieri (pbarbieri@fi.uba.ar)
% Version: V103
% Date: 30-01-2019
% -------------------------------------------------------------------------
% USAGE: 
%   main: Simulates blast vibrations
%   
% -------------------------------------------------------------------------
% INPUT:
%       XlsFile         Excel file with detonation sequence
%       fo              Predominant frequency of indvidual blast [Hz]
%       xi              Damping of indvidual blast [Hz]   
%       sigmaDelay      Standar deviation of the delay [ms]
%       Vw              Velocity of wave propagation
%       PPV             PPV [m/s]
%       SiteX           Site X coordinate
%       SiteY           Site Y coordiante
% -------------------------------------------------------------------------
% OUTPUT:
%       
% -------------------------------------------------------------------------
% EXAMPLE:
%   
% -------------------------------------------------------------------------
% BIBLIO:
%   [1] D. Blair, "Blast vibration control in the presence of delay scatter
%       and random fluctuations between blastholes," International journal 
%       for numerical and analytical methods in geomechanics
%   [2] D. Blair, "Statistical models for ground vibration and airblast,"
%       Fragblast, vol. 3, no. 4, pp. 335--364, 1999
% -------------------------------------------------------------------------
% VALIDATE:
% Version:
% Date:
% Validated by:
% -------------------------------------------------------------------------
% LOG
%   V101    30/01/2019      First version
% =========================================================================

clc
% Build blast table
[BlastSeqTable] = build_blast_sequence_table(XlsFile,fo,xi,SiteX,SiteY,sigmaDelay,Vw);

% Build blast sequence
[UT,VT,AT,t] = get_regular_blast_sequence(BlastSeqTable,PPV,R);

save([OutputFileName,'.mat'],'BlastSeqTable','UT','VT','AT','t');
write_csv_output(AT,t,OutputFileName);

% % Plots
% close all
% hfig = figure(1);
% set(hfig,'Color',[1 1 1],'Position',[50,50,470*1.5,400*1.5]);
% hold on
% scatter(BlastSeqTable.X,BlastSeqTable.Y,'ok','filled');
% scatter([SiteX,SiteX],[SiteY,SiteY],'^r','filled');
% hold off
% grid on
% legend({'Blast Holes','Site'},'location','southwest','FontSize',12);
% xlabel('X [m]','FontSize',12);
% ylabel('Y [m]','FontSize',12);
% 
% hfig = figure(2);
% set(hfig,'Color',[1 1 1],'Position',[100,100,1000,300]);
% plot(t,UT)
% grid on
% xlabel('t [s]');
% ylabel('U [m]');
% set(gca,'Position',[0.07,0.14,0.85,0.8]);
% 
% hfig = figure(3);
% set(hfig,'Color',[1 1 1],'Position',[150,150,1000,300]);
% plot(t,VT)
% grid on
% xlabel('t [s]');
% ylabel('V [m/s]');
% set(gca,'Position',[0.07,0.14,0.85,0.8]);
% 
% hfig = figure(4);
% set(hfig,'Color',[1 1 1],'Position',[200,200,1000,300]);
% plot(t,AT)
% grid on
% xlabel('t [s]');
% ylabel('A [m/s/s]');
% set(gca,'Position',[0.07,0.14,0.85,0.8]);

end

function [BlastSeqTable] = build_blast_sequence_table(XlsFile,fo,xi,SiteX,SiteY,sigmaDelay,Vw)

% Read data
[~,~,BlastSeqTable] = xlsread(XlsFile);
% Build Blast sequence table
BlastSeqTable = cell2struct(BlastSeqTable(2:end,:),BlastSeqTable(1,:),2);
BlastSeqTable = struct2table(BlastSeqTable);
NBlast = size(BlastSeqTable,1);
Yblast = BlastSeqTable.Y;
Xblast = BlastSeqTable.X;
T = BlastSeqTable.T;
W = BlastSeqTable.W;

D = sqrt((SiteX-Xblast).^2+(SiteY-Yblast).^2);
R = D./sqrt(W);
fo = ones(NBlast,1)*fo;
xi = ones(NBlast,1)*xi;
RandomDelay = random(makedist('normal',0,sigmaDelay),NBlast,1);
RandomDelay(1) = 0; % Avoid T<0 in first blast
T = (T + RandomDelay)/1000;
to = T+D/Vw;

a  = -1.6863;
b = 6.2693;
muLnPPV = @(r)a*log(r)+b;
sigmaLnPPV = 0.4880;
epsilon = random(makedist('normal',0,1),1,1);
PPV = exp(muLnPPV(R)+epsilon*sigmaLnPPV);

BlastSeqTable.T = BlastSeqTable.T/1000;
BlastSeqTable.to = to;
BlastSeqTable.D = D;
BlastSeqTable.R = R;
BlastSeqTable.PPV = PPV;
BlastSeqTable.fo = fo;
BlastSeqTable.xi = xi;

end

function [UT,VT,AT,t] = get_regular_blast_sequence(BlastSeqTable,PPV,R)
%% ========================================================================
% Copyright SRK/FIUBA (C) 2018
% Coded By: P. Barbieri (pbarbieri@fi.uba.ar)
% Version: V103
% Date: 19-01-2019
% -------------------------------------------------------------------------
% USAGE: 
%   blast_vibration_simulation: Simulates blast vibrations
%   
% -------------------------------------------------------------------------
% INPUT:
%       BlastTable      Table with the blast sequence with fields:
%                       - X & Y [m] coordinates of the blast
%                       - T [s] time of explotion
%                       - fo [Hz] natural frec of each blast
%                       - xi []   dampling coefficient of each blast
%       PPV             PPV [m/s]
%       R               % of randomization
% -------------------------------------------------------------------------
% OUTPUT:
%       
% -------------------------------------------------------------------------
% EXAMPLE:
%   
% -------------------------------------------------------------------------
% BIBLIO:
%   [1] D. Blair, "Blast vibration control in the presence of delay scatter
%       and random fluctuations between blastholes," International journal 
%       for numerical and analytical methods in geomechanics
%   [2] D. Blair, "Statistical models for ground vibration and airblast,"
%       Fragblast, vol. 3, no. 4, pp. 335--364, 1999
% -------------------------------------------------------------------------
% VALIDATE:
% Version:
% Date:
% Validated by:
% -------------------------------------------------------------------------
% LOG
%   V101    19/01/2019      First version
% =========================================================================

NBlast = size(BlastSeqTable,1);
to = BlastSeqTable.to;
fo = BlastSeqTable.fo;
xi = BlastSeqTable.xi;
PVV = BlastSeqTable.PPV;


% Set timeseries parameters
Tmax = max(BlastSeqTable.to)*1.2;
dt = 10^(floor(log10(1/(8*min(fo(fo>0))))));
NPo = ceil(Tmax/dt)+1;
t = linspace(0,(NPo-1)*dt,NPo).';
Tpulse = log(0.001)/(-2*pi*min(fo)*min(xi));
NPpulse = ceil(Tpulse/dt);
tpulse = linspace(0,(NPpulse-1)*dt,NPpulse).';

% Single blast generatorn
Delta = 0.01; % remaining amplitud at -t1
t1 = 10*dt;
k = log(1/Delta-1)/t1;
L = @(t) 1./(1+exp(-k*t));
Upulse = @(t,to,fo,xi) sin(2*pi*fo*(t-to-t1)).*exp(-xi*2*pi*fo*(t-to-t1)).*L(t-to-t1);

Vpulse = @(t,to,fo,xi) (2*pi*fo*cos(2*pi*fo*(t-to-t1)).*min(exp(-xi*2*pi*fo*(t-to-t1)),1000).*L(t-to-t1)...
                       - xi*2*pi*fo*sin(2*pi*fo*(t-to-t1)).*min(exp(-xi*2*pi*fo*(t-to-t1)),1000).*L(t-to-t1)...
                       + sin(2*pi*fo*(t-to-t1)).*min(exp(-xi*2*pi*fo*(t-to-t1)),1000).*k.*L(t-to-t1).*(1-L(t-to-t1)));

Apulse = @(t,to,fo,xi) -(2*pi*fo).^2*sin(2*pi*fo*(t-to-t1)).*exp(-xi*2*pi*fo*(t-to-t1)).*L(t-to-t1)...
                       -(2*pi*fo).^2*xi*cos(2*pi*fo*(t-to-t1)).*exp(-xi*2*pi*fo*(t-to-t1)).*L(t-to-t1)...
                       + 2*pi*fo*cos(2*pi*fo*(t-to-t1)).*exp(-xi*2*pi*fo*(t-to-t1)).*k.*L(t-to-t1).*(1-L(t-to-t1))...
                       - xi*(2*pi*fo)^2*cos(2*pi*fo*(t-to-t1)).*exp(-xi*2*pi*fo*(t-to-t1)).*L(t-to-t1)...
                       + (xi*2*pi*fo)^2*sin(2*pi*fo*(t-to-t1)).*exp(-xi*2*pi*fo*(t-to-t1)).*L(t-to-t1)...
                       - xi*2*pi*fo* sin(2*pi*fo*(t-to-t1)).*exp(-xi*2*pi*fo*(t-to-t1)).*k.*L(t-to-t1).*(1-L(t-to-t1))... 
                       + 2*pi*fo*cos(2*pi*fo*(t-to-t1)).*exp(-xi*2*pi*fo*(t-to-t1)).*k.*L(t-to-t1).*(1-L(t-to-t1))...
                       - 2*pi*fo*xi*sin(2*pi*fo*(t-to-t1)).*exp(-xi*2*pi*fo*(t-to-t1)).*k.*L(t-to-t1).*(1-L(t-to-t1))...
                       + sin(2*pi*fo*(t-to-t1)).*exp(-xi*2*pi*fo*(t-to-t1)).*k.^2.*(L(t-to-t1).*(1-L(t-to-t1)).^2-L(t-to-t1).^2.*(1-L(t-to-t1)));        

% Build blast sequecnce
VT = zeros(NPo,NBlast);
if R>0
    % Randomizer
    r = zeros(NBlast*NPo,1);
    r(1) = 0.52; for k = 2:NBlast*NPo, r(k) = 2*r(k-1)^2-1; end
    r = reshape(r,NPo,NBlast);
    for k = 1:NBlast
        vo = Vpulse(t,to(k),fo(k),xi(k));
        v1 = conv(r(:,k),vo);
        VT(:,k) = (1-R)*vo+ R/sum(r(:,k).^2)*v1(1:NPo);
        VT(:,k) = PVV(k)*VT(:,k)/max(abs(VT(:,k)));
    end
    VT = sum(VT,2);
    UT = zeros(NPo,1);
    AT = zeros(NPo,1);
    for k = 2:NPo-1
        AT(k) = (VT(k+1)-VT(k-1))/(2*dt);
        UT(k) = UT(k-1)+(1-1/2)*dt*VT(k-1)+1/2*dt*VT(k);
    end
    AT(NPo) = (VT(end)-VT(end-1))/dt;
    UT(NPo) = UT(NPo-1)+(1-1/2)*dt*VT(NPo-1)+1/2*dt*VT(NPo);
else
    UT = zeros(NPo,NBlast);
    AT = zeros(NPo,NBlast);
    for k = 1:NBlast
        VT(:,k) = Vpulse(t,to(k),fo(k),xi(k));
        FS = PPV(k)./max(abs(VT(:,k)));
        VT(:,k) = FS*VT(:,k);
        UT(:,k) = FS*Upulse(t,to(k),fo(k),xi(k));
        AT(:,k) = FS*Apulse(t,to(k),fo(k),xi(k));
    end 
    UT = sum(UT,2);
    VT = sum(VT,2);
    AT = sum(AT,2);
end
FS = PPV/max(abs(VT));
VT = VT*FS;
AT = AT*FS;
UT = UT*FS;

% Remove unncesary pre-pading
tini = min(to)*0.8;
UT = UT(t>tini);
VT = VT(t>tini);
AT = AT(t>tini);
t = t(t>tini); t = t-t(1);
end

function [] = write_csv_output(AT,t,OutputFileName)
idx = find(round(AT,5)~=0,1);
if idx>1
    AT = AT(idx-1:end);
    t = t(1:numel(AT));
end
NP = size(AT,1);
AT(1) = 0;
AT = AT/9.81;
AT(round(AT,6)==0) = 0;

FileName = [OutputFileName,'.txt'];
fiad = fopen(FileName,'w');
for k = 1:NP
    fprintf(fiad,'%f\t%f\n',t(k),AT(k));
end
fclose(fiad);

end



