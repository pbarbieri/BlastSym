function [VT,AT,t] = get_regular_blast_sequence(BlastSeqTable,SiteX,SiteY,Vw,PPV,R)
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
%       SiteX           Site X coordinate
%       SiteY           Site Y coordiante
%       Vw              Velocity of wave propagation
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
Yblast = BlastSeqTable.Y;
Xblast = BlastSeqTable.X;
T = BlastSeqTable.T;
fo = BlastSeqTable.fo;
xi = BlastSeqTable.xi;
D = sqrt((SiteX-Xblast).^2+(SiteY-Yblast).^2);
to = T+D/Vw;
BlastSeqTable.D = D;


% Set timeseries parameters
Tmax = (max(BlastSeqTable.T)+max(BlastSeqTable.D)/Vw)*1.2;
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

% Build blast sequecnce
VT = zeros(NPo,NBlast);
AT = zeros(NPo,NBlast);
% % Randomizer
% r = zeros(NBlast*NPo,1);
% r(1) = 0.52; for k = 2:NBlast*NPo, r(k) = 2*r(k-1)^2-1; end
% r = reshape(r,NPo,NBlast);
% tr = linspace(0,(2*NPo-1)*dt,2*NPo).'-NPo*dt;
for k = 1:NBlast
%     vo = Vpulse(tr,to(k),fo(k),xi(k));
%     v = zeros(NPo,NPo);
%     for j = 1:NPo
%         v(:,j) = vo(NPo-j+1:2*NPo-j).*r(:,k);
%     end
%     VT(:,k) = (1-R)*vo(NPo+1:2*NPo)+ R/sum(r(:,k))*sum(v,2);
    VT(:,k) = Vpulse(t,to(k),fo(k),xi(k));
    FS = 1/max(abs(VT(:,k)));
    VT(:,k) = FS*VT(:,k);
    AT(:,k) = FS*Apulse(t,to(k),fo(k),xi(k));
end
VT = sum(VT,2);
AT = sum(AT,2);
FS = PPV/max(abs(VT));
VT = VT*FS;
AT = AT*FS;

% Remove unncesary pre-pading
tini = min(to)*0.8;
VT = VT(t>tini);
AT = AT(t>tini);
t = t(t>tini); t = t-t(1);
end
