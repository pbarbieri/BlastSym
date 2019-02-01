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
    tr = round(linspace(0,(2*NPo-1)*dt,2*NPo).'-NPo*dt,abs(ceil(log10(dt))));
    for k = 1:NBlast
        vo = Vpulse(tr,to(k),fo(k),xi(k));
        v = zeros(NPo,1);
        for j = 1:NPo
            v(j) = flipud(vo(j+1:NPo+j)).'*r(:,k);
        end
        VT(:,k) = (1-R)*vo(NPo+1:2*NPo)+ R/sum(r(:,k).^2)*v;
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
        FS = 1./max(abs(VT(:,k)));
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
