function [VT,AT,VF,AF,t,f] = get_regular_blast_sequence(BlastTable,SiteX,SiteY,Vw,PPV)
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


pulse = @(t,to,fo,xi) sin(2*pi*fo*(t-to)).*exp(-xi*2*pi*fo*(t-to)).*(t>=to);


NBlast = size(BlastTable,1);
BlastTable = table2struct(BlastTable);
for k = 1:NBlast
    Xblast = BlastTable(k).X;
    Yblast = BlastTable(k).Y;
    R = sqrt((SiteX-Xblast)^2+(SiteY-Yblast)^2);
    BlastTable(k).R=R;
end

% Set timeseries parameters
Tmax = max([BlastTable.T])+max([BlastTable.R])/Vw*1.2;
fo = min([BlastTable.fo]);
dt = 10^(floor(log10(1/(4*fo))));
NPo = ceil(Tmax/dt)+1;
t = linspace(0,(NPo-1)*dt,NPo);
% Build blast sequecnce
VT = zeros(NPo,NBlast);
for k = 1:NBlast
    fo = BlastTable(k).fo;
    xi = BlastTable(k).xi;
    T = BlastTable(k).T;
    R = BlastTable(k).R;
    to = T+R/Vw;
    VT(:,k) = pulse(t,to,fo,xi);
    VT(:,k) = VT(:,k)/max(abs(VT(:,k)));
end
VT = sum(VT,2);
VT = VT/max(abs(VT))*PPV;
VT(isnan(VT)) = 0;

NFFT = pow2(nextpow2(NPo)+1);
df = 1/(NFFT*dt);
while df>0.1
    NFFT = NFFT*2;
    df = 1/(NFFT*dt);
end
NUP = NFFT/2+1;
VT(NUP) = 0;
t = linspace(0,(NUP-1)*dt,NUP);

[VF,f] = Get_FS(VT,t);
AF = 2*pi*f.*VF;
[AT,t] = Get_TS(AF,f);

t = t(1:NPo);
VT = VT(1:NPo);
AT = AT(1:NPo);
end
