function [VT,VF,t,f] = get_regular_blast_sequence(BlastTable,SiteX,SiteY,Vw,PPV)
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
[TSPar] = set_timeseries_parameters(Tmax,8/1000);
t = TSPar.t;
f = TSPar.f;
NUP = TSPar.NUP;

% Build blast sequecnce
VF = zeros(NUP,NBlast);
for k = 1:NBlast
    fo = BlastTable(k).fo;
    xi = BlastTable(k).xi;
    T = BlastTable(k).T;
    R = BlastTable(k).R;
    VF(:,k) = get_damped_armonic(f,T+R/Vw,fo,xi);
    [VTo,~] = Get_TS(VF(:,k),f);
    VF(:,k) = VF(:,k)/max(abs(VTo));
end
VF = sum(VF,2);
[VT,~] = Get_TS(VF,f);
FS = PPV/max(abs(VT));
VT = VT*FS;
VF = VF*FS;
end


function [TSPar] = set_timeseries_parameters(Tmax,delay)

dt = 10^floor(log10(min(delay)))/1000;
NP = ceil((1.2*Tmax)/dt);   


NFFT = pow2(nextpow2(NP)+1);
NUP = NFFT/2+1;
df = 1/(2*NUP*dt);
while df>0.1 % Si el df es muy chico hay ruido en la señal
    NP = 2*NP;
    NFFT = pow2(nextpow2(NP)+1);
    NUP = NFFT/2+1;
    df = 1/(2*NUP*dt);
end
f = 1/(2*dt)*linspace(0,1,NUP).';
TSPar.NFFT = NFFT;
TSPar.NUP = NUP;
TSPar.f = f;
TSPar.df = f(2)-f(1);
TSPar.t = linspace(0,(NUP-1)*dt,NUP);
TSPar.tmax = (NUP-1)*dt;
TSPar.dt = dt;
end
