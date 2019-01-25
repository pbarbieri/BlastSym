function [BlastTable] = BlastSimulation(BlastTable,BlastAttModel,SiteX,SiteY,Nsim,OutputFolder)
%% ========================================================================
% Copyright SRK/FIUBA (C) 2018
% Coded By: P. Barbieri (pbarbieri@fi.uba.ar)
% Version: V103
% Date: 19-01-2019
% -------------------------------------------------------------------------
% USAGE: 
%   BlastSimulation: Simulates blast vibrations
%   
% -------------------------------------------------------------------------
% INPUT:
%       BlastTable      Table with the blast sequence with fields:
%                       - X & Y [m] coordinates of the blast
%                       - W [kg] explosive mass
%                       - T [ms] time of explotion
%                       - S [ ]  secrrening flag (1 to consider screening)
%                       - D [m] depth of blast
%       BlastAttModel   Blast attenuation model defined by muLnPPV(R) and
%                       sigmalLnPPV(R) in [mm/s2] where R [m/kg^-05] is the
%                       adimetional distance.
%       SiteX           Site X coordinate
%       SiteY           Site Y coordiante
%       Nsim            Number of scenario simulations
%       OutputFolder    Folder where the simulated records will be saved
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

if ~isdir(OutputFolder), mkdir(OutputFolder); end
% Determine screening, muLnPPV, sigmaLnPPV, direction cosines of each blast
BlastTable = screening(BlastTable,SiteX,SiteY,BlastAttModel);

% Simulation of site properties
[SiteTable] = site_simulation(Nsim);

% Simulation of wave existation for each site
[BlastTable] = wave_simulation(BlastTable,SiteTable,Nsim,OutputFolder);



end

function [BlastTable,NBlast] = screening(BlastTable,SiteX,SiteY,AttenuationModel)

BlastTable = table2struct(BlastTable);
muLnPPV = AttenuationModel.muLnPPV;
sigmaLnPPV = AttenuationModel.sigmaLnPPV;

NBlast = size(BlastTable,1);
X = [BlastTable.X];
Y = [BlastTable.Y];
T = [BlastTable.T];
W = [BlastTable.W];

% Screening widht
D = zeros(NBlast,NBlast);
for k = 1:NBlast
    for j = 1:NBlast
        D(k,j) = sqrt((X(k)-X(j))^2+(Y(k)-Y(j))^2);
    end
end
D = triu(D);
D = D(:);
D = D(D>0);
ScreeningWidth = mode(D);

% Number of blast
for k = 1:NBlast
        
    % Remove future blasts
    X_filtered = X; X_filtered = X_filtered(T<T(k));
    Y_filtered = Y; Y_filtered = Y_filtered(T<T(k));
    % Remove blasts from same borehole
    idx = and(X_filtered~=X(k),Y_filtered~=Y(k));
    X_filtered = X_filtered(idx);
    Y_filtered = Y_filtered(idx);
    % unitary direction vector from site to blast
    l = [SiteX-X(k); SiteY-Y(k)]; % Versor blast 2 site
    l = l/dot(l,l);
    n = [-l(2),l(1)];    
    % directing cosines for radial vibration
    alpha = atan(l(2)/l(1));
    BlastTable(k).alpha = alpha;
    BlastTable(k).Rx = cos(alpha);
    BlastTable(k).Ry = sin(alpha);
    % directing cosines for tangential displacement
    BlastTable(k).Tx = cos(alpha+pi/2);
    BlastTable(k).Ty = sin(alpha+pi/2);
    
    % Screening polygon
    Coord = zeros(4,2);
    Coord(1,:) = [SiteX,SiteY]-n*ScreeningWidth;
    Coord(2,:) = [X(k),Y(k)]-n*ScreeningWidth;
    Coord(3,:) = [X(k),Y(k)]+n*ScreeningWidth;
    Coord(4,:) = [SiteX,SiteY]+n*ScreeningWidth;
    [in,on] = inpolygon(X_filtered,Y_filtered,Coord(:,1),Coord(:,2));
    idx = or(in,on);
    Ns = sum(idx);
    Dist = sqrt((SiteX-X(k))^2+(SiteY-Y(k))^2);
    BlastTable(k).S = 1/(1+23337*sqrt(Ns)/Dist^2);
    BlastTable(k).Ns = Ns;
    BlastTable(k).Dist = Dist;
    
    % Blast attenuation parameters
    R = Dist./sqrt(W);
    BlastTable(k).R = R;
    BlastTable(k).muLnPPV = muLnPPV(R);
    BlastTable(k).sigmaLnPPV = sigmaLnPPV(R);
end

BlastTable = struct2table(BlastTable);
end

function [TSPar] = set_timeseries_parameters(BlastTable)

T = [BlastTable.T]/1000; % ms 2 s
T = T(T>0);
if ~isempty(T)
    dt = 10^floor(log10(min(T)))/1000;
    NP = ceil(1.2*max(T)/dt);   
else
    dt = 1/1000^2;
    NP = ceil(20/1000/dt); % 20 ms
end

NFFT = pow2(nextpow2(NP)+1);
NUP = NFFT/2+1;
f = 1/(2*dt)*linspace(0,1,NUP).';
TSPar.NFFT = NFFT;
TSPar.NUP = NUP;
TSPar.f = f;
TSPar.df = f(2)-f(1);
TSPar.t = linspace(0,(NUP-1)*dt,NUP);
TSPar.tmax = (NUP-1)*dt;

end

function [SiteTable] = site_simulation(Nsim)

% Pressure wave velocity  [m/s]
Vpmin = 1800;
Vppeak = 2000;
Vpmax = 2200;
VpDist = makedist('triangular','a',Vpmin,'b',Vppeak,'c',Vpmax);
SiteTable.Vp = random(VpDist,Nsim,1);

% Shear wave velocity [m/s]
Vsmin = 800;
Vspeak = 1000;
Vsmax = 1200;
VsDist = makedist('triangular','a',Vsmin,'c',Vsmax,'b',Vspeak);
Vs = random(VsDist,Nsim,1);
SiteTable.Vs = Vs;

% Poisson coefficient  [ ]
numin = 0.18;
mumax = 0.22;
nupeak = 0.20;
nuDist = makedist('triangular','a',numin,'c',mumax,'b',nupeak);
nu = random(nuDist,Nsim,1);
SiteTable.nu = nu;

% Raleigh wave velocity [m/s]
SiteTable.Vr = (0.054.*nu.^4-0.08*nu.^3-0.038*nu.^2+0.195*nu+0.874).*Vs;
SiteTable = struct2table(SiteTable);

% Energy divition between waves
CompWaveEnergyMin = 0.05;
CompWaveEnergyMax = 0.15;
CompWaveEnergy = random(makedist('uniform',CompWaveEnergyMin,CompWaveEnergyMax),Nsim,1);
ShearWaveEnergyMin = 0.2;
ShearWaveEnergyMax = 0.3;
ShearWaveEnergy = random(makedist('uniform',ShearWaveEnergyMin,ShearWaveEnergyMax),Nsim,1);
SurfaceWaveEnergy = 1 - CompWaveEnergy - ShearWaveEnergy;

SiteTable.CompWaveEnergy = CompWaveEnergy;
SiteTable.ShearWaveEnergy = ShearWaveEnergy;
SiteTable.SurfaceWaveEnergy = SurfaceWaveEnergy;
end

function [SiteTable] = wave_simulation(BlastTable,SiteTable,Nsim,OutputFolder)
% Set timeseries parameters
[TSPar] = set_timeseries_parameters(BlastTable);

NBlast = size(BlastTable,1);

% Attenuaiton
muLnPPV = [BlastData.muLnPPV];
sigmaLnPPV = [BlastData.sigmaLnPPV];
NsigmaPPV = random(makedist('normal',0,1),Nsim,1);
PPV = zeros(Nsim,Nblast);
for k = 1:Nsim
    for j = 1:Nblast
        PPV(k,j) = exp(muLnPPV(j)+NsigmaPPV(k)*sigmaLnPPV(j));
    end 
end

% Chaotic number generator
r = zeros(Nsim,NBlast);
chaos = 0.54;
for j = 1:Nsim
    for k = 1:Nblast
        r(j,k) = chaos;
        chaos = 2*chaos^2-1;
    end
end

% Delay simulation
muDelay = 0; %[s]
sigmaDelay = 0.1/1000; %[s]
Delay = random(makedist('normal',muDelay,sigmaDelay),Nsim,NBlast) + repmat(BlastData.T.',Nsim,1);

% Arrival P wave time
PArrivalTime = repmat(BlastData.Dist.',Nsim,1)./repmat(SiteTable.Vp,1,NBlast) + Delay;
% Arrival S wave time
SArrivalTime = repmat(BlastData.Dist.',Nsim,1)./repmat(SiteTable.Vs,1,NBlast) + Delay;
% Arrival R wave time
RArrivalTime = repmat(BlastData.Dist.',Nsim,1)./repmat(SiteTable.Vr,1,NBlast) + Delay;

% Single blast P-wave frequency range [Hz] 
PWFmin = 100;
PWFmax = 200;
PWFreqDist = makedist('uniform',PWFmin,PWFmax);
PWFreq = random(PWFreqDist,Nsim,1);


% Single blast S-wave frequency range [Hz] 
SWFmin = 40;
SWFmax = 60;
SWFreqDist = makedist('uniform',SWFmin,SWFmax);
SWFreq = random(SWFreqDist,Nsim,1);

% Single blast R-wave frequency range [Hz] 
RWFmin = 10;
RWFmax = 20;
RWFreqDist = makedist('uniform',RWFmin,RWFmax);
RWFreq = random(RWFreqDist,Nsim,1);

% P-wave Single Blast Damping
PWXimin = 100;
PWXimax = 200;
PWXiDist = makedist('uniform',PWXimin,PWXimax);
PWXi = random(PWXiDist,Nsim,1);

% S-wave Single Blast Damping
SWXimin = 100;
SWXimax = 200;
SWXiDist = makedist('uniform',SWXimin,SWXimax);
SWXi = random(SWXiDist,Nsim,1);

% R-wave Single Blast Damping
RWXimin = 100;
RWXimax = 200;
RWXiDist = makedist('uniform',RWXimin,RWXimax);
RWXi = random(RWXiDist,Nsim,1);


for j = 1:size(SiteTable,1)
    [IM] = get_blast_sequence(BlastTable,PPV(j,:),chaos(j,:),PArrivalTime(j,:),SArrivalTime(j,:),RArrivalTime(j,:),PWFreq(k),SWFreq(j),RWFreq(j),PWXi(j),SWXi(j),RWXi(j),TSPar);
end


end

function [IM,AT,AF,VT,VF,UT,UF] = get_blast_sequence(BlastTable,PPV,chaos,PArrivalTim,SArrivalTime,RArrivalTime,PWFreq,SWFreq,RWFreq,PWXi,SWXi,RWXi,CompWaveEnergy,ShearWaveEnergy,SurfaceWaveEnergy,TSPar)


f = TSPar.f;
NUP = TSPar.NUP;
w = 2*pi*f;

% Build seed signal for P, S y R
PVFo = get_damped_armonic(f,0,PWFreq,PWXi);
SVFo = get_damped_armonic(f,0,SWFreq,SWXi);
RVFo = get_damped_armonic(f,0,RWFreq,RWXi);

% Energy of each signal
PRMS = sqrt(2*(PVFo(1:NUP).*w)*(PVFo(1:NUP).*w).');
SRMS = sqrt(2*(SVFo(1:NUP).*w)*(SVFo(1:NUP).*w).');
RRMS = sqrt(2*(RVFo(1:NUP).*w)*(RVFo(1:NUP).*w).');

M = repmat([PRMS SRMS RRMS],3,1).*[(1-CompWaveEnergy) CompWaveEnergy CompWaveEnergy; ShearWaveEnergy (1-ShearWaveEnergy) ShearWaveEnergy ;SurfaceWaveEnergy SurfaceWaveEnergy (1-SurfaceWaveEnergy)];
FS = M\[0;0;0];

PVFo = PVFo*FS(1);
SVFo = SVFo*FS(2);
RVFo = RVFo*FS(3);


end

function [F] = get_damped_armonic(f,to,fo,xi)
w = 2*pi*f;
wo = 2*pi*fo;
F = 1./(-w.^2+2*1i*xi*w*wo+wo^2).*exp(-1i*2*pi*f*to);
end

function [] = get_gaussian_noise(f,to,fo,cov)
% Normal noise with complex mean and deviation
mu = fo;
sigma = mu*cov;
AF = exp(-(f-mu).^2/(2*sigma^2))/sigma.*exp(-1i*2*pi*f*to);
end