function [BlastTable] = blast_vibration_simulation(BlastTable,BlastAttModel,SiteProp,SiteX,SiteY,Nsim,OutputFolder)
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
%                       - W [kg] explosive mass
%                       - T [ms] time of explotion
%                       - S [ ]  secrrening flag (1 to consider screening)
%                       - D [m] depth of blast
%       BlastAttModel   Blast attenuation model defined by muLnPPV(R) and
%                       sigmalLnPPV(R) in [mm/s2] where R [m/kg^-05] is the
%                       adimetional distance.
%       SiteProp        Sturct with site properties
%       SiteX           Site X coordinate
%       SiteY           Site Y coordiante
%       Nsim            Number of scenario simulations
%       Seed            
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
BlastTable = screening(BlastTable,SiteX,SiteY,BlastAttModel,SiteProp);

% Simulation of wave existation for each site
[BlastTable] = wave_simulation(BlastTable,SiteProp,Nsim,OutputFolder);



end

function [BlastTable,NBlast] = screening(BlastTable,SiteX,SiteY,AttenuationModel,SiteProp)

% Blast detonation data
NBlast = size(BlastTable,1);
X = BlastTable.X;
Y = BlastTable.Y;
T = BlastTable.T;
W = BlastTable.W;
BlastTable = table2struct(BlastTable);

% Attenuation model
muLnPPV = AttenuationModel.muLnPPV;
sigmaLnPPV = AttenuationModel.sigmaLnPPV;

% Site properties
Vp = SiteProp.Vp;
Vs = SiteProp.Vs;
Vr = SiteProp.Vr;


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
    R = Dist./sqrt(W(k));
    BlastTable(k).R = R;
    BlastTable(k).muLnPPV = muLnPPV(R);
    BlastTable(k).sigmaLnPPV = sigmaLnPPV(R);
    
    % Wave travel times
    BlastTable(k).tp = Dist/Vp;
    BlastTable(k).ts = Dist/Vs;
    BlastTable(k).tr = Dist/Vr;
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
TSPar.dt = dt;
end

function [SiteTable] = wave_simulation(BlastTable,SiteTable,Nsim,OutputFolder)
% Set timeseries parameters
[TSPar] = set_timeseries_parameters(BlastTable);

NBlast = size(BlastTable,1);

% Attenuaiton
muLnPPV = BlastTable.muLnPPV;
sigmaLnPPV = [BlastTable.sigmaLnPPV];
NsigmaPPV = random(makedist('normal',0,1),Nsim,1);
PPV = zeros(Nsim,NBlast);
for k = 1:Nsim
    for j = 1:NBlast
        PPV(k,j) = exp(muLnPPV(j)+NsigmaPPV(k)*sigmaLnPPV(j));
    end 
end

% Delay simulation
muDelay = 0; %[s]
sigmaDelay = 0.1/1000; %[s]
Delay = random(makedist('normal',muDelay,sigmaDelay),Nsim,NBlast) + repmat(BlastTable.T.',Nsim,1);

% Arrival P wave time
PArrivalTime = repmat(BlastTable.Dist.',Nsim,1)./repmat(SiteTable.Vp,1,NBlast) + Delay;
% Arrival S wave time
SArrivalTime = repmat(BlastTable.Dist.',Nsim,1)./repmat(SiteTable.Vs,1,NBlast) + Delay;
% Arrival R wave time
RArrivalTime = repmat(BlastTable.Dist.',Nsim,1)./repmat(SiteTable.Vr,1,NBlast) + Delay;

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
PWXimin = 0.1;
PWXimax = 0.2;
PWXiDist = makedist('uniform',PWXimin,PWXimax);
PWXi = random(PWXiDist,Nsim,1);

% S-wave Single Blast Damping
SWXimin = 0.05;
SWXimax = 0.1;
SWXiDist = makedist('uniform',SWXimin,SWXimax);
SWXi = random(SWXiDist,Nsim,1);

% R-wave Single Blast Damping
RWXimin = 0.01;
RWXimax = 0.05;
RWXiDist = makedist('uniform',RWXimin,RWXimax);
RWXi = random(RWXiDist,Nsim,1);




for j = 1:size(SiteTable,1)
    get_blast_sequence(BlastTable,PPV(j,:),PArrivalTime(j,:),SArrivalTime(j,:),RArrivalTime(j,:),PWFreq(k),SWFreq(j),RWFreq(j),PWXi(j),SWXi(j),RWXi(j),SiteTable.CompWaveEnergy,SiteTable.ShearWaveEnergy,SiteTable.SurfaceWaveEnergy,TSPar);
end


end

function [ ] = get_blast_sequence(BlastTable,PPV,PArrivalTim,SArrivalTime,RArrivalTime,PWFreq,SWFreq,RWFreq,PWXi,SWXi,RWXi,CompWaveEnergy,ShearWaveEnergy,SurfaceWaveEnergy,TSPar)

f = TSPar.f;
NUP = TSPar.NUP;
w = 2*pi*f;

% Build seed signal for P, S y R
PVFo = get_damped_armonic(f,0,PWFreq,PWXi);
SVFo = get_damped_armonic(f,0,SWFreq,SWXi);
RVFo = get_damped_armonic(f,0,RWFreq,RWXi);

[PATo,~] = Get_TS(PVFo.*w,f);
[SATo,~] = Get_TS(SVFo.*w,f);
[RATo,~] = Get_TS(RVFo.*w,f);

[FSp,FSs,FSr] = get_scale_parameters(PATo,SATo,RATo,CompWaveEnergy,ShearWaveEnergy,SurfaceWaveEnergy);


% Energy of each signal
PRMS = sqrt(2*dot((PVFo(1:NUP).*w),(PVFo(1:NUP).*w).'));
SRMS = sqrt(2*dot((SVFo(1:NUP).*w),(SVFo(1:NUP).*w).'));
RRMS = sqrt(2*dot((RVFo(1:NUP).*w),(RVFo(1:NUP).*w).'));

M = repmat([PRMS SRMS RRMS],3,1).*[(1-CompWaveEnergy) CompWaveEnergy CompWaveEnergy; ShearWaveEnergy (1-ShearWaveEnergy) ShearWaveEnergy ;SurfaceWaveEnergy SurfaceWaveEnergy (1-SurfaceWaveEnergy)];
FS = M\[0;0;0];

PVFo = PVFo*FS(1); [PVTo,~] = Get_TS(PVFo,f);
SVFo = SVFo*FS(2); [SVTo,~] = Get_TS(SVFo,f);
RVFo = RVFo*FS(3); [RVTo,~] = Get_TS(RVFo,f);

% Scale to unit amplitude
Vmax = max([max(abs(PVTo)) max(abs(SVTo)) max(abs(RVTo))]);
PVFo = PVFo/Vmax;
SVFo = SVFo/Vmax;
RVFo = RVFo/Vmax;


NBlast = size(BlastTable,1);
% Tangential wave amplitude (% of longitudinal wave)
TangAmp = random(makedist('uniform',-0.5,0.5),NBlast,1);


% Chaotic number generator
r = zeros(NUP,NBlast);
chaos = 0.54;
for j = 1:NBlast
    for k = 1:NUP
        r(k,j) = chaos;
        chaos = 2*chaos^2-1;
    end
end
Rj = sum(r,1);

R = 0.8; % Relative ammount of noise energy




end

function [FSp,FSs,FSr] = get_scale_parameters(PA,SA,RA,CompWaveEnergy,ShearWaveEnergy,SurfaceWaveEnergy)

PA = PA/sqrt(dot(PA,PA));
SA = SA/sqrt(dot(SA,SA));
RA = RA/sqrt(dot(RA,RA));

PP = dot(PA,PA);
PS = dot(PA,SA);
PR = dot(PA,RA);
SS = dot(SA,SA);
SR = dot(SA,RA);
RR = dot(RA,RA);

RSMo = dot(PA+SA+RA,PA+SA+RA);

G = @(x) PP*x(1)^2 + SS*x(2)^2 + RR*x(3)^2 + 2*PS*x(1)*x(2) + 2*PR*x(1)*x(3) + 2*SR*x(2)*x(3);
F = @(x) [PP*x(1)^2/CompWaveEnergy^2 ; SS*x(2)^2/ShearWaveEnergy^2 ; RR*x(3)^2/SurfaceWaveEnergy^2] - ones(3,1)*G(x);
J = @(x) 2*[PP*x(1)/CompWaveEnergy^2 ; SS*x(2)/ShearWaveEnergy^2 ; RR*x(3)/SurfaceWaveEnergy] - 2*[PP,PS,PR;PS,SS,SR;PR,SR,RR]*[x(1);x(2);x(3)];

x = NaN(3,100);
a = [PP*(CompWaveEnergy^2-1) ; SS*(ShearWaveEnergy^2-1); RR*(SurfaceWaveEnergy^2-1)];
b = [2*CompWaveEnergy^2*(PS+PR); 2*ShearWaveEnergy^2*(PS+SR); 2*SurfaceWaveEnergy^2*(PR+SR) ] ;
c = [CompWaveEnergy^2*(SS+RR+2*SR);ShearWaveEnergy^2*(PP+RR+2*PR); SurfaceWaveEnergy^2*(PP+SS+2*PS) ];
x(:,1) = max([(-b+sqrt(b.^2-4.*a.*c))./(2.*a),(-b-sqrt(b.^2-4.*a.*c))./(2.*a)],[],2);
for k = 2:100
    x(:,k) = x(:,k-1)-J(x(:,k-1))\F(x(:,k-1));
end

FSp = x(1);
FSs = x(2);
FSr = x(3);

close all
hold on
plot(x(1,:),'-r');
plot(x(2,:),'-b');
plot(x(3,:),'-k');
hold off
grid on


(CompWaveEnergy - sqrt(x(1,end)^2*PP/G(x(:,end))))/CompWaveEnergy
(ShearWaveEnergy - sqrt(x(2,end)^2*SS/G(x(:,end))))/ShearWaveEnergy
(SurfaceWaveEnergy - sqrt(x(3,end)^2*RR/G(x(:,end))))/SurfaceWaveEnergy




end

function [F] = get_damped_armonic(f,to,fo,xi)
w = 2*pi*f;
wo = 2*pi*fo;
F = 1./(w.^2-2*1i*xi*w*wo-wo^2).*exp(-1i*w*to);
end

function [] = get_gaussian_noise(f,to,fo,cov)
% Normal noise with complex mean and deviation
mu = fo;
sigma = mu*cov;
AF = exp(-(f-mu).^2/(2*sigma^2))/sigma.*exp(-1i*2*pi*f*to);
end