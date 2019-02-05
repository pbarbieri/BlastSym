%% ATTENUATION MODEL
a  = -1.6863;
b = 6.2693;
muLnPPV = @(r)a*log(r)+b;
sigmaLnPPV = @(r) 0.4880;
AttModel.muLnPPV = muLnPPV;
AttModel.sigmaLnPPV = sigmaLnPPV;

%% SITE CHARACTERISTICS
% Pressure wave veloity
Site.Cp = 4000; %[m/s]
% Shear wave velocity
Site.Cs = 2500; %[m/s]
% Poisson coeficient
Site.nu = 0.2;
% Rayleigh wave veloicty

%% SEED GENERATOR CHARACTERISTICS

% Seed function
SeedGen.fun = 'sine'; % 'gauss' %'file'
% Generator mode
SeedGen.mode = 'single wave'; % 'psr'
% Generator
SeedGen.ND = '2D'; % '1D'

% P-wave main frequency
SeedGen.fp = 800; % [Hz]
% P-wave damping
SeedGen.xip = 0.25; % []
% P-wave energy participation
SeedGen.Ep = 0.5; % [% over the total energy transmited by the ground in a single explotion]

% S-wave main frequency
SeedGen.fs = 600; % [Hz]
% S-wave damping
SeedGen.xis = 0.15; % []
% S-wave energy participation
SeedGen.Es = 28; % [% over the total energy transmited by the ground in a single explotion]

% R-wave main frequency
SeedGen.fr = 400; % [Hz]
% R-wave damping
SeedGen.xir = 0.1; % []
% R-wave energy participation
SeedGen.Er = 67; % [% over the total energy transmited by the ground in a single explotion]

%% BLASTING MODEL
% Blast sequence
BlastModel.SequenceFile = 'SEQ02.xlsx';
% Site Coordiantes
BlastModel.SiteX = 48.2;
BlastModel.SiteY = -145;
% Delay standard deviation 
BlastModel.SigmaDelay = 2; %[ms]
% Ratio of chaos in the signal
BlastModel.R = 0.8;
% Screening flag
BlastModel.Sflag = 0;
% Screening function
BlastModel.Sfun = @(Ns,D) 1/(1+23337*sqrt(Ns)/D^2);

%% RUN
Run.ID = '';
% Output folder
Run.OutPutFolder = '';
% plotting flag
Run.Plot = 1;
% Number of simulations
Run.Nsim = 100;
% Run
main_r2(Run,AttModel,Site,SeedGen,BlastModel);
