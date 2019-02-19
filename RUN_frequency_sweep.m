clc
clear all

%% ATTENUATION MODEL
a  = -1.6863;
b = 6.2693;
muLnPPV = @(r)a*log(r)+b;
sigmaLnPPV = @(r) 0.4880;
AttModel.muLnPPV = muLnPPV;
AttModel.sigmaLnPPV = sigmaLnPPV;

%% SITE CHARACTERISTICS
% Pressure wave veloity
Site.Cp.mu = 4000; %[m/s]
Site.Cp.sigma = 100; %[m/s] 0 for constant value = mean
% Shear wave velocity
Site.Cs.mu = 2500; %[m/s]
Site.Cs.sigma = 100; %[m/s] 0 for constant value = mean
% Poisson coeficient
Site.nu.mu = 0.2;
Site.nu.sigma = 0.01; %[ ] 0 for constant value = mean
% Rayleigh wave veloicty

%% SEED GENERATOR CHARACTERISTICS

% Seed function
SeedGen.fun = 'sin'; %'sin' 'blair' %'gauss' %'file'
% Generator mode
SeedGen.mode = 's'; % 'psr' %p' %s' %r'
% Generator
SeedGen.ND = '1D'; % '1D' '2D'

% P-wave main frequency
SeedGen.fp.mu =  20; % [Hz]
SeedGen.fp.sigma = 0; % [Hz] 0 for constant value = mean
% P-wave damping
SeedGen.xip.mu = .01; % [ ]
SeedGen.xip.sigma = 0; % [ ] 0 for constant value = mean
% P-wave energy participation
SeedGen.Ep = 0.05; % [% over the total energy transmited by the ground in a single explotion]

% S-wave main frequency
SeedGen.fs.mu = 20; % [Hz]
SeedGen.fs.sigma = 0; % [Hz] 0 for constant value = mean
% S-wave damping
SeedGen.xis.mu = .1; % [ ]
SeedGen.xis.sigma = 0; % [ ] 0 for constant value = mean
% S-wave energy participation
SeedGen.Es = 0.28; % [% over the total energy transmited by the ground in a single explotion]

% R-wave main frequency
SeedGen.fr.mu = 400; % [Hz]
SeedGen.fr.sigma = 20; % [Hz]
% R-wave damping
SeedGen.xir.mu = 0.1; % [ ]
SeedGen.xir.sigma = 0.03; % [ ] 0 for constant value = mean
% R-wave energy participation
SeedGen.Er = 0.67; % [% over the total energy transmited by the ground in a single explotion]

%% BLASTING MODEL
% Blast sequence
BlastModel.SequenceFile = 'SEQ02.xlsx';
% Site Coordiantes
BlastModel.SiteX = 48.2;
BlastModel.SiteY = -145;
% Delay standard deviation 
BlastModel.SigmaDelay = 2; %[ms]
% Ratio of chaos in the signal
BlastModel.R.mu = 0.8; % [ ]
BlastModel.R.sigma = 0; % [ ] 0 for constant value = mean
% Screening flag
BlastModel.Sflag = 0;
% Screening function
BlastModel.Sfun = @(Ns,D) 1/(1+23337*sqrt(Ns)/D^2);
% Global Scale Factor
BlastModel.GSF = 30; % [mm/s] 0 or ~=0. if GSF=0 no global PPV scaling is applied


%% RUN
fs = [5 7 9 13 15 17 19 21 23 25];
Nfs = numel(fs);
for j = 1:Nfs
    Run.ID = ['Solomon TSF1 fo ',num2str(fs(j)),' sin'];
    SeedGen.fs.mu = fs(j);
    % Seed function
    SeedGen.fun = 'sin';
    % Output folder
    Run.OutPutFolder = fullfile(pwd,Run.ID);
    % Number of simulations
    Run.Nsim = 1;
    % Export to slide flag
    Run.ETS = 1;
    % plotting flag
    Run.Plot = 1;
    % Run
    [~,~] = main(Run,AttModel,Site,SeedGen,BlastModel);
end
for j = 1:Nfs
    Run.ID = ['Solomon TSF1 fo ',num2str(fs(j)),' blair'];
    SeedGen.fs.mu = fs(j);
    % Seed function
    SeedGen.fun = 'blair';
    % Output folder
    Run.OutPutFolder = fullfile(pwd,Run.ID);
    % Number of simulations
    Run.Nsim = 1;
    % Export to slide flag
    Run.ETS = 1;
    % plotting flag
    Run.Plot = 1;
    % Run
    [~,~] = main(Run,AttModel,Site,SeedGen,BlastModel);
end

get_incremental_plots;




