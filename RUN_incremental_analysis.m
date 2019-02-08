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
SeedGen.mode = 'p'; % 'psr' %p' %s' %r'
% Generator
SeedGen.ND = '1D'; % '1D' '2D'

% P-wave main frequency
SeedGen.fp.mu =  800; % [Hz]
SeedGen.fp.sigma = 0; % [Hz] 0 for constant value = mean
% P-wave damping
SeedGen.xip.mu = .6; % [ ]
SeedGen.xip.sigma = 0; % [ ] 0 for constant value = mean
% P-wave energy participation
SeedGen.Ep = 0.05; % [% over the total energy transmited by the ground in a single explotion]

% S-wave main frequency
SeedGen.fs.mu = 600; % [Hz]
SeedGen.fs.sigma = 20; % [Hz] 0 for constant value = mean
% S-wave damping
SeedGen.xis.mu = .4; % [ ]
SeedGen.xis.sigma = 0.05; % [ ] 0 for constant value = mean
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
BlastModel.R.sigma = 0.05; % [ ] 0 for constant value = mean
% Screening flag
BlastModel.Sflag = 1;
% Screening function
BlastModel.Sfun = @(Ns,D) 1/(1+23337*sqrt(Ns)/D^2);
% Global Scale Factor
BlastModel.GSF = 0; % [mm/s] 0 or ~=0. if GSF=0 no global PPV scaling is applied




%% RUN
PPV = [1 2 3 4 5 10 15 20 25 30 40 50 60 70 80 90 100 125 150 175 200 225 250 275 300 325 350];
NPPV = numel(PPV);
for j = 1:NPPV
    Run.ID = ['Solomon TSF1 PPV ',num2str(PPV(j)),' sin'];
    BlastModel.GSF = PPV(j); % [mm/s] 0 or ~=0. if GSF=0 no global PPV scaling is applied
    % Seed function
    SeedGen.fun = 'sin';
    % Output folder
    Run.OutPutFolder = fullfile(pwd,Run.ID);
    % Number of simulations
    Run.Nsim = 10;
    % Export to slide flag
    Run.ETS = 1;
    % plotting flag
    Run.Plot = 1;
    % Run
    [~,~] = main(Run,AttModel,Site,SeedGen,BlastModel);
end
for j = 1:NPPV
    Run.ID = ['Solomon TSF1 PPV ',num2str(PPV(j)),' blair'];
    BlastModel.GSF = PPV(j); % [mm/s] 0 or ~=0. if GSF=0 no global PPV scaling is applied
    % Seed function
    SeedGen.fun = 'blair';
    % Output folder
    Run.OutPutFolder = fullfile(pwd,Run.ID);
    % Number of simulations
    Run.Nsim = 10;
    % Export to slide flag
    Run.ETS = 1;
    % plotting flag
    Run.Plot = 1;
    % Run
    [~,~] = main(Run,AttModel,Site,SeedGen,BlastModel);
end



