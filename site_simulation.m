function [SiteTable] = site_simulation(Nsim)
%% ========================================================================
% Copyright SRK/FIUBA (C) 2018
% Coded By: P. Barbieri (pbarbieri@fi.uba.ar)
% Version: V103
% Date: 19-01-2019
% -------------------------------------------------------------------------
% USAGE: 
%   site_simulation: Simulates site properties
%   
% -------------------------------------------------------------------------
% INPUT:
%       Nsim    Number of site to simulate
% -------------------------------------------------------------------------
% OUTPUT:
%       
% -------------------------------------------------------------------------
% EXAMPLE:
%   
% -------------------------------------------------------------------------
% BIBLIO:
%   [1] Miller -  On the partition of energy between elastic waves in a 
%       semi-infinite solid
%   [2] Kausel - Fundamental soutions of elasto-dynamics
% -------------------------------------------------------------------------
% VALIDATE:
% Version:
% Date:
% Validated by:
% -------------------------------------------------------------------------
% LOG
%   V101    24/01/2019      First version
% =========================================================================

% Pressure wave velocity  [m/s]
Vpmin = 2800;
Vppeak = 3000;
Vpmax = 3200;
VpDist = makedist('triangular','a',Vpmin,'b',Vppeak,'c',Vpmax);
SiteTable.Vp = random(VpDist,Nsim,1);

% Shear wave velocity [m/s]
Vsmin = 1800;
Vspeak = 2000;
Vsmax = 2200;
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

% Raleigh wave velocity [m/s] see ref. [2]
SiteTable.Vr = (0.054.*nu.^4-0.08*nu.^3-0.038*nu.^2+0.195*nu+0.874).*Vs;
SiteTable = struct2table(SiteTable);

% Energy divition between waves - see ref. [1]
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