function [] = blast_sequence_simulation(Nsim,OutputFolder)
%% ========================================================================
% Copyright SRK/FIUBA (C) 2018
% Coded By: P. Barbieri (pbarbieri@fi.uba.ar)
% Version: V103
% Date: 19-01-2019
% -------------------------------------------------------------------------
% USAGE: 
%   blast_sequence_simulation: Simulates single hole blast sequences
%   
% -------------------------------------------------------------------------
% INPUT:
%       Nsim            Number of blast sequences to simulate
%       OutputFolder    Folder where the blast sequences will be stored
% -------------------------------------------------------------------------
% OUTPUT:
%       
% -------------------------------------------------------------------------
% EXAMPLE:
%   
% -------------------------------------------------------------------------
% BIBLIO:
%
% -------------------------------------------------------------------------
% VALIDATE:
% Version:
% Date:
% Validated by:
% -------------------------------------------------------------------------
% LOG
%   V101    24/01/2019      First version
% =========================================================================

if isdir(OutputFolder), rmdir(OutputFolder,'s'); end
mkdir(OutputFolder);

Namin = 2;
Namax = 10;
Na = randi([Namin,Namax],Nsim,1);
Nbmin = 2;
Nbmax = 20;
Nb = randi([Nbmin,Nbmax],Nsim,1);
NBlast = Na.*Nb;
NBlastMax = max(NBlast);

% Distance between blasts
Rmin = 5;
Rpeak = 10;
Rmax = 15;
R = random(makedist('triangular','a',Rmin,'b',Rpeak,'c',Rmax),Nsim,1);


% Delay [ms]
DelayMin = 4;
DelayPeak = 8;
DelayMax = 12;
Delay = floor(random(makedist('triangular','a',DelayMin,'b',DelayPeak,'c',DelayMax),Nsim,1));

% Charge Multipliyer
Px(2) = 0.25; % probability of x2 multiplier
Px(3) = 0.1; % probability of x3 multiplier
Px(1) = 1-sum(Px);
Px = cumsum(Px);
ChargeMultiSeed = rand(Nsim,NBlastMax);
ChargeMulti = NaN(Nsim,NBlastMax);
ChargeMulti(ChargeMultiSeed<=Px(1)) = 1;
for k = 2:numel(Px)
    ChargeMulti(and(Px(k-1)<ChargeMultiSeed,ChargeMultiSeed<=Px(k))) = k;
end

% Charge 
ChargeMin = 10;
ChargePeak = 45;
ChargeMax = 120;
Charge = random(makedist('triangular','a',ChargeMin,'b',ChargePeak,'c',ChargeMax),Nsim,1);

for k = 1:Nsim
    BlastSeq = struct;
    X = linspace(0,Na(k)-1,Na(k)).'*R(k); X= X-(Na(k)-1)*R(k)/2; X = repmat(X,1,Nb(k));
    Y = linspace(0,Nb(k)-1,Nb(k))*R(k); Y = Y-(Nb(k)-1)*R(k)/2; Y = repmat(Y,Na(k),1);
    BlastSeq.X = X(:);
    BlastSeq.Y = Y(:);
    BlastSeq.D = zeros(NBlast(k),1);
    Multiplier = ChargeMulti(k,:).';
    Multiplier = Multiplier(1:NBlast(k));
    BlastSeq.W = Charge(k)*Multiplier;
    BlastSeq.T = linspace(0,NBlast(k)-1,NBlast(k)).'*Delay(k);
    BlastSeq.S = false(NBlast(k),1);
    BlastSeq = struct2table(BlastSeq);
    save(fullfile(OutputFolder,['BlastSeq_',num2str(k,'%05.f')]),'BlastSeq');
    clear BlastSeq
end


end