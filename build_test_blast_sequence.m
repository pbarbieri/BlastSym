function [BlastSeq] = build_test_blast_sequence(Nx,Ny,S,Delay,Px,Charge)
%% ========================================================================
% Copyright SRK/FIUBA (C) 2019
% Coded By: P. Barbieri (pbarbieri@fi.uba.ar)
% Version: V103
% Date: 19-01-2019
% -------------------------------------------------------------------------
% USAGE: 
%   build_test_blast_sequence: builds an speceific blasts sequence
%   
% -------------------------------------------------------------------------
% INPUT:
%       Nx          Number of blast-holes in the x-x axis
%       Ny          Number of blast-holes in the y-y axis
%       S           Distance between blast holes
%       D           Depth of charges
%       Delaty      Delay between blasts
%       Charge      Base charge
%       Px          Probability of charge multipliyer
% -------------------------------------------------------------------------
% OUTPUT:
%       BlastSeq    Blast sequence table
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


% total number of blasts
NBlast = Nx*Ny;

% Charge per blast
Px(1) = 1-sum(Px(2:end));
Px = cumsum(Px);
ChargeMultiSeed = rand(NBlast,1);
ChargeMulti = NaN(NBlast,1);
ChargeMulti(ChargeMultiSeed<=Px(1)) = 1;
for k = 2:numel(Px)
    ChargeMulti(and(Px(k-1)<ChargeMultiSeed,ChargeMultiSeed<=Px(k))) = k;
end

% Blast Sequence table
BlastSeq = struct;
X = linspace(0,Nx-1,Nx).'*S; X= X-(Nx-1)*S/2; X = repmat(X,1,Ny);
Y = linspace(0,Ny-1,Ny)*S; Y = Y-(Ny-1)*S/2; Y = repmat(Y,Nx,1);
BlastSeq.X = X(:);
BlastSeq.Y = Y(:);
BlastSeq.D = zeros(NBlast,1);
BlastSeq.W = Charge*ChargeMulti;
BlastSeq.T = linspace(0,NBlast-1,NBlast).'*Delay;
BlastSeq.S = false(NBlast,1);
BlastSeq = struct2table(BlastSeq);


end