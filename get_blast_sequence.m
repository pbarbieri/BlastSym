function [ ] = get_blast_sequence(BlastTable,PPV,PArrivalTim,SArrivalTime,RArrivalTime,PWFreq,SWFreq,RWFreq,PWXi,SWXi,RWXi,CompWaveEnergy,ShearWaveEnergy,SurfaceWaveEnergy,TSPar)

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
TangAmp = random(makedist('uniform',-0.5,0.5).NBlast,1);


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

function [F] = get_damped_armonic(f,to,fo,xi)
w = 2*pi*f;
wo = 2*pi*fo;
F = 1./(-w.^2+2*1i*xi*w*wo+wo^2).*exp(-1i*2*pi*f*to);
end