function [AF,f] = Get_FS(AT,t)
% FOURIER TRANSFORM OF RECORD(S)
%   AT: Timeseries of record(s) by column.
%   t:  Time vector by column.
%   AF: Frequency content of record(s) by column.
%   f:  Frequency vector by column.

NP = size(AT,1);
dt = t(2)-t(1);
NFFT = pow2(nextpow2(NP)+1);
NUP = NFFT/2+1;
f = 1/(2*dt)*linspace(0,1,NUP).';
AF = fft(AT,NFFT,1);
AF = AF(1:NUP,:);
end
