function [AT,t] = Get_TS(AF,f)
% INVERSE FOURIER TRANSFORM OF RECORD(S)
%   AF: Frequency content of record(s) by column.
%   f:  Frequency vector by column.
%   AT: Timeseries of record(s) by column.
%   t:  Time vector by column.

df = f(2)-f(1);
NUP = size(AF,1);
dt = 1/(2*df*(NUP-1));
t = linspace(0,(NUP-1)*dt,NUP).';
AF = cat(1,AF,fliplr(AF));
AT = ifft(AF,'symmetric');
AT = AT(1:NUP,:);
end