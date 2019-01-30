function [AT] = Filter_AT(AT,t,HPFilterObj,LPFilterObj)
% FILTER OF ACELERATION TIMESERIE(s)
%   AT: Timeseries of record(s) by column.
%   t:  Time vector by column.
%   HPFilterObj: High-pass filter object
%   LPFilter: : Low-pass filter object

NP =  size(AT,1);
dt = t(2) - t(1);
fhp = HPFilterObj.HalfPowerFrequency;
nfilt = HPFilterObj.FilterOrder;
% Padding 
NPad = ceil(1.5/2*nfilt/fhp*1/dt);
AT(NPad+1:NPad+NP,:) = AT(:,:);
AT(1:NPad,:) = 0;
AT(NPad+NP+1:2*NPad+NP,:) = 0;
% Filtering
for c = 1:size(AT,2)
    AT(:,c) = filtfilt(HPFilterObj,AT(:,c));
    AT(:,c) = filtfilt(LPFilterObj,AT(:,c));
end
% Unpaddig
AT = AT(NPad+1:NPad+NP);
end