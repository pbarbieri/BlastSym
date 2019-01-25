function [VT,UT] = Get_VUT(AT,t)
% NEWMARK INTEGRATION OF ACELERATION TIMESERIE(s)
%   AT: Acceleration timeseries of record(s) by column.
%   t:  Time vector by column.
%   VT: Velocity timeseries of record(s) by column.
%   UT: Displacement timeseries of record(s) by column.

% Default Newmark integration Parameters
gamma = 1/2;
beta = 1/4;
    
NP = size(AT,1);
dt = t(2) - t(1);
VT = zeros(NP,size(AT,2));
UT = VT;
for j = 2:NP
    VT(j,:) = VT(j-1,:)+(1-gamma)*dt*AT(j-1,:)+gamma*dt*AT(j,:);
    UT(j,:) = UT(j-1,:)+dt*VT(j-1,:)+(1/2-beta)*dt^2*AT(j-1,:)+beta*dt^2*AT(j,:);
end

end