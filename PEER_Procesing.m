function [AT,VT,UT] = PEER_Procesing(AT,t,hpf)
% PEER PROCESSING OF ACELERATION TIMESERIE(s) TO OBTAIN VELOCITY AND
% DISPLACEMENTS
%   AT: Timeseries of record(s) by column.
%   t:  Time vector by column.
%   VT: Velocity timeseries of record(s) by column.
%   UT: Displacement timeseries of record(s) by column.

dt = t(2)-t(1);
[HPFilterObj,LPFilterObj] = Build_Filters(dt,hpf);
% Filtered acceleration
[AT] = Filter_AT(AT,t,HPFilterObj,LPFilterObj);
% Baseline correction
[AT,VT,UT] = BaselineCorrection(AT,t);
end

function [HPFilterObj,LPFilterObj] = Build_Filters(dt,hpf)
% BUILDS FILTERS FOR RECORD PROCESSING
%   dt: Time step of the record
%   HPFilterObj: High-pass filter object
%   LPFilter: : Low-pass filter object

HPFilterObj = designfilt('highpassiir', 'FilterOrder', 6, ...
                 'HalfPowerFrequency',hpf, 'SampleRate',1/dt, ...
                 'DesignMethod', 'butter');
LPFilterObj = 1;
% if 1/(2*dt)<1000
%     lpf = 1/(2*dt);
% else
%     lpf = 1000;
% end
% LPFilterObj =  designfilt('lowpassiir', 'PassbandFrequency', lpf-5, ...
%             'StopbandFrequency', lpf, 'PassbandRipple', 1, ...
%             'StopbandAttenuation', lpf-10, 'SampleRate', 1/dt, 'MatchExactly', 'passband');
end

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
%     AT(:,c) = filtfilt(LPFilterObj,AT(:,c));
end
% Unpaddig
AT = AT(NPad+1:NPad+NP);
end

function [AT,VT,UT] = BaselineCorrection(AT,t)
% Newmark beta integration
[~,UT] = Get_VUT(AT,t);
% Base-line correction
for r = 1:size(AT,2)
    % 6th order polynom fit to fit displacement with ao=a1=0
    p = PolyFit(t,UT(:,r),[1 1 1 1 0 0]);
    % Acceleration correction
    AT(:,r) = AT(:,r) - polyval(polyder(polyder(p)),t);
end
% Newmark beta integration
[VT,UT] = Get_VUT(AT,t);
end


function [p] = PolyFit(x,y,coefs)
% returns the coefficients for a polynomial p(x) that is a ...
% best fit (in a least-squares sense) for the data in y. coefs is a vector
% of 1 or 0 such that the polynom is defined by:
%
% n = numel(coef)
% p(x) = coef(1)*p(1)*x^(n-1)+coef(2)*p(2)*x^(n-2)+....

NP = numel(x);
N = numel(coefs); % order of the polynom
n = sum(coefs); % number of nonzero coef
M = zeros(n,n);
b = zeros(n,1);
p = zeros(N,1); % coeficients

fil = 0;
for r = 1:N
if coefs(r)==1
    fil = fil+1;
    for j = 1:NP
        b(fil) = b(fil) + y(j) * x(j)^(N-fil);
    end
    col = 0;
    for k = 1:N
    if coefs(k)==1
       col = col+1;
       pwr_fil = x.^(N-k);
       pwr_col = x.^(N-r);
       M(fil,col) = sum(pwr_fil.*pwr_col);
    end
    end
end
end

a = M\b;
k = 0;
for j = 1:N
    if coefs(j)==1
        k = k+1;
        p(j) = a(k);
    end
end
end