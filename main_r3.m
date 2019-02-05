function [] = main_r2(Run,AttModel,Site,SeedGen,BlastModel)
%% ========================================================================
% Copyright SRK/FIUBA (C) 2019
% Coded By: P. Barbieri (pbarbieri@fi.uba.ar)
% Version: V103
% Date: 30-01-2019
% -------------------------------------------------------------------------
% USAGE: 
%   main: Simulates blast vibrations
%   
% -------------------------------------------------------------------------
% INPUT:
%       Run         
%       AttModel    
%       Site         
%       SeedGen     
%       BlastModel   
% -------------------------------------------------------------------------
% OUTPUT:
%       
% -------------------------------------------------------------------------
% EXAMPLE:
%   
% -------------------------------------------------------------------------
% BIBLIO:
%   [1] D. Blair, "Blast vibration control in the presence of delay scatter
%       and random fluctuations between blastholes," International journal 
%       for numerical and analytical methods in geomechanics
%   [2] D. Blair, "Statistical models for ground vibration and airblast,"
%       Fragblast, vol. 3, no. 4, pp. 335--364, 1999
% -------------------------------------------------------------------------
% VALIDATE:
% Version:
% Date:
% Validated by:
% -------------------------------------------------------------------------
% LOG
%   V101    30/01/2019      First version
% =========================================================================

clc
% Build blast table
[BlastSeqTable] = build_blast_sequence_table(BlastModel,AttModel,Site,Run.Nsim);

% Build blast sequence
[UT,VT,AT,t] = bluid_blast_sequence(BlastSeqTable,PPV,R);

% save([OutputFileName,'.mat'],'BlastSeqTable','UT','VT','AT','t');
% write_csv_output(AT,t,OutputFileName);



end

function [BlastSeqTable] = build_blast_sequence_table(BlastModel,AttModel,Site,Nsim)

SiteX = BlastModel.SiteX;
SiteY = BlastModel.SiteY;
Sfun = BlastModel.Sfun;
Cp = Site.Cp;
Cs = Site.Cs;
Cr = Site.Cr;
muLnPPv = AttModel.muLnPPV;
sigmaLnPPv = AttModel.sigmaLnPPv;

% Read data
[~,~,BlastSeqTable] = xlsread(BlastModel.SequenceFile);
% Build Blast sequence table
BlastSeqTable = cell2struct(BlastSeqTable(2:end,:),BlastSeqTable(1,:),2);
BlastSeqTable = struct2table(BlastSeqTable);
NBlast = size(BlastSeqTable,1);
Yblast = BlastSeqTable.Y;
Xblast = BlastSeqTable.X;
T = BlastSeqTable.T;
W = BlastSeqTable.W;

% Distances and geometry
D = sqrt((SiteX-Xblast).^2+(SiteY-Yblast).^2);
R = D./sqrt(W);
% unitary direction vector from site to blast
l = [ones(Nblast,1)*SiteX, ones(Nblast,1)*SiteX]-[Xblast,Yblast]; % Versor blast 2 site
l = l./sqrt(l(:,1).^2+l(:,2).^2);
n = [-l(2,:),l(1,:)];    
% directing cosines for radial vibration
alpha = atan(l(:,2)/l(:,1));
Rx = cos(alpha);
Ry = sin(alpha);
% directing cosines for tangential displacement
Tx = cos(alpha+pi/2);
Ty = sin(alpha+pi/2);

% Arrival times
RandomDelay = random(makedist('normal',0,BlastModel.sigmaDelay),NBlast,Nsim);
RandomDelay(1,:) = 0; % Avoid T<0 in first blast
T = (repmat(T,1,Nsim) + RandomDelay)/1000;
tp = T+repmat(D,1,Nsim)/Cp;
ts = T+repmat(D,1,Nsim)/Cs;
tr = T+repmat(D,1,Nsim)/Cr;

% Screening widht
Dmat = zeros(NBlast,NBlast);
for k = 1:NBlast
    for j = 1:NBlast
        Dmat(k,j) = sqrt((X(k)-X(j))^2+(Y(k)-Y(j))^2);
    end
end
Dmat = triu(Dmat);
Dmat = Dmat(:);
Dmat = Dmat(Dmat>0);
ScreeningWidth = min(Dmat); % minumun distnace between contiguos blasholes

Ns = NaN(NBlast,1);
S = ones(NBlast,1);
if BlastModel.Sflag
    for k = 1:NBlast
        % Remove future blasts
        X_filtered = Xblast; X_filtered = X_filtered(T<T(k));
        Y_filtered = Yblast; Y_filtered = Y_filtered(T<T(k));
        % Remove blasts from same borehole
        idx = and(X_filtered~=X(k),Y_filtered~=Y(k));
        X_filtered = X_filtered(idx);
        Y_filtered = Y_filtered(idx);
        % Screening polygon
        Coord = zeros(4,2);
        Coord(1,:) = [SiteX,SiteY]-n*ScreeningWidth;
        Coord(2,:) = [X(k),Y(k)]-n*ScreeningWidth;
        Coord(3,:) = [X(k),Y(k)]+n*ScreeningWidth;
        Coord(4,:) = [SiteX,SiteY]+n*ScreeningWidth;
        % Count blasts
        [in,on] = inpolygon(X_filtered,Y_filtered,Coord(:,1),Coord(:,2));
        idx = or(in,on);
        Ns(k) = sum(idx);
        % Screening
        S(k) = Sfun(Ns(k),D(k));
    end
end

% Ratio between the amplitude of the transversal and longitudinal vibration
Tranv2LongRatio = random(makedist('uniform',-0.5,0.5),NBlast,1);


epsilon = random(makedist('normal',0,1),1,Nsim);
epsilon = repmat(epsilon,NBlast,1);
PPV = exp(muLnPPV(repmat(R,1,Nsim))+epsilon*sigmaLnPPV(repmat(R,1,Nsim)));


BlastSeqTable.T = BlastSeqTable.T/1000;
BlastSeqTable.tp = tp;
BlastSeqTable.ts = ts;
BlastSeqTable.tr = tr;
BlastSeqTable.D = D;
BlastSeqTable.R = R;
BlastSeqTable.l = l;
BlastSeqTable.n = n;
BlastSeqTable.alpha = alpha;
BlastSeqTable.Rx = Rx;
BlastSeqTable.Ry = Ry;
BlastSeqTable.Tx = Tx;
BlastSeqTable.Ty = Ty;
BlastSeqTable.S = S;
BlastSeqTable.Tranv2LongRatio = Tranv2LongRatio;
BlastSeqTable.PPV = PPV;

% fo = ones(NBlast,1)*fo;
% xi = ones(NBlast,1)*xi;
% a  = -1.6863;
% b = 6.2693;
% muLnPPV = @(r)a*log(r)+b;
% sigmaLnPPV = 0.4880;
% epsilon = random(makedist('normal',0,1),1,1);
% PPV = exp(muLnPPV(R)+epsilon*sigmaLnPPV);
% BlastSeqTable.fo = fo;
% BlastSeqTable.xi = xi;
% BlastSeqTable.PPV = PPV;

end

function [UT,VT,AT,t] = bluid_blast_sequence(BlastSeqTable,PPV,R)
%% ========================================================================
% Copyright SRK/FIUBA (C) 2018
% Coded By: P. Barbieri (pbarbieri@fi.uba.ar)
% Version: V103
% Date: 19-01-2019
% -------------------------------------------------------------------------
% USAGE: 
%   blast_vibration_simulation: Simulates blast vibrations
%   
% -------------------------------------------------------------------------
% INPUT:
%       BlastTable      Table with the blast sequence with fields:
%                       - X & Y [m] coordinates of the blast
%                       - T [s] time of explotion
%                       - fo [Hz] natural frec of each blast
%                       - xi []   dampling coefficient of each blast
%       PPV             PPV [m/s]
%       R               % of randomization
% -------------------------------------------------------------------------
% OUTPUT:
%       
% -------------------------------------------------------------------------
% EXAMPLE:
%   
% -------------------------------------------------------------------------
% BIBLIO:
%   [1] D. Blair, "Blast vibration control in the presence of delay scatter
%       and random fluctuations between blastholes," International journal 
%       for numerical and analytical methods in geomechanics
%   [2] D. Blair, "Statistical models for ground vibration and airblast,"
%       Fragblast, vol. 3, no. 4, pp. 335--364, 1999
% -------------------------------------------------------------------------
% VALIDATE:
% Version:
% Date:
% Validated by:
% -------------------------------------------------------------------------
% LOG
%   V101    19/01/2019      First version
% =========================================================================

NBlast = size(BlastSeqTable,1);
to = BlastSeqTable.to;
fo = BlastSeqTable.fo;
xi = BlastSeqTable.xi;
PVV = BlastSeqTable.PPV;


% Set timeseries parameters
Tmax = max(BlastSeqTable.to)*1.2;
dt = 10^(floor(log10(1/(8*min(fo(fo>0))))));
NPo = ceil(Tmax/dt)+1;
t = linspace(0,(NPo-1)*dt,NPo).';
Tpulse = log(0.001)/(-2*pi*min(fo)*min(xi));
NPpulse = ceil(Tpulse/dt);
tpulse = linspace(0,(NPpulse-1)*dt,NPpulse).';

% Single blast generatorn
Delta = 0.01; % remaining amplitud at -t1
t1 = 10*dt;
k = log(1/Delta-1)/t1;
L = @(t) 1./(1+exp(-k*t));
Upulse = @(t,to,fo,xi) sin(2*pi*fo*(t-to-t1)).*exp(-xi*2*pi*fo*(t-to-t1)).*L(t-to-t1);

Vpulse = @(t,to,fo,xi) (2*pi*fo*cos(2*pi*fo*(t-to-t1)).*min(exp(-xi*2*pi*fo*(t-to-t1)),1000).*L(t-to-t1)...
                       - xi*2*pi*fo*sin(2*pi*fo*(t-to-t1)).*min(exp(-xi*2*pi*fo*(t-to-t1)),1000).*L(t-to-t1)...
                       + sin(2*pi*fo*(t-to-t1)).*min(exp(-xi*2*pi*fo*(t-to-t1)),1000).*k.*L(t-to-t1).*(1-L(t-to-t1)));

Apulse = @(t,to,fo,xi) -(2*pi*fo).^2*sin(2*pi*fo*(t-to-t1)).*exp(-xi*2*pi*fo*(t-to-t1)).*L(t-to-t1)...
                       -(2*pi*fo).^2*xi*cos(2*pi*fo*(t-to-t1)).*exp(-xi*2*pi*fo*(t-to-t1)).*L(t-to-t1)...
                       + 2*pi*fo*cos(2*pi*fo*(t-to-t1)).*exp(-xi*2*pi*fo*(t-to-t1)).*k.*L(t-to-t1).*(1-L(t-to-t1))...
                       - xi*(2*pi*fo)^2*cos(2*pi*fo*(t-to-t1)).*exp(-xi*2*pi*fo*(t-to-t1)).*L(t-to-t1)...
                       + (xi*2*pi*fo)^2*sin(2*pi*fo*(t-to-t1)).*exp(-xi*2*pi*fo*(t-to-t1)).*L(t-to-t1)...
                       - xi*2*pi*fo* sin(2*pi*fo*(t-to-t1)).*exp(-xi*2*pi*fo*(t-to-t1)).*k.*L(t-to-t1).*(1-L(t-to-t1))... 
                       + 2*pi*fo*cos(2*pi*fo*(t-to-t1)).*exp(-xi*2*pi*fo*(t-to-t1)).*k.*L(t-to-t1).*(1-L(t-to-t1))...
                       - 2*pi*fo*xi*sin(2*pi*fo*(t-to-t1)).*exp(-xi*2*pi*fo*(t-to-t1)).*k.*L(t-to-t1).*(1-L(t-to-t1))...
                       + sin(2*pi*fo*(t-to-t1)).*exp(-xi*2*pi*fo*(t-to-t1)).*k.^2.*(L(t-to-t1).*(1-L(t-to-t1)).^2-L(t-to-t1).^2.*(1-L(t-to-t1)));        

% Build blast sequecnce
VT = zeros(NPo,NBlast);
if R>0
    % Randomizer
    r = zeros(NBlast*NPo,1);
    r(1) = 0.52; for k = 2:NBlast*NPo, r(k) = 2*r(k-1)^2-1; end
    r = reshape(r,NPo,NBlast);
    for k = 1:NBlast
        vo = Vpulse(t,to(k),fo(k),xi(k));
        v1 = conv(r(:,k),vo);
        VT(:,k) = (1-R)*vo+ R/sum(r(:,k).^2)*v1(1:NPo);
        VT(:,k) = PVV(k)*VT(:,k)/max(abs(VT(:,k)));
    end
    VT = sum(VT,2);
    UT = zeros(NPo,1);
    AT = zeros(NPo,1);
    for k = 2:NPo-1
        AT(k) = (VT(k+1)-VT(k-1))/(2*dt);
        UT(k) = UT(k-1)+(1-1/2)*dt*VT(k-1)+1/2*dt*VT(k);
    end
    AT(NPo) = (VT(end)-VT(end-1))/dt;
    UT(NPo) = UT(NPo-1)+(1-1/2)*dt*VT(NPo-1)+1/2*dt*VT(NPo);
else
    UT = zeros(NPo,NBlast);
    AT = zeros(NPo,NBlast);
    for k = 1:NBlast
        VT(:,k) = Vpulse(t,to(k),fo(k),xi(k));
        FS = PPV(k)./max(abs(VT(:,k)));
        VT(:,k) = FS*VT(:,k);
        UT(:,k) = FS*Upulse(t,to(k),fo(k),xi(k));
        AT(:,k) = FS*Apulse(t,to(k),fo(k),xi(k));
    end 
    UT = sum(UT,2);
    VT = sum(VT,2);
    AT = sum(AT,2);
end
FS = PPV/max(abs(VT));
VT = VT*FS;
AT = AT*FS;
UT = UT*FS;

% Remove unncesary pre-pading
tini = min(to)*0.8;
UT = UT(t>tini);
VT = VT(t>tini);
AT = AT(t>tini);
t = t(t>tini); t = t-t(1);
end

function [] = write_csv_output(AT,t,OutputFileName)
idx = find(round(AT,5)~=0,1);
if idx>1
    AT = AT(idx-1:end);
    t = t(1:numel(AT));
end
NP = size(AT,1);
AT(1) = 0;
AT = AT/9.81;
AT(round(AT,6)==0) = 0;

FileName = [OutputFileName,'.txt'];
fiad = fopen(FileName,'w');
for k = 1:NP
    fprintf(fiad,'%f\t%f\n',t(k),AT(k));
end
fclose(fiad);

end



