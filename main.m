function [RCTable,BlastSeqTable] = main(Run,AttModel,Site,SeedGen,BlastModel)
%% ========================================================================
% Copyright SRK/FIUBA (C) 1984
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
OutPutFolder = Run.OutPutFolder;
fprintf(1,'Blast vibration simulation ID =  %s \n',Run.ID);
if isdir(OutPutFolder), rmdir(OutPutFolder,'s'); end
mkdir(OutPutFolder);
% Scenario simulation
fprintf(1,'\t> Simulation of properties for %i scenarios. \n',Run.Nsim);
[RCTable] = build_record_table(Site,SeedGen,BlastModel,Run);

% Build blast table
[BlastSeqTable] = build_blast_sequence_table(BlastModel,AttModel,RCTable,Run.Nsim);

% Build blast sequence
if strcmpi(SeedGen.ND,'1D')
    [RCTable] = build_1D_sequence(BlastSeqTable,RCTable,OutPutFolder,SeedGen.mode,SeedGen.fun,BlastModel.GSF,Run.ETS);
elseif strcmpi(SeedGen.ND,'2D')
    [RCTable] = build_2D_sequence(BlastSeqTable,RCTable,OutPutFolder,SeedGen.mode,SeedGen.fun,BlastModel.GSF,Run.ETS);
else
    keyboard;
end

% Save output
fprintf(1,'\t> Saving data.... \n');
save(fullfile(OutPutFolder,'Run.mat'),'Run');
save(fullfile(OutPutFolder,'RCTable.mat'),'RCTable');
save(fullfile(OutPutFolder,'BlastSeqTable.mat'),'BlastSeqTable');
save(fullfile(OutPutFolder,'Site.mat'),'Site');
save(fullfile(OutPutFolder,'SeedGen.mat'),'SeedGen');
save(fullfile(OutPutFolder,'BlastModel.mat'),'BlastModel');

fprintf(1,'Done! \n');


end

function [RCTable] = build_record_table(Site,SeedGen,BlastModel,Run)
Nsim = Run.Nsim;

RCID = [repmat([Run.ID,' - RC '],Nsim,1), num2str(linspace(1,Nsim,Nsim).','%.5i')];

% Site geotechnical properties
if Site.Cp.sigma>0
    Cp = random(makedist('normal',Site.Cp.mu,Site.Cp.sigma),Nsim,1);
else
    Cp = Site.Cp.mu*ones(Nsim,1);
end

if Site.Cs.sigma>0
    Cs = random(makedist('normal',Site.Cs.mu,Site.Cs.sigma),Nsim,1);
else
    Cs = Site.Cs.mu*ones(Nsim,1);
end

if Site.nu.sigma>0
    nu = random(makedist('normal',Site.nu.mu,Site.nu.sigma),Nsim,1);
else
    nu = Site.nu.mu*ones(Nsim,1);
end
Cr = (0.054.*nu.^4-0.08*nu.^3-0.038*nu.^2+0.195*nu+0.874).*Cs;

% Seed generator porperties
if SeedGen.fp.sigma>0
    fpo = random(makedist('normal',SeedGen.fp.mu,SeedGen.fp.sigma),Nsim,1);
else
    fpo = SeedGen.fp.mu*ones(Nsim,1);
end


if SeedGen.xip.sigma>0
    xip = random(makedist('normal',SeedGen.xip.mu,SeedGen.xip.sigma),Nsim,1);
else
    xip = SeedGen.xip.mu*ones(Nsim,1);
end


if SeedGen.fs.sigma>0
    fso = random(makedist('normal',SeedGen.fs.mu,SeedGen.fs.sigma),Nsim,1);
else
    fso = SeedGen.fs.mu*ones(Nsim,1);
end

if SeedGen.xis.sigma>0
    xis = random(makedist('normal',SeedGen.xis.mu,SeedGen.xis.sigma),Nsim,1);
else
    xis = SeedGen.xis.mu*ones(Nsim,1);
end


if SeedGen.fr.sigma>0
    fro = random(makedist('normal',SeedGen.fr.mu,SeedGen.fr.sigma),Nsim,1);
else
    fro = SeedGen.fr.mu*ones(Nsim,1);
end

if SeedGen.xir.sigma>0
    xir = random(makedist('normal',SeedGen.xir.mu,SeedGen.xir.sigma),Nsim,1);
else
    xir = SeedGen.xir.mu*ones(Nsim,1);
end

% Ratio of chaos in the signal
if BlastModel.R.sigma>0
    R = random(makedist('normal',BlastModel.R.mu,BlastModel.R.sigma),Nsim,1);
else
    R = BlastModel.R.mu*ones(Nsim,1);
end
R = max(R,0);

% Attenuation model
epsilonPPV = random(makedist('normal',0,1),Nsim,1);


RCTable = struct;
RCTable.RCID = RCID;
RCTable.Cp = Cp;
RCTable.Cs = Cs;
RCTable.Cr = Cr;
RCTable.nu = nu;

RCTable.fpo = fpo;
RCTable.xip = xip;
RCTable.fso = fso;
RCTable.xis = xis;
RCTable.fro = fro;
RCTable.xir = xir;

RCTable.R = R;
RCTable.epsilonPPV = epsilonPPV;

RCTable = struct2table(RCTable);
end

function [BlastSeqTable] = build_blast_sequence_table(BlastModel,AttModel,RCTable,Nsim)

SiteX = BlastModel.SiteX;
SiteY = BlastModel.SiteY;
Sfun = BlastModel.Sfun;

muLnPPV = AttModel.muLnPPV;
sigmaLnPPV = AttModel.sigmaLnPPV;

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


Cp = repmat(RCTable.Cp.',NBlast,1);
Cs = repmat(RCTable.Cs.',NBlast,1);
Cr = repmat(RCTable.Cr.',NBlast,1);

% Distances and geometry
D = sqrt((SiteX-Xblast).^2+(SiteY-Yblast).^2);
SD = D./sqrt(W);
% unitary direction vector from site to blast
l = [ones(NBlast,1)*SiteX, ones(NBlast,1)*SiteX]-[Xblast,Yblast]; % Versor blast 2 site
l = l./sqrt(l(:,1).^2+l(:,2).^2);
n = [-l(:,2),l(:,1)];    
% directing cosines for radial vibration
alpha = atan(l(:,2)./l(:,1));
Rx = cos(alpha);
Ry = sin(alpha);
% directing cosines for tangential displacement
Tx = cos(alpha+pi/2);
Ty = sin(alpha+pi/2);

% Arrival times
RandomDelay = random(makedist('normal',0,BlastModel.SigmaDelay),NBlast,Nsim);
RandomDelay(1,:) = 0; % Avoid T<0 in first blast
T = (repmat(T,1,Nsim) + RandomDelay)/1000;
tp = T+repmat(D,1,Nsim)./Cp;
ts = T+repmat(D,1,Nsim)./Cs;
tr = T+repmat(D,1,Nsim)./Cr;

% Screening widht
Dmat = zeros(NBlast,NBlast);
for k = 1:NBlast
    for j = 1:NBlast
        Dmat(k,j) = sqrt((Xblast(k)-Xblast(j))^2+(Yblast(k)-Yblast(j))^2);
    end
end
Dmat = triu(Dmat);
Dmat = Dmat(:);
ScreeningWidth = min(Dmat(Dmat>0)); % minumun distnace between contiguos blasholes

Ns = NaN(NBlast,Nsim);
S = ones(NBlast,Nsim);
if BlastModel.Sflag
    for s = 1:Nsim
        for k = 1:NBlast
                % Remove future blasts
                X_filtered = Xblast; X_filtered = X_filtered(T(:,s)<T(k,s));
                Y_filtered = Yblast; Y_filtered = Y_filtered(T(:,s)<T(k,s));
                % Remove blasts from same borehole
                idx = and(X_filtered~=Xblast(k),Y_filtered~=Yblast(k));
                if ~isempty(idx)
                X_filtered = X_filtered(idx);
                Y_filtered = Y_filtered(idx);
                % Screening polygon
                Coord = zeros(4,2);
                Coord(1,:) = [SiteX,SiteY]-n(k,:)*ScreeningWidth;
                Coord(2,:) = [Xblast(k),Yblast(k)]-n(k,:)*ScreeningWidth;
                Coord(3,:) = [Xblast(k),Yblast(k)]+n(k,:)*ScreeningWidth;
                Coord(4,:) = [SiteX,SiteY]+n(k,:)*ScreeningWidth;
                % Count blasts
                [in,on] = inpolygon(X_filtered,Y_filtered,Coord(:,1),Coord(:,2));
                idx = or(in,on);
                Ns(k,s) = sum(idx);
                % Screening
                S(k,s) = Sfun(Ns(k,s),D(k));
            else
                Ns(k,s) = 0;
                S(k,s) = 1;
            end
        end
    end
end

% Ratio between the amplitude of the transversal and longitudinal vibration
Tranv2LongRatio = random(makedist('uniform',-0.5,0.5),NBlast,Nsim);

% Attenuation model
epsilonPPV = repmat(RCTable.epsilonPPV.',NBlast,1);
PPV = exp(muLnPPV(repmat(SD,1,Nsim))+epsilonPPV.*sigmaLnPPV(repmat(SD,1,Nsim)));

BlastSeqTable.T = BlastSeqTable.T;
BlastSeqTable.tp = tp;
BlastSeqTable.ts = ts;
BlastSeqTable.tr = tr;
BlastSeqTable.D = D;
BlastSeqTable.SD = SD;
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
end

function [Vpulse] = get_pulse_generator(GenFun,dt)

switch lower(GenFun)
    case {'sine','sin'}
        Delta = 0.01; % remaining amplitud at -t1
        t1 = 10*dt;
        k = log(1/Delta-1)/t1;
        L = @(t) 1./(1+exp(-k*t));
        Vpulse = @(t,to,fo,xi) (2*pi*fo*cos(2*pi*fo*(t-to-t1)).*min(exp(-xi*2*pi*fo*(t-to-t1)),1000).*L(t-to-t1)...
                               - xi*2*pi*fo*sin(2*pi*fo*(t-to-t1)).*min(exp(-xi*2*pi*fo*(t-to-t1)),1000).*L(t-to-t1)...
                               + sin(2*pi*fo*(t-to-t1)).*min(exp(-xi*2*pi*fo*(t-to-t1)),1000).*k.*L(t-to-t1).*(1-L(t-to-t1)));
    case {'gauss'}

    case {'file'}

    case {'blair'}
        L = @(t) (sign(t)+1)/2;
        Vpulse = @(t,to,fo,xi) exp(-max((t-to),0)*fo*10).*(((t-to)*fo*10).^6-12*((t-to)*fo*10).^5+30*((t-to)*fo*10).^4).*L((t-to));
    otherwise
        keyboard;
end 


end

function [RCTable] = build_1D_sequence(BlastSeqTable,RCTable,OutputFoler,Mode,GenFun,GSF,ETS)

% Flag for wave types
PFlag = 1;
SFlag = 1;
RFlag = 1;
if isempty(find(lower(Mode)=='p',1)), PFlag = 0; end
if isempty(find(lower(Mode)=='s',1)), SFlag = 0; end
if isempty(find(lower(Mode)=='r',1)), RFlag = 0; end

NBlast = size(BlastSeqTable,1);
tp = BlastSeqTable.tp;
ts = BlastSeqTable.ts;
tr = BlastSeqTable.tr;
PVV = BlastSeqTable.PPV/1000; % mm2m


Nsim = size(RCTable,1);
RCID = RCTable.RCID;
fpo = RCTable.fpo;
xip = RCTable.xip;
% Ep = RCTable.Ep;
fso = RCTable.fso;
xis = RCTable.xis;
% Es = RCTable.Ep;
fro = RCTable.fro;
xir = RCTable.xir;
% Er = RCTable.Er;
R = RCTable.R;

% Set timeseries parameters
fo = min([min(fpo)/PFlag,min(fso)/SFlag,min(fro)/RFlag]);
Tmax = max(max(BlastSeqTable.tr))*1.2;
dt = 10^(floor(log10(1/(8*min(fo(fo>0))))));
NPo = ceil(Tmax/dt)+1;
t = linspace(0,(NPo-1)*dt,NPo).';


% Select pulse generator function
[Vpulse] = get_pulse_generator(GenFun,dt);

% Blast generator
PGA = zeros(Nsim,1);
PGV = zeros(Nsim,1);
PGD = zeros(Nsim,1);
IA = zeros(Nsim,1);
RMSA = zeros(Nsim,1);
RMSV = zeros(Nsim,1);
RMSD = zeros(Nsim,1);
ky = logspace(-3,1,20);
Newmark1 = zeros(Nsim,20);
Newmark2 = zeros(Nsim,20);
fprintf(1,'\t> Simulation of %i 1D %s records with %s seed.\n',Nsim,Mode,GenFun);
for s = 1:Nsim
    fprintf(1,'\t\t- Building record %i of %i.\n',s,Nsim);
    % Build blast sequecnce
    VT = zeros(NPo,NBlast);
    UT = zeros(NPo,1);
    AT = zeros(NPo,1);
    % Randomizer
    r = zeros(NBlast*NPo,1);
    r(1) = 0.52; for k = 2:NBlast*NPo, r(k) = 2*r(k-1)^2-1; end
    r = reshape(r,NPo,NBlast);
    for k = 1:NBlast
        Vp = Vpulse(t,tp(k,s),fpo(s),xip(s)); Vp = PFlag*Vp/max(abs(Vp));
        Vs = Vpulse(t,ts(k,s),fso(s),xis(s)); Vs = SFlag*Vs/max(abs(Vs));
        Vr = Vpulse(t,tr(k,s),fro(s),xir(s)); Vr = RFlag*Vr/max(abs(Vr));
        Vo = Vp+Vs+Vr;
        V1 = conv(r(:,k),Vo);
        VT(:,k) = (1-R(s))*Vo+ R(s)/sum(r(:,k).^2)*V1(1:NPo);
        VT(:,k) = PVV(k,s)*VT(:,k)/max(abs(VT(:,k)));
    end
    VT = sum(VT,2);
    if GSF~=0, VT = GSF/(1000*max(abs(VT)))*VT;  end
    for k = 2:NPo-1
        AT(k) = (VT(k+1)-VT(k-1))/(2*dt);
        UT(k) = UT(k-1)+(1-1/2)*dt*VT(k-1)+1/2*dt*VT(k);
    end
    AT(NPo) = (VT(end)-VT(end-1))/dt;
    UT(NPo) = UT(NPo-1)+(1-1/2)*dt*VT(NPo-1)+1/2*dt*VT(NPo);
    % Intensity measures
    [IM] = Get_IM(t,AT,VT,UT,ky);
    PGA(s) = IM.PGA;
    PGV(s) = IM.PGV;
    PGD(s) = IM.PGD;
    IA(s) = IM.IA;
    RMSA(s) = IM.RMSA;
    RMSV(s) = IM.RMSV;
    RMSD(s) = IM.RMSD;
    Newmark1(s,:) = IM.Newmark1;
    Newmark2(s,:) = IM.Newmark2;
    FileName = fullfile(OutputFoler,[RCID(s,:),'.mat']);
    save(FileName,'AT','VT','UT','t');
    if ETS
        FileName = fullfile(OutputFoler,[RCID(s,:),'.txt']);
        export2slide(t,AT,FileName);
    end
end
% Save IM in RCTable
RCTable.NP =  ones(Nsim,1)*NPo;
RCTable.dt =  ones(Nsim,1)*dt;
RCTable.Tmax =  ones(Nsim,1)*Tmax;
RCTable.PGA = PGA;
RCTable.PGV = PGV;
RCTable.PGD = PGD;
RCTable.IA = IA;
RCTable.RMSA = RMSA;
RCTable.RMSV = RMSV;
RCTable.RMSD = RMSD;
RCTable.ky = repmat(ky,Nsim,1);
RCTable.Newmark1 = Newmark1;
RCTable.Newmark2 = Newmark2;

end

function [RCTable] = build_2D_sequence(BlastSeqTable,RCTable,OutputFoler,Mode,GenFun,GSF,ETS)
% Flag for wave types
PFlag = 1;
SFlag = 1;
RFlag = 1;
if isempty(find(lower(Mode)=='p',1)), PFlag = 0; end
if isempty(find(lower(Mode)=='s',1)), SFlag = 0; end
if isempty(find(lower(Mode)=='r',1)), RFlag = 0; end

NBlast = size(BlastSeqTable,1);
Rx = BlastSeqTable.Rx;
Ry = BlastSeqTable.Ry;
Tx = BlastSeqTable.Tx;
Ty = BlastSeqTable.Ty;
tp = BlastSeqTable.tp;
ts = BlastSeqTable.ts;
tr = BlastSeqTable.tr;
PVV = BlastSeqTable.PPV/1000; % mm2m
Tranv2LongRatio = BlastSeqTable.Tranv2LongRatio;

Nsim = size(RCTable,1);
RCID = RCTable.RCID;
fpo = RCTable.fpo;
xip = RCTable.xip;
% Ep = RCTable.Ep;
fso = RCTable.fso;
xis = RCTable.xis;
% Es = RCTable.Ep;
fro = RCTable.fro;
xir = RCTable.xir;
% Er = RCTable.Er;
R = RCTable.R;

% Set timeseries parameters
fo = min([min(fpo)/PFlag,min(fso)/SFlag,min(fro)/RFlag]);
Tmax = max(max(BlastSeqTable.tr))*1.2;
dt = 10^(floor(log10(1/(8*min(fo(fo>0))))));
NPo = ceil(Tmax/dt)+1;
t = linspace(0,(NPo-1)*dt,NPo).';

% Proyection matrices
Rx = repmat(Rx.',NPo,1);
Ry = repmat(Ry.',NPo,1);
Tx = repmat(Tx.',NPo,1);
Ty = repmat(Ty.',NPo,1);

% Select pulse generator function
[Vpulse] = get_pulse_generator(GenFun,dt);

% Blast generator
PGAx = zeros(Nsim,1);
PGVx = zeros(Nsim,1);
PGDx = zeros(Nsim,1);
IAx = zeros(Nsim,1);
RMSAx = zeros(Nsim,1);
RMSVx = zeros(Nsim,1);
RMSDx = zeros(Nsim,1);
ky = logspace(-3,1,20);
Newmark1x = zeros(Nsim,20);
Newmark2x = zeros(Nsim,20);
PGAy = zeros(Nsim,1);
PGVy = zeros(Nsim,1);
PGDy = zeros(Nsim,1);
IAy = zeros(Nsim,1);
RMSAy = zeros(Nsim,1);
RMSVy = zeros(Nsim,1);
RMSDy = zeros(Nsim,1);
Newmark1y = zeros(Nsim,20);
Newmark2y = zeros(Nsim,20);
fprintf(1,'\t> Simulation of %i 1D %s records with %s seed.\n',Nsim,Mode,GenFun);
for s = 1:Nsim
    fprintf(1,'\t\t- Building record %i of %i.\n',s,Nsim);
    % Build blast sequecnce
    VTlong = zeros(NPo,NBlast);
    VTtrans = zeros(NPo,NBlast);
    ATx = zeros(NPo,1);
    ATy = zeros(NPo,1);  
    UTx = zeros(NPo,1);
    UTy = zeros(NPo,1);
    % Randomizer
    r = zeros(NBlast*NPo,1);
    r(1) = 0.52; for k = 2:NBlast*NPo, r(k) = 2*r(k-1)^2-1; end
    r = reshape(r,NPo,NBlast);
    for k = 1:NBlast
        Vp = Vpulse(t,tp(k,s),fpo(s),xip(s)); Vp = PFlag*Vp/max(abs(Vp));
        Vs = Vpulse(t,ts(k,s),fso(s),xis(s)); Vs = SFlag*Vs/max(abs(Vs));
        Vr = Vpulse(t,tr(k,s),fro(s),xir(s)); Vr = RFlag*Vr/max(abs(Vr));
        Vo = Vp+Vs+Vr;
        V1 = conv(r(:,k),Vo);
        VTlong(:,k) = (1-R(s))*Vo+ R(s)/sum(r(:,k).^2)*V1(1:NPo);
        VTlong(:,k) = PVV(k,s)*VTlong(:,k)/max(abs(VTlong(:,k)));
        VTtrans(:,k) = Tranv2LongRatio(k,s)*VTlong(:,k);
    end
    VTx = sum(Rx.*VTlong+Tx.*VTtrans,2);
    VTy = sum(Ry.*VTlong+Ty.*VTtrans,2);
    if GSF~=0
        VTx = GSF/(1000*max(abs(VTx)))*VTx;
        VTy = GSF/(1000*max(abs(VTy)))*VTy;
    end
    for k = 2:NPo-1
        ATx(k) = (VTx(k+1)-VTx(k-1))/(2*dt);
        UTx(k) = UTx(k-1)+(1-1/2)*dt*VTx(k-1)+1/2*dt*VTx(k);
        ATy(k) = (VTy(k+1)-VTy(k-1))/(2*dt);
        UTy(k) = UTy(k-1)+(1-1/2)*dt*VTy(k-1)+1/2*dt*VTy(k);
    end
    ATx(NPo) = (VTx(end)-VTx(end-1))/dt;
    UTx(NPo) = UTx(NPo-1)+(1-1/2)*dt*VTx(NPo-1)+1/2*dt*VTx(NPo);
    ATy(NPo) = (VTy(end)-VTy(end-1))/dt;
    UTy(NPo) = UTy(NPo-1)+(1-1/2)*dt*VTy(NPo-1)+1/2*dt*VTy(NPo);
    
    % Intensity measures
    [IMx] = Get_IM(t,ATx,VTx,UTx,ky);
    [IMy] = Get_IM(t,ATy,VTy,UTy,ky);
    PGAx(s) = IMx.PGA;
    PGVx(s) = IMx.PGV;
    PGDx(s) = IMx.PGD;
    IAx(s) = IMx.IA;
    RMSAx(s) = IMx.RMSA;
    RMSVx(s) = IMx.RMSV;
    RMSDx(s) = IMx.RMSD;
    Newmark1x(s,:) = IMx.Newmark1;
    Newmark2x(s,:) = IMx.Newmark2;
    PGAy(s) = IMy.PGA;
    PGVy(s) = IMy.PGV;
    PGDy(s) = IMy.PGD;
    IAy(s) = IMy.IA;
    RMSAy(s) = IMy.RMSA;
    RMSVy(s) = IMy.RMSV;
    RMSDy(s) = IMy.RMSD;  
    Newmark1y(s,:) = IMy.Newmark1;
    Newmark2y(s,:) = IMy.Newmark2;
    FileName = fullfile(OutputFoler,[RCID(s,:),'.mat']);
    save(FileName,'ATx','VTx','UTx','ATy','VTy','UTy','t');
    if ETS
        FileName = fullfile(OutputFoler,[RCID(s,:),'_x.txt']);
        export2slide(t,ATx,FileName);
        FileName = fullfile(OutputFoler,[RCID(s,:),'_y.txt']);
        export2slide(t,ATx,FileName);
    end
end
% Save IM in RCTable
RCTable.NP =  ones(Nsim,1)*NPo;
RCTable.dt =  ones(Nsim,1)*dt;
RCTable.Tmax =  ones(Nsim,1)*Tmax;
RCTable.PGAx = PGAx;
RCTable.PGVx = PGVx;
RCTable.PGDx = PGDx;
RCTable.IAx = IAx;
RCTable.RMSAx = RMSAx;
RCTable.RMSVx = RMSVx;
RCTable.RMSDx = RMSDx;
RCTable.PGAy = PGAy;
RCTable.PGVy = PGVy;
RCTable.PGDy = PGDy;
RCTable.IAy = IAy;
RCTable.RMSAy = RMSAy;
RCTable.RMSVy = RMSVy;
RCTable.RMSDy = RMSDy;
RCTable.ky = repmat(ky,Nsim,1);
RCTable.Newmark1x = Newmark1x;
RCTable.Newmark2x = Newmark2x;
RCTable.Newmark1y = Newmark1y;
RCTable.Newmark2y = Newmark2y;
end

function [IM] = Get_IM(t,AT,VT,UT,ac)
g = 9.81;
NKy = numel(ac);
Newmark1 = zeros(NKy,1);
Newmark2 = zeros(NKy,1);
NP = numel(AT);
dt = t(2)-t(1);

IM.PGA = max(abs(AT));
IM.RMSA = sqrt(dot(AT,AT)/NP);
IM.IA = pi/2*g*dot(AT,AT)*dt;
IM.ky = ac;
for k = 1:NKy
    AT1 = AT-ac(k); AT1(AT<ac(k)) = 0;
    AT2 = AT+ac(k); AT2(AT>-ac(k)) = 0;
    [UT1,~] = Get_NVUT(AT1,t,ac(k));
    [UT2,~] = Get_NVUT(AT2,t,ac(k));
    Newmark1(k) = UT1(end);
    Newmark2(k) = UT2(end);
end
IM.Newmark1 = Newmark1;
IM.Newmark2 = Newmark2;

IM.PGV = max(abs(VT));
IM.RMSV = sqrt(dot(VT,VT)/NP);

IM.PGD = max(abs(UT));
IM.RMSD = sqrt(dot(UT,UT)/NP);

end

function [] = export2slide(t,AT,OutputFileName)
idx = find(round(AT,5)~=0,1);
if idx>1
    AT = AT(idx-1:end);
    t = t(1:numel(AT));
end
NP = size(AT,1);
AT(1) = 0;
AT = AT/9.81; %m/s2 2 g
% slide only reads 6 digits, negative accelerations <10^-6 results in a
% -0.00000 which crashes Slide.
AT(round(AT,6)==0) = 0; 

FileName = [OutputFileName,'.txt'];
fiad = fopen(FileName,'w');
for k = 1:NP
    fprintf(fiad,'%f\t%f\n',t(k),AT(k));
end
fclose(fiad);

end

function [VT,UT,t] = Get_NVUT(AT,t,ac)
% NEWMARK INTEGRATION OF ACELERATION TIMESERIE(s) FOR SLINDING BLOCK WITH
% CRICITAL ACCELERATION
%   AT: Acceleration timeseries of record(s) by column.
%   t:  Time vector by column.
%   VT: Velocity timeseries of record(s) by column.
%   UT: Displacement timeseries of record(s) by column.
%   ac: Critical acceleration

% Default Newmark integration Parameters
gamma = 1/2;
beta = 1/4;

AT = abs(AT);

NP = size(AT,1);
dt = t(2) - t(1);
VT = zeros(NP,size(AT,2));
UT = VT;
for j = 2:NP
    if AT(j)>0
        VT(j,:) = VT(j-1,:)+(1-gamma)*dt*AT(j-1,:)+gamma*dt*AT(j,:);
        UT(j,:) = UT(j-1,:)+dt*VT(j-1,:)+(1/2-beta)*dt^2*AT(j-1,:)+beta*dt^2*AT(j,:);
    else
        VT(j,:) = VT(j-1,:)+(1-gamma)*dt*(-ac)+gamma*dt*(-ac);
        if VT(j,:)<0, VT(j,:) = 0;  end
        UT(j,:) = UT(j-1,:)+dt*VT(j-1,:)+(1/2-beta)*dt^2*(-ac)+beta*dt^2*(-ac);
    end
    
end
j = NP;
while VT(j,:)>0
    j = j+1;
    VT(j,:) = VT(j-1,:)+(1-gamma)*dt*(-ac)+gamma*dt*(-ac);
    if VT(j,:)<0, VT(j,:) = 0;  end
    UT(j,:) = UT(j-1,:)+dt*VT(j-1,:)+(1/2-beta)*dt^2*(-ac)+beta*dt^2*(-ac);
    t(j) = t(j-1) + dt;
end

end
