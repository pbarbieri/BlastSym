function [] = main(XlsFile,fo,xi,sigmaDelay,Vw,PPV,SiteX,SiteY,R,OutputFileName)
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
%       XlsFile         Excel file with detonation sequence
%       fo              Predominant frequency of indvidual blast [Hz]
%       xi              Damping of indvidual blast [Hz]   
%       sigmaDelay      Standar deviation of the delay [ms]
%       Vw              Velocity of wave propagation
%       PPV             PPV [m/s]
%       SiteX           Site X coordinate
%       SiteY           Site Y coordiante
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
% Read data
[~,~,BlastSeqTable] = xlsread(XlsFile);
% Build Blast sequence table
BlastSeqTable = cell2struct(BlastSeqTable(2:end,:),BlastSeqTable(1,:),2);
BlastSeqTable = struct2table(BlastSeqTable);
NBlast = size(BlastSeqTable,1);
Yblast = BlastSeqTable.Y;
Xblast = BlastSeqTable.X;
T = BlastSeqTable.T;

D = sqrt((SiteX-Xblast).^2+(SiteY-Yblast).^2);
fo = ones(NBlast,1)*fo;
xi = ones(NBlast,1)*xi;
RandomDelay = random(makedist('normal',0,sigmaDelay),NBlast,1);
RandomDelay(1) = 0; % Avoid T<0 in first blast
T = (T + RandomDelay)/1000;
to = T+D/Vw;

BlastSeqTable.T = BlastSeqTable.T/1000;
BlastSeqTable.to = to;
BlastSeqTable.D = D;
BlastSeqTable.fo = fo;
BlastSeqTable.xi = xi;

% Build blast sequence
[UT,VT,AT,t] = get_regular_blast_sequence(BlastSeqTable,PPV,R);

save([OutputFileName,'.mat'],'BlastSeqTable','UT','VT','AT','t');

% Plots
close all
hfig = figure(1);
set(hfig,'Color',[1 1 1],'Position',[50,50,470*1.5,400*1.5]);
hold on
scatter(BlastSeqTable.X,BlastSeqTable.Y,'ok','filled');
scatter([SiteX,SiteX],[SiteY,SiteY],'^r','filled');
hold off
grid on
legend({'Blast Holes','Site'},'location','southwest','FontSize',12);
xlabel('X [m]','FontSize',12);
ylabel('Y [m]','FontSize',12);

hfig = figure(2);
set(hfig,'Color',[1 1 1],'Position',[100,100,1000,300]);
plot(t,UT)
grid on
xlabel('t [s]');
ylabel('U [m]');
set(gca,'Position',[0.07,0.14,0.85,0.8]);

hfig = figure(3);
set(hfig,'Color',[1 1 1],'Position',[150,150,1000,300]);
plot(t,VT)
grid on
xlabel('t [s]');
ylabel('V [m/s]');
set(gca,'Position',[0.07,0.14,0.85,0.8]);

hfig = figure(4);
set(hfig,'Color',[1 1 1],'Position',[200,200,1000,300]);
plot(t,AT)
grid on
xlabel('t [s]');
ylabel('A [m/s/s]');
set(gca,'Position',[0.07,0.14,0.85,0.8]);

end


