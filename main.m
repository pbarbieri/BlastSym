function [] = main(XlsFile,fo,xi,sigmaDelay,Vw,PPV,SiteX,SiteY,R)
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
BlastSeqTable = cell2struct(BlastSeqTable(2:end,:),BlastSeqTable(1,:),2);
NBlast = numel(BlastSeqTable);
RandomDelay = random(makedist('normal',0,sigmaDelay),NBlast,1);
for k = 1:NBlast
    BlastSeqTable(k).fo = fo;
    BlastSeqTable(k).xi = xi;
    if k>1
        BlastSeqTable(k).T = (BlastSeqTable(k).T + RandomDelay(k));
    end
end
BlastSeqTable = struct2table(BlastSeqTable);
BlastSeqTable.T = BlastSeqTable.T/1000;

% build blast sequence
[UT,VT,AT,t] = get_regular_blast_sequence(BlastSeqTable,SiteX,SiteY,Vw,PPV,R);

close all
figure(1)
hold on
scatter(BlastSeqTable.X,BlastSeqTable.Y,'ok','filled');
scatter([SiteX,SiteX],[SiteY,SiteY],'^r','filled');
hold off
grid on


figure(2)
plot(t,UT)
grid on
xlabel('t [s]');
ylabel('U [m]');

figure(3)
plot(t,VT)
grid on
xlabel('t [s]');
ylabel('V [m/s]');


figure(4)
plot(t,AT)
grid on
xlabel('t [s]');
ylabel('A [m/s/s]');

end


