clear all
%% Single blast parameters
fo = 800;
xi = 0.25;
Delay = 8/1000;
sigmaDelay = 0.05*Delay;

%% Blast sequence table
% Number of blast in x dir
Nx = 10;
% Number of blast in y dir
Ny = 10;
% Distance between blasts [m]
S = 10;

BlastSeqTable = build_test_blast_sequence(Nx,Ny,S,Delay,[1],50);

BlastSeqTable = table2struct(BlastSeqTable);
NBlast = numel(BlastSeqTable);
RandomDelay = random(makedist('normal',0,sigmaDelay),NBlast,1);
for k = 1:NBlast
    BlastSeqTable(k).fo = fo;
    BlastSeqTable(k).xi = xi;
    if k>1
        BlastSeqTable(k).T = BlastSeqTable(k).T + RandomDelay(k);
    end
end
BlastSeqTable = struct2table(BlastSeqTable);

%% Wave propagation velocity
Vw = 2500;

%% PPV [m/s]
PPV = 25/1000;

%% Site coordiantes
SiteX = 257;
SiteY = 136;
close all
figure(1)
hold on
scatter(BlastSeqTable.X,BlastSeqTable.Y,'ok','filled');
scatter([SiteX,SiteX],[SiteY,SiteY],'^r','filled');
hold off
grid on

%% Blast sequence
[VT,AT,VF,AF,t,f] = get_regular_blast_sequence(BlastSeqTable,SiteX,SiteY,Vw,PPV);

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


% AT = AT(t<0.14);
% t = t(t<0.14);
% AT = AT(t>0.1);
% t = t(t>0.1);
% t = t-t(1);


