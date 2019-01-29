clear all
%% Single blast parameters
fo = 800;
xi = 0.25;
Delay = 8/1000;
sigmaDelay = 0.05*Delay;

%% Blast sequence table
% Number of blast in x dir
Nx = 1;
% Number of blast in y dir
Ny = 2;
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
SiteX = 300;
SiteY = 0;
close all
figure(1)
hold on
scatter(BlastSeqTable.X,BlastSeqTable.Y,'ok','filled');
scatter([SiteX,SiteX],[SiteY,SiteY],'^r','filled');
hold off
grid on

%% Blast sequence
[VT,VF,t,f] = get_regular_blast_sequence(BlastSeqTable,SiteX,SiteY,Vw,PPV);
AF = VF.*(2*pi*f);
[AT,~] = Get_TS(AF,f);

figure(2)
plot(f,abs(VF))
xlim([0 1000])
grid on
xlabel('f [Hz]');
ylabel('V [m/s]');


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



