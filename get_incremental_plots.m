close all
fs = [5 7 9 13 15 17 19 21];
Nfs = numel(fs);
hfig = figure(1);
set(hfig,'color',[1 1 1]);
ha = tightPlots(Nfs, 2, 18, [4.5 1], [0.5 1], 1, 2, 'centimeters');
for j = 1:Nfs
    Run.ID = ['Solomon TSF1 fo ',num2str(fs(j)),' sin'];
    OutPutFolder = fullfile(pwd,Run.ID,[Run.ID,' - RC 00001.mat']);
    load(OutPutFolder);
    
    axes(ha(j+(j-1)));
    plot(t,VT,'-r');
    grid on
    xlim([0,2]);
    if j~=Nfs,xticklabels(''); else, xlabel('t [s]'); end
    
    Run.ID = ['Solomon TSF1 fo ',num2str(fs(j)),' blair'];
    OutPutFolder = fullfile(pwd,Run.ID,[Run.ID,' - RC 00001.mat']);
    load(OutPutFolder);
    
    axes(ha(2*j));
    plot(t,VT,'-r');
    grid on
    xlim([0,2]);
    if j~=Nfs,xticklabels(''); else, xlabel('t [s]'); end
    
end
Pos = hfig.Position;
Pos(2) = 20;
set(hfig,'Position',Pos);


hfig = figure(2);
set(hfig,'color',[1 1 1]);
ha = tightPlots(Nfs, 2, 18, [4.5 1], [0.5 1], 1, 2, 'centimeters');
for j = 1:Nfs
    Run.ID = ['Solomon TSF1 fo ',num2str(fs(j)),' sin'];
    OutPutFolder = fullfile(pwd,Run.ID,[Run.ID,' - RC 00001.mat']);
    load(OutPutFolder);
    
    axes(ha(j+(j-1)));
    [VF,f] = Get_FS(VT,t);
    plot(f,abs(VF),'-r');
    grid on
    xlim([0,200]);
    if j~=Nfs,xticklabels(''); else, xlabel('f [Hz]'); end
    
    Run.ID = ['Solomon TSF1 fo ',num2str(fs(j)),' blair'];
    OutPutFolder = fullfile(pwd,Run.ID,[Run.ID,' - RC 00001.mat']);
    load(OutPutFolder);
    
    axes(ha(2*j));
    [VF,f] = Get_FS(VT,t);
    plot(f,abs(VF),'-r');
    grid on
    xlim([0,200]);
    if j~=Nfs,xticklabels(''); else, xlabel('f [Hz]'); end
    
end
Pos = hfig.Position;
Pos(2) = 20;
set(hfig,'Position',Pos);


hfig = figure(3);
set(hfig,'color',[1 1 1]);
ha = tightPlots(Nfs, 2, 18, [4.5 1], [0.5 1], 1, 2, 'centimeters');
for j = 1:Nfs
    Run.ID = ['Solomon TSF1 fo ',num2str(fs(j)),' sin'];
    OutPutFolder = fullfile(pwd,Run.ID,[Run.ID,' - RC 00001.mat']);
    load(OutPutFolder);
    
    axes(ha(j+(j-1)));
    [VF,f] = Get_FS(VT,t);
    plot(f,abs(VF),'-r');
    grid on
    xlim([0,50]);
    if j~=Nfs,xticklabels(''); else, xlabel('f [Hz]'); end
    
    Run.ID = ['Solomon TSF1 fo ',num2str(fs(j)),' blair'];
    OutPutFolder = fullfile(pwd,Run.ID,[Run.ID,' - RC 00001.mat']);
    load(OutPutFolder);
    
    axes(ha(2*j));
    [VF,f] = Get_FS(VT,t);
    plot(f,abs(VF),'-r');
    grid on
    xlim([0,50]);
    if j~=Nfs,xticklabels(''); else, xlabel('f [Hz]'); end
    
end
Pos = hfig.Position;
Pos(2) = 20;
set(hfig,'Position',Pos);