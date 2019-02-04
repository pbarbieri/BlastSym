clc
clear all
xi = 0.25;
XlsFile = 'SEQ02.xlsx';
Vw = 2500;
SiteX = 48.2;
SiteY = -145;
R = 0.8;


%% PPV 25
PPV = 25;

fo = 200;
OutputFileName = ['PPV',num2str(PPV,'%i'),' fo',num2str(fo,'%i'),' uniform'];
main(XlsFile,fo,xi,2,Vw,PPV/1000,SiteX,SiteY,R,OutputFileName);

OutputFileName = ['PPV',num2str(PPV,'%i'),' fo',num2str(fo,'%i'),' attenuated'];
main_r2(XlsFile,fo,xi,2,Vw,PPV/1000,SiteX,SiteY,R,OutputFileName);

fo = 800;

OutputFileName = ['PPV',num2str(PPV,'%i'),' fo',num2str(fo,'%i'),' uniform'];
main(XlsFile,fo,xi,2,Vw,PPV/1000,SiteX,SiteY,R,OutputFileName);

OutputFileName = ['PPV',num2str(PPV,'%i'),' fo',num2str(fo,'%i'),' attenuated'];
main_r2(XlsFile,fo,xi,2,Vw,PPV/1000,SiteX,SiteY,R,OutputFileName);

%% PPV 50
PPV = 50;

fo = 200;
OutputFileName = ['PPV',num2str(PPV,'%i'),' fo',num2str(fo,'%i'),' uniform'];
main(XlsFile,fo,xi,2,Vw,PPV/1000,SiteX,SiteY,R,OutputFileName);

OutputFileName = ['PPV',num2str(PPV,'%i'),' fo',num2str(fo,'%i'),' attenuated'];
main_r2(XlsFile,fo,xi,2,Vw,PPV/1000,SiteX,SiteY,R,OutputFileName);

fo = 800;

OutputFileName = ['PPV',num2str(PPV,'%i'),' fo',num2str(fo,'%i'),' uniform'];
main(XlsFile,fo,xi,2,Vw,PPV/1000,SiteX,SiteY,R,OutputFileName);

OutputFileName = ['PPV',num2str(PPV,'%i'),' fo',num2str(fo,'%i'),' attenuated'];
main_r2(XlsFile,fo,xi,2,Vw,PPV/1000,SiteX,SiteY,R,OutputFileName);

clear all