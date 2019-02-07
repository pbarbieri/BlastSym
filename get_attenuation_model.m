function [Model] = get_attenuation_model(XlsFile,PlotFlag)
%% ========================================================================
% Copyright SRK/FIUBA (C) 2018
% Coded By: P. Barbieri (pbarbieri@fi.uba.ar)
% Version: V103
% Date: 21-01-2019
% -------------------------------------------------------------------------
% USAGE: 
%   get_attenuation_model: computes an attenuation model from PPV
%                          measurmentes
%   
% -------------------------------------------------------------------------
% INPUT:
%   XlsFile     Excel file with R, W and PPV data
%   PlotFlag    Flag, 1 for ploting.
%
% -------------------------------------------------------------------------
% OUTPUT:
%       
% -------------------------------------------------------------------------
% EXAMPLE:
%   
% -------------------------------------------------------------------------
% BIBLIO:
%   
% -------------------------------------------------------------------------
% VALIDATE:
% Version:
% Date:
% Validated by:
% -------------------------------------------------------------------------
% LOG
%   V101    19/01/2019      First version
% =========================================================================

[~,~,PPVData] =  xlsread(XlsFile);
PPVData = cell2struct(PPVData(2:end,:),PPVData(1,:),2);
for k = 1:size(PPVData,1)
    PPVData(k).R = PPVData(k).D/sqrt(PPVData(k).MIC);
end


X = [PPVData.R].';
Y = [PPVData.PPV].';
f=fit(log10(X),log10(Y),'poly1');
Int = confint(f);
a = f.p1;
b = log(10.^f.p2);

Nsigma = icdf(makedist('normal',0,1),1-0.025);
muLnPPV = @(r) a*log(r)+b;
sigmaLnPPV = (log(10.^Int(2,2))-log(10.^Int(1,2)))/Nsigma;

PPVdist = @(r) makedist('lognormal',muLnPPV(r),sigmaLnPPV);

Model.a = a;
Model.b = b;
Model.muLnPPV = muLnPPV;
Model.sigmaLnPPV = sigmaLnPPV;
Model.PPVdist = PPVdist;
Model.fit = f;


if PlotFlag
    X = linspace(min([PPVData.R]),max([PPVData.R]),1000).';
    Z = zeros(1000,3);
    for k = 1:1000
        Z(k,1) = icdf(PPVdist(X(k)),0.5);
        Z(k,2) = icdf(PPVdist(X(k)),0.975);
        Z(k,3) = icdf(PPVdist(X(k)),0.025);
    end
    close all
    hfig = figure(1);
    set(hfig,'Color',[1 1 1]);
    hp = gobjects(4,1);
    hold on
    hp(1) = scatter([PPVData.R],[PPVData.PPV],15,'filled','blue');
    hp(2) = plot(X,X.^f.p1*10^f.p2,'-r','linewidth',2);
    hp(3) = plot(X,X.^Int(1,1)*10^Int(1,2),'--r','linewidth',2);
    plot(X,X.^Int(2,1)*10^Int(2,2),'--r','linewidth',2);
    plot(X,Z(:,1),'--k','linewidth',2);
    hp(4) = plot(X,Z(:,2),'-.k','linewidth',2);
    plot(X,Z(:,3),'-.k','linewidth',2);
    hold off
    grid on
    xlim([5,120]);
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    legend(hp,{'Measured data','Median fit','95% confidence fit','Proposed model'});
    xlabel('R/M^{-0.5} [m/kg^{-0.5}]');
    ylabel('PPV [mm/s]');
end


end