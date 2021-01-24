% script to plot different parameters between identified sources and True sources.

% Author: QYQ
% 1/7/2021
clear;
%% Load data
simDataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11';
idDataDir = '/Users/qyq/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands/SuperNarrow/Union2_xMBLT/Union2_iMBLT_after20';
Filename = 'GWBsimDataSKASrlz1Nrlz3';
identifyFilename = 'identifiedSrc';
ext = '.mat';

simFile = [simDataDir,filesep,Filename,ext];
idFile = [idDataDir,filesep,identifyFilename,ext];
id2true = [idDataDir,filesep,'Matched_Sources.mat'];
load(simFile);
load(idFile);
load(id2true);

%% Sky location
simRA_nm = []; % not matched true sources
simDec_nm = [];

idRA = [];
idDec = [];
idSNR = [];

idBand1 = sum(~cellfun('isempty',idsrc(1,:)));
idBand2 = sum(~cellfun('isempty',idsrc(2,:)));

for b = 1:2 % for 2 bands
    switch b
        case 1
            N = idBand1;
        case 2
            N = idBand2;
    end
    idx = setdiff(1:length(SrcAlpha{b}),id_max(:,b)); % get the complementary index
    for i = 1:N
        idRA = [idRA idsrc{b,i}.alpha];
        idDec = [idDec idsrc{b,i}.delta];
        [idSNR_tmp,~] = Amp2Snr(idsrc{b,i},simParams,yr);
        idSNR = [idSNR idSNR_tmp];
        simRA_nm = [simRA_nm SrcAlpha{b}(idx)];
        simDec_nm = [simDec_nm SrcDelta{b}(idx)];
    end
end

%% plot
figure
% plot(simRA,simDec,'o','MarkerEdgeColor','##D4D8D9') % plot true sources
scatter(simRA_nm,simDec_nm,'o','MarkerEdgeColor','#D4D8D9') % use scatter to make marker size change

hold on

% plot(idRA,idDec,'rs') % plot identified sources
% scatter(idRA,idDec,idSNR,'rs') % with different size
scatter(idRA,idDec,[],idSNR, 'filled') % with different color

% plot(matched_alpha,matched_dec,'ob') % plot truly matched true sources
% scatter(matched_alpha,matched_dec,matched_snr,'ob')
scatter(matched_alpha,matched_dec,[],matched_snr,'s')

for j = 1:length(idRA)
    plot([idRA(j),matched_alpha(j)],[idDec(j),matched_dec(j)],'Color','m') % connect identified and matched sources
end

xlabel('RA')
ylabel('DEC')
title('Sky Location')
colormap autumn
c = colorbar;
ylabel(c,'SNR','FontSize',14)
legend('True Sources', 'Identified Sources','Matched True Source','Matched & Identi.','Location','bestoutside')
saveas(gcf,[idDataDir,filesep,'fig',filesep,'SkyLocationC.png'])
savefig([idDataDir,filesep,'fig',filesep,'SkyLocationC'])

%END