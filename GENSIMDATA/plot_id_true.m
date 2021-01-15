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
simRA = alpha;
simDec = delta;

idRA = [];
idDec = [];

idBand1 = sum(~cellfun('isempty',idsrc(1,:)));
idBand2 = sum(~cellfun('isempty',idsrc(2,:)));

for b = 1:2 % for 2 bands
    switch b
        case 1
            N = idBand1;
        case 2
            N = idBand2;
    end
    for i = 1:N
        idRA = [idRA idsrc{b,i}.alpha];
        idDec = [idDec idsrc{b,i}.delta];
    end
end

%% plot
figure
plot(simRA,simDec,'bo') % plot true sources
hold on
plot(idRA,idDec,'rs') % plot identified sources
plot(matched_alpha,matched_dec,'g*') % plot truly matched true sources

for j = 1:length(idRA)
    plot([idRA(j),matched_alpha(j)],[idDec(j),matched_dec(j)],'Color','m') % connect identified and matched sources
end

xlabel('RA')
ylabel('DEC')
title('Sky Location')
legend('True Sources', 'Identified Sources','Matched True Source','Line bewteen identified and turly matched source','Location','bestoutside')
saveas(gcf,[idDataDir,filesep,'fig',filesep,'SkyLocation.png'])
savefig([idDataDir,filesep,'fig',filesep,'SkyLocation'])
