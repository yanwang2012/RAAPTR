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
load(simFile);
load(idFile);

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
plot(simRA,simDec,'bo',idRA,idDec,'rs');
xlabel('RA')
ylabel('DEC')
title('Sky Location')
legend('True Sources', 'Identified Sources','Location','bestoutside')
saveas(gcf,[idDataDir,filesep,'fig',filesep,'SkyLocation.png'])
savefig([idDataDir,filesep,'fig',filesep,'SkyLocation'])
