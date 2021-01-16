% Script to plot frequency gaps between estimated sources in order to place
% the band edge more properly.

% Author: QYQ
% 1/7/2021
%% Data Dir
clear

estSrcDir = '/Users/qyq/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands/SuperNarrow/Results_supNar_rand1/GWBsimDataSKASrlz1Nrlz3_xMBLT/results/Union2';
Filename = 'GWBsimDataSKASrlz1Nrlz3';
ext = '.mat';

estSrcFile = dir([estSrcDir,filesep,'*',Filename,'*',ext]);
estSrcFilename = sort_nat({estSrcFile.name});
Nfiles = length(estSrcFile);

%% load estimated sources
Np = 1000; % # of pulsars used in simulation
estSrc = {};
estFreq = zeros(1,Nfiles);

for i = 1:Nfiles
    path_to_estSrc = [estSrcDir,filesep,estSrcFilename{i}];
    estSrc{i} = ColSrcParams(path_to_estSrc,Np);
    estFreq(i) = estSrc{i}.omega/(2*pi*365*24*3600); % convert rad/year -> Hz
end

estFreq = sort(estFreq); % sort frequencies in ascending order

% calculate joint freq difference
difreq = zeros(1,Nfiles);
for j = 1:Nfiles - 1
        difreq(j) = abs(estFreq(j+1) - estFreq(j))/estFreq(j); % relative frequency gap
end

%% plot
figure
plot(estFreq,difreq,'ro')
xlabel('Frequency [Hz]')
ylabel('Relative Freq. Diff')
title('Frequency Gap')
saveas(gcf,[estSrcDir,filesep,'fig',filesep,'FreqGap.png'])
savefig([estSrcDir,filesep,'fig',filesep,'FreqGap'])

%EOF