% Script to plot frequency gaps between estimated sources in order to place
% the band edge more properly.

% Author: QYQ
% 1/7/2021
%% Data Dir
clear

estSrcDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands/SuperNarrow/Results_supNar_rand1/GWBsimDataSKASrlz1Nrlz3_xMBLT/results/Union2';
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

threshold = 0.4;
idx = find(difreq > 0.4 & difreq < 0.9);
bandedge = (1 + 1/2 * difreq(idx)) .* estFreq(idx);
bandedge_av = bandedge * 365*24*3600*2*pi; % convert Hz to rad/year

%% Generate search file
filename = 'Nyquist.json';
outDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/searchParams/Band_opt/New';
mkdir(outDir);
searchParams = jsondecode(fileread(filename));
NumBands = length(idx);%5; %10; % total number of bands
% bandwidth = bandedge_av;
FreqRange = searchParams.angular_velocity;
for i = 1:NumBands + 1 % i band edges can split whole band into i+1 segments
    if i == NumBands + 1
        searchParams.angular_velocity(1) = FreqRange(1);
        save([outDir,filesep,'searchParams_Nyquist',num2str(i),'.mat'],'searchParams','NumBands','FreqRange');
    else
        searchParams.angular_velocity(1) = bandedge_av(i); % update upper limit
        searchParams.band_num = i;
        save([outDir,filesep,'searchParams_Nyquist',num2str(i),'.mat'],'searchParams','NumBands','FreqRange');
        tmp = searchParams.angular_velocity(1); % save upper limit for prev. band
        searchParams.angular_velocity(2) = tmp; % update lower limit
    end
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