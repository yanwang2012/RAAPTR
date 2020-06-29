% Script to plot different parameters over different estimated results.
% Yiqian 6.29 2020

clear;
tic

%% Dir settings
simParamsDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/searchParams/2bands/superNarrow';
simdataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands';
estSrc1Dir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands/Results_supNar/GWBsimDataSKASrlz1Nrlz3_xMBLT/results';
estSrc2Dir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands/supperNarrow_iMBLT20';
Filename = 'GWBsimDataSKASrlz1Nrlz3';
ext = '.mat';

%% Files
paraFile = dir([simParamsDir,filesep,'searchParams','*.mat']);
Nband = length(paraFile);
simFile = [simdataDir,filesep,Filename,ext];
estSrc1File = dir([estSrc1Dir,filesep,'*',Filename,'*',ext]);
estSrc2File = dir([estSrc2Dir,filesep,'*',Filename,'*',ext]);
Nestsrc = length(estSrc2File);
paraFilename = sort_nat({paraFile.name});
estSrc2Filename = sort_nat({estSrc2File.name});
estSrc1Filename = sort_nat({estSrc1File.name});
load(simFile,'simParams','yr');

%% Get estimated sources info
NestsrcBand = Nestsrc/Nband; % number of sources in a band.
EstSrc2 = {};
EstSrc1 = {};
for band = 1:Nband
    for k = 1:NestsrcBand
        path_to_estimatedDataestSrc2 = [estSrc2Dir,filesep,char(estSrc2Filename((band - 1) * NestsrcBand + k))];
        path_to_estimatedDataestSrc1 = [estSrc1Dir,filesep,char(estSrc1Filename((band - 1) * NestsrcBand + k))];
        EstSrc2{band,k} = ColSrcParams(path_to_estimatedDataestSrc2);
        EstSrc1{band,k} = ColSrcParams(path_to_estimatedDataestSrc1);
    end
end

%% Plotting
prefix = [estSrc2Dir,filesep,'EstvsEst'];
mkdir(prefix)
% SNR
figname1 = 'xMBLT-iMBLT-SNR';
snr1 = zeros(Nband,NestsrcBand);
snr2 = zeros(Nband,NestsrcBand);
for i = 1:Nband
    for j = 1:NestsrcBand
        [snr1(i,j),~] = Amp2Snr(EstSrc1{i,j},simParams,yr);
        [snr2(i,j),~] = Amp2Snr(EstSrc2{i,j},simParams,yr);
    end
end

subplot(1,2,1)
plot(snr1(1,:),snr2(1,:),'ob')
title('band 1')
xlabel('xMBLT')
ylabel('iMBLT')
subplot(1,2,2)
plot(snr1(2,:),snr2(2,:),'or')
title('band 2')
xlabel('xMBLT')
ylabel('iMBLT')
sgtitle('xMBLT vs. iMBLT SNR')
saveas(gcf,[prefix,filesep,figname1],'png')
savefig([prefix,filesep,figname1])

% Frequency
figname2 = 'xMBLT-iMBLT-Freq';
freq1 = zeros(Nband,NestsrcBand);
freq2 = zeros(Nband,NestsrcBand);
for i = 1:Nband
    freq1(i,:) = arrayfun(@(x) EstSrc1{i,x}.omega/(2*pi*365*24*3600), 1:NestsrcBand);
    freq2(i,:) = arrayfun(@(x) EstSrc2{i,x}.omega/(2*pi*365*24*3600), 1:NestsrcBand);
end

subplot(1,2,1)
plot(freq1(1,:),freq2(1,:),'ob')
title('band 1')
xlabel('xMBLT')
ylabel('iMBLT')
subplot(1,2,2)
plot(freq1(2,:),freq2(2,:),'or')
title('band 2')
xlabel('xMBLT')
ylabel('iMBLT')
sgtitle('xMBLT vs. iMBLT Freq.')
saveas(gcf,[prefix,filesep,figname2],'png')
savefig([prefix,filesep,figname2])




toc

%END