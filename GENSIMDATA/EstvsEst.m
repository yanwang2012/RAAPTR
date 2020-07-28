% Script to plot different parameters over different estimated results.
% Yiqian 6/29/2020

clear;
tic

%% Dir settings
simParamsDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/searchParams/2bands/superNarrow';
simdataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands';
estSrc1Dir = '/Users/qyq/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands/Results_supNar';
estSrc2Dir = '/Users/qyq/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands/SupNar_xMBLT_iMBLT20/iMBLT20';
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
baseName = 'initial-iMBLT-after-xMBLT';
% SNR
figname1 = [baseName,'-SNR'];
snr1 = zeros(Nband,NestsrcBand);
snr2 = zeros(Nband,NestsrcBand);
for i = 1:Nband
    for j = 1:NestsrcBand
        [snr1(i,j),~] = Amp2Snr(EstSrc1{i,j},simParams,yr);
        [snr2(i,j),~] = Amp2Snr(EstSrc2{i,j},simParams,yr);
    end
end

figure(1)
subplot(1,2,1)
plot(snr1(1,:),snr2(1,:),'ob')
title('band 1')
xlabel('initial')
ylabel('iMBLT after xMBLT')
subplot(1,2,2)
plot(snr1(2,:),snr2(2,:),'or')
title('band 2')
xlabel('initial')
ylabel('iMBLT after xMBLT')
sgtitle('initial vs. iMBLT after xMBLT SNR')
saveas(gcf,[prefix,filesep,figname1],'png')
savefig([prefix,filesep,figname1])

% Frequency
figname2 = [baseName,'-Freq'];
freq1 = zeros(Nband,NestsrcBand);
freq2 = zeros(Nband,NestsrcBand);
for i = 1:Nband
    freq1(i,:) = arrayfun(@(x) EstSrc1{i,x}.omega/(2*pi*365*24*3600), 1:NestsrcBand);
    freq2(i,:) = arrayfun(@(x) EstSrc2{i,x}.omega/(2*pi*365*24*3600), 1:NestsrcBand);
end

figure(2)
subplot(1,2,1)
plot(freq1(1,:),freq2(1,:),'ob')
title('band 1')
xlabel('initial')
ylabel('iMBLT after xMBLT')
xlim([0 4.5e-8])
% ylim([0 4.5e-7])
subplot(1,2,2)
plot(freq1(2,:),freq2(2,:),'or')
title('band 2')
xlabel('initial')
ylabel('iMBLT after xMBLT')
sgtitle('initial vs. iMBLT after xMBLT Freq.')
saveas(gcf,[prefix,filesep,figname2],'png')
savefig([prefix,filesep,figname2])

% SNR vs. Freq.
figname3 = [baseName,'-SNRvsFreq'];
for i = 1:2
    figure(3+i-1)
    plot(snr1(i,:),freq1(i,:),'ob',snr2(i,:),freq2(i,:),'sr');
    text(snr1(i,:) + 2,freq1(i,:),num2str((1:numel(snr1(i,:)))'),'Color','#0072BD');
    text(snr2(i,:) - 2,freq2(i,:),num2str((1:numel(snr2(i,:)))'),'Color','#D95319','HorizontalAlignment','right');
    title([figname3,'Band ',num2str(i)]);
    xlabel('SNR');
    ylabel('Frequency');
    legend('initial','iMBLTafterxMBLT');
    saveas(gcf,[prefix,filesep,figname3,'Band ',num2str(i)],'png')
    savefig([prefix,filesep,figname3,'Band ',num2str(i)])
end

toc

%END