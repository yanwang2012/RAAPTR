clear;
tic
%% Extract parameters of sources in frequency bin X (Mauritius Poster)
% Load the frequency bin edges from the search parameter file for bin X.
simParamsDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/searchParams/2bands';
simDataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands';
estDataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands/MBLT/GWBsimDataSKASrlz1Nrlz3_MBLT1/Results20';
inputFileName = 'GWBsimDataSKASrlz1Nrlz3';
% Load the simulated source parameters.
load([simDataDir,filesep,inputFileName,'.mat'],'omega','alpha','delta',...
    'timingResiduals_tmp', 'yr','snr_chr','simParams');

%% setting fig axis
ybin_up = [];
ybin_low = [];
x = [];
y = [];

sy = [];
sx = [];

% Esti. Sky location
ra = [];
dec = [];
% simulated sky location
sra = [];
sdec = [];

%% reading the files
inParamsList = dir([simParamsDir,filesep,'searchParams','*.mat']);
inDataList = dir([estDataDir,filesep,'*',inputFileName,'*.mat']);

%%%%%%%%%%%%%%%%%%% DON'T FORGET %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nFile = dir([estDataDir,filesep,inputFileName,'band1','*.mat']); % count how many iterations are used.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_ite = length(nFile);
% inParamNames = {};
% inDataNames = {};
N = length(inParamsList);% number of bands
% % get the name and sort it
% for n = 1:N
%     inParamNames = [inParamNames inParamsList(n).name];
% end
%
% for m = 1:N*num_ite
%     inDataNames = [inDataNames inDataList(m).name];
% end

inParamNames = sort_nat({inParamsList.name});
inDataNames = sort_nat({inDataList.name});

%% data pre-processing
bx = zeros(N,50);
by = zeros(N,50);
etyband = 0; % band doesn't have sources inside.
for i = 1:N
    % load bands and estimated data
    load([simParamsDir,filesep,char(inParamNames(i))]);
    ybin_up = [ybin_up searchParams.angular_velocity(1)];% saving the band boundary
    ybin_low = [ybin_low searchParams.angular_velocity(2)];% Find the sources with frequencies in specific band
    Indx = find(omega >= searchParams.angular_velocity(2) & ...
        omega <= searchParams.angular_velocity(1));
    
    
    if isempty(Indx)
        disp(["There's no signal injected in band",num2str(i)]);
        etyband = i;
        continue
    end
    % Get their frequencies in Hz
    binsrcOmega = omega(Indx);
    y = [y binsrcOmega];
    size = length(binsrcOmega);
    by(i,1:size) = binsrcOmega;
    %binsrcF = (binsrcOmega/(2*pi))/(24*3600*365);
    % Get their SNR
    binsrcSNR = snr_chr(Indx);
    x = [x binsrcSNR];
    bx(i,1:size) = binsrcSNR;
    % Get the sky location
    sra = [sra alpha(Indx)];
    sdec = [sdec delta(Indx)];
    for j = 1:num_ite
        load([estDataDir,filesep,char(inDataNames(j + num_ite * (i-1)))],'bestRealLoc');
        disp(['File: ',char(inDataNames(j + num_ite * (i-1))),' loaded']);
        path_to_estimatedData = [estDataDir,filesep,char(inDataNames(j + num_ite * (i-1)))];
        
        %% Estimated source
        %path_to_simulationData = '~/Research/PulsarTiming/SimDATA/Mauritius/GWBsimDataSKA/GWBsimDataSKASrlz1Nrlz9.mat';
        %     path_to_pulsar_catalog = 'survey_ska.mat';
        estFreq = bestRealLoc(3)/(2*pi*365*24*3600);% convert unit from yr-1 to Hz
        [sourceParams]=ColSrcParams(path_to_estimatedData);
        [estSNR,estTimRes] = Amp2Snr(sourceParams,simParams,yr);
        sy = [sy estFreq];
        sx = [sx estSNR];
        ra = [ra bestRealLoc(1)];
        dec = [dec bestRealLoc(2)];
        %%
        % Plot the FFT of the timing residuals of the sources in bin 4
        %         timingResFFt = fft(timingResiduals_tmp');
        %         timingResPdg = timingResFFt(1:(floor(130/2)+1),:);
        %         timingResPdg = abs(timingResPdg);
        %         freqVec = (0:(floor(130/2)))*(1/((yr(end)-yr(1))*365*24*3600));
        % figure(1)
        % plot(freqVec,timingResPdg);
        % hold on
        % plot(bin5srcF,6e-6,'.')
        %%
        % Reproduce the Mauritius poster figure (Source frequency vs SNR)
        %     figure(1)
        %     plot(binsrcSNR,binsrcF,'o',estSNR,estFreq,'r*');
        %     xlabel('Network Signal to Noise Ratio');
        %     ylabel('Source Frequency (Hz)');
        %     legend('Injected source','Estimated source');
        %
        %     title(['Realization#9 bin',num2str(i)]);
        %     figname = ['Realization#9 bin',num2str(i)];
        %saveas(gcf,figname,'png')
        % Spacing of source frequencies relative to Fourier frequency spacing
        % [bin5srcFsort,bin5srcFsortIndx] = sort(bin5srcF,'descend');
        % figure(3)
        % plot(bin5srcSNR(bin5srcFsortIndx),[0,diff(sort(bin5srcFsort))/(freqVec(2)-freqVec(1))],'o');
        % xlabel('Network Signal to Noise Ratio');
        % ylabel(' Frequency spacing\times Data duration');
    end
end

y = y/(2*pi*365*24*3600);
uplim = max(max(x),max(sx))+50;
binSNR = 0:1:uplim;
binSNR_log = logspace(-5,3,length(binSNR));
ybin_up = ybin_up/(2*pi*365*24*3600);
ybin_up = repmat(ybin_up,length(binSNR),1);% stack itself vertically to broadcast to the dimension of x
ybin_low = ybin_low/(2*pi*365*24*3600);
ybin_low = repmat(ybin_low,length(binSNR),1);
stage = 1:1:num_ite;

%% Noise processing
% noisedir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test10/noise/Results/';
% simNoiseDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test10/noise';
% disp('Processing noise-only data');
% [avgnoise] = noiseprocess(noisedir,simNoiseDir,10,5);
% noisefile = dir([noisedir,filesep,'*.mat']);
% numNoise = length(noisefile); % number of noise files.
% simNoiseFileName = 'noise1.mat';
% noise = load([simNoiseDir,filesep,simNoiseFileName]);
% noiseFileName = {};
% nx = []; % blank array for noise SNR.

% for i = 1:numNoise
%     noiseFileName = [noiseFileName noisefile(i).name];
% end
%
% for i = 1:numNoise
%     path_noise_file = [noisedir,filesep,char(noiseFileName(i))];
%     noiseParams = ColSrcParams(path_noise_file);
%     [noiseSNR,~] = Amp2Snr(noiseParams,noise.simParams,yr);
%     nx = [nx noiseSNR];
% end
%
% avgnx = sum(reshape(nx,5,10),1)/N;

%% plot the entire map
close all;
figname = '2 bands 200 SRC 20 ITE 1 MBLT';
figure(1)
% yyaxis right
% loglog(x,y,'o',sx,sy,'kd','MarkerSize',10);
plot(x,y,'o',sx,sy,'s')
% semilogx(x,y,'o',sx,sy,'s');
% disp(sy)

hold on
% plot grid
for k=1:N
    %         semilogx(binSNR_log,ybin_up(:,k),'b-');
    plot(binSNR,ybin_up(:,k),'b-');
    %         semilogx(binSNR_log,ybin_low(:,k),'b--');
    plot(binSNR,ybin_low(:,k),'b--');
end
hold off

xlabel('SNR');
xlim([0 uplim]);
ylabel('Frequency');
legend('True','Estimated','Location','northeast');
title(figname);
saveas(gcf,figname,'png');
savefig(figname);
% save('estTimRes01.mat','estTimRes');

% figure(2)
% Legend = {N,1};
% hold on
% for i = 1:N - etyband
%     %plot(stage,sx(10*(i-1)+1:10*i));
%     semilogy(stage,sx(num_ite*(i-1)+1:num_ite*i));
%     Legend{i} = ['Band ', num2str(i)+etyband];
%     if i == (N - etyband)
%         %plot(stage, avgnx,'--k');
%         semilogy(stage,avgnoise,'--k');
%         Legend{i+1} = 'Noise';
%     end
% end
% hold off
% legend(Legend);
% title([figname,' SNR vs. Stage']);
% xlabel('Stage');
% ylabel('SNR');
% saveas(gcf,[figname,' SNR-Stage'],'png');
% savefig([figname,' SNR-Stage']);

figure(3)
plot(sra,sdec,'ob',ra,dec,'sr');
title(['Sky location for ',figname]);
legend('Simulated source','Estimated source');
xlabel('RA');
ylabel('Dec');
saveas(gcf,[figname,' skyloc'],'png')
savefig([figname,' skyloc']);

figure(4)
plot(sra,y,'ob',ra,sy,'sr');
title(['RA vs. Freq for ',figname]);
xlabel('RA');
ylabel('Freq.');
legend('True','Estimated')
saveas(gcf,[figname,' RA'],'png')
savefig([figname,' RA']);


figure(5)
plot(sdec,y,'ob',dec,sy,'sr');
title(['DEC vs. Freq for ',figname]);
xlabel('DEC');
ylabel('Freq.');
legend('True','Estimated')
saveas(gcf,[figname,' DEC'],'png')
savefig([figname,' DEC']);

%% set up cutoff
SNRcut = 50;
Idx = find(x >= SNRcut); % set SNR cutoff for simulated sources
Sy = y(Idx);
sra = sra(Idx);
sdec = sdec(Idx);

Sidx = find(sx >= SNRcut);% SNR cutoff for est. sources
sx = sx(Sidx);
sy = sy(Sidx);
ra = ra(Sidx);
dec = dec(Sidx);

figure(6)
plot(sra,Sy,'ob',ra,sy,'sr');
title(['RA vs. Freq for ',figname]);
xlabel('RA');
ylabel('Freq.');
legend('True','Estimated')
saveas(gcf,[figname,' RA-SNRcutoff ',num2str(SNRcut)],'png')
savefig([figname,' RA-SNRcutoff ',num2str(SNRcut)]);


figure(7)
plot(sdec,Sy,'ob',dec,sy,'sr');
title(['DEC vs. Freq for ',figname]);
xlabel('DEC');
ylabel('Freq.');
legend('True','Estimated')
saveas(gcf,[figname,' DEC-SNRcutoff ',num2str(SNRcut)],'png')
savefig([figname,' DEC-SNRcutoff ',num2str(SNRcut)]);

figure(8)
% yyaxis right
% loglog(x,y,'o',sx,sy,'kd','MarkerSize',10);
plot(x,y,'o',sx,sy,'s')
% semilogx(x,y,'o',sx,sy,'s');
% disp(sy)

hold on
% plot grid
for k=1:N
    %         semilogx(binSNR_log,ybin_up(:,k),'b-');
    plot(binSNR,ybin_up(:,k),'b-');
    %         semilogx(binSNR_log,ybin_low(:,k),'b--');
    plot(binSNR,ybin_low(:,k),'b--');
end
hold off

xlabel('SNR');
xlim([0 uplim]);
ylabel('Frequency');
legend('True','Estimated','Location','northeast');
title([figname,' SNRcutoff ',num2str(SNRcut)]);
saveas(gcf,[figname,' SNRcutoff ',num2str(SNRcut)],'png');
savefig([figname,' SNRcutoff ',num2str(SNRcut)]);

%% freq. only plot
yb = Sy(1e-7 < Sy & Sy <= 2e-7);% simulated sources blow the boundary, 2e-7 is the freq. separates two bands
yb2 = Sy(Sy <= 1e-7);
szyb = length(yb);
szyb2 = length(yb2);
cstb = zeros(szyb,1);
cstb2 = zeros(szyb2,1);
yu = Sy(Sy > 2e-7); 
szyu = length(yu);
cstu = zeros(szyu,1);

syb = sy(1e-7< sy & sy <= 2e-7); % est. sources...
syb2 = sy(sy <= 1e-7);
szsyb = length(syb);
szsyb2 = length(syb2);
scstb = zeros(szsyb,1);
scstb2 = zeros(szsyb2,1);
syu = sy(sy > 2e-7);
szsyu = length(syu);
scstu = zeros(szsyu,1);



figure(9)
subplot(3,1,2)
plot(yb,cstb,'ob',syb,scstb,'sr');
xlabel('Frequency');
ylabel('Constant');
title([figname,' Freq-only SNRcutoff ',num2str(SNRcut),' Band 1 upper']);
legend('Ture','Estimated');
subplot(3,1,3)
plot(yu,cstu,'ob',syu,scstu,'sr')
xlabel('Frequency');
ylabel('Constant');
title([figname,' Freq-only SNRcutoff ',num2str(SNRcut),' Band 2']);
legend('Ture','Estimated');
subplot(3,1,1)
plot(yb2,cstb2,'ob',syb2,scstb2,'sr')
xlabel('Frequency');
ylabel('Constant');
title([figname,' Freq-only SNRcutoff ',num2str(SNRcut),' Band 1 lower']);
legend('Ture','Estimated');
saveas(gcf,[figname,' Freq-only SNRcutoff ',num2str(SNRcut)],'png')
savefig([figname,' Freq-only SNRcutoff ',num2str(SNRcut)]);

toc

%END