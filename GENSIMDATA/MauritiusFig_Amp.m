clear;
tic
%% Extract parameters of sources in frequency bin X (Mauritius Poster)
% Load the frequency bin edges from the search parameter file for bin X.
simParamsDir = '~/Research/PulsarTiming/SimDATA/Mauritius/searchParams_GWBsimDataSKA';
simDataDir = '~/Research/PulsarTiming/SimDATA/Mauritius/GWBsimDataSKA';
estDataDir = '~/Research/PulsarTiming/SimDATA/Mauritius/MauritiusNew_results_mat';
% Load the source parameters across the entire frequency range
load([simDataDir,filesep,'GWBsimDataSKASrlz1Nrlz9.mat'],'omega',...
    'timingResiduals_tmp','yr','Amp');

%% setting fig axis
y = [];
ybin_up = [];
ybin_low = [];
x = [];

sy = [];
sx = [];

%% reading the files and sort it naturely
inParamsList = dir([simParamsDir,filesep,'*.mat']);
inDataList = dir([estDataDir,filesep,'*Nrlz9.mat']);
inParamNames = {};
inDataNames = {};
% get the name and sort it
for n = 1:10
    inParamNames = [inParamNames inParamsList(n).name];
    inDataNames = [inDataNames inDataList(n).name];
end
inParamNames = sort_nat(inParamNames);
inDataNames = sort_nat(inDataNames);

%% data pre-processing
for i = 1:length(inDataList)
    % load bands and estimated source
    load([simParamsDir,filesep,char(inParamNames(i))]);
    load([estDataDir,filesep,char(inDataNames(i))],'bestRealLoc');
    %path_to_estimatedData = [estDataDir,filesep,char(inDataNames(i))];
    
    % Find the sources with frequencies in bin 4
    binsrcOmgIndx = find(omega >= xmaxmin(3,2) & omega <= xmaxmin(3,1));
    ybin_up = [ybin_up xmaxmin(3,1)];% saving the band boundary
    ybin_low = [ybin_low xmaxmin(3,2)];
    
    if isempty(binsrcOmgIndx)
        disp(["There's no signal injected in band",num2str(i)]);
        continue
    end
    % Get their frequencies in Hz
    binsrcOmega = omega(binsrcOmgIndx);
    y = [y binsrcOmega];
    %binsrcF = (binsrcOmega/(2*pi))/(24*3600*365);
    % Get their SNR
    %binsrcSNR = snr_chr(binsrcOmgIndx);
    % Get their amplitude
    binsrcAmp = Amp(binsrcOmgIndx);
    x = [x binsrcAmp];
    %% Estimated source
    %path_to_simulationData = '~/Research/PulsarTiming/SimDATA/Mauritius/GWBsimDataSKA/GWBsimDataSKASrlz1Nrlz9.mat';
    %path_to_pulsar_catalog = 'survey_ska.mat';
    estFreq = bestRealLoc(3)/(2*pi*365*24*3600);
    %estSNR = convertAmp2snr(path_to_estimatedData,path_to_pulsar_catalog);
    estAmp = bestRealLoc(5);
    sy = [sy estFreq];
    sx = [sx estAmp];
    %%
    % Plot the FFT of the timing residuals of the sources in bin 4
    timingResFFt = fft(timingResiduals_tmp');
    timingResPdg = timingResFFt(1:(floor(130/2)+1),:);
    timingResPdg = abs(timingResPdg);
    freqVec = (0:(floor(130/2)))*(1/((yr(end)-yr(1))*365*24*3600));
    % figure(1)
    % plot(freqVec,timingResPdg);
    % hold on
    % plot(bin5srcF,6e-6,'.')
    %%
    % Reproduce the Mauritius poster figure (Source frequency vs Amp)
%     figure(1)
%     plot(binAmp,binsrcF,'o',estAmp,estFreq,'r*');
%     xlabel('Network Signal to Noise Ratio');
%     ylabel('Source Frequency (Hz)');
%     legend('Injected source','Estimated source');
%     if i == 1
%         title(['Realization#9 bin',num2str(1)]);
%         figname = ['Realization#9 bin',num2str(1)];
%     elseif i == 2
%         title(['Realization#9 bin',num2str(10)]);
%         figname = ['Realization#9 bin',num2str(10)];
%     else
%         title(['Realization#9 bin',num2str(i-1)]);
%         figname = ['Realization#9 bin',num2str(i-1)];
%     end
    %pause
    %saveas(gcf,figname,'png')
    % Spacing of source frequencies relative to Fourier frequency spacing
    % [bin5srcFsort,bin5srcFsortIndx] = sort(bin5srcF,'descend');
    % figure(3)
    % plot(bin5srcSNR(bin5srcFsortIndx),[0,diff(sort(bin5srcFsort))/(freqVec(2)-freqVec(1))],'o');
    % xlabel('Network Signal to Noise Ratio');
    % ylabel(' Frequency spacing\times Data duration');
end
y = y/(2*pi*365*24*3600);
%binSNR = 0:1:450;
binAmp = 0:10^(-9):10^(-7);
ybin_up = ybin_up/(2*pi*365*24*3600);
ybin_up = repmat(ybin_up,length(binAmp),1);% repeat itself to broadcast to the dimension of x
ybin_low = ybin_low/(2*pi*365*24*3600);
ybin_low = repmat(ybin_low,length(binAmp),1);
%% plot the entire map
figure(2)
yyaxis left
semilogy(log10(x),y,'o',log10(sx),sy,'*','MarkerSize',10);
%hold on
% for j=1:10
%     semilogy(binAmp,ybin_up(:,j),'b-');
%     %plot(binSNR,ybin_up(:,j),'b-');
%     semilogy(binSNR,ybin_low(:,j),'b--');
%     %plot(binSNR,ybin_low(:,j),'b--');
% 
% end
xlabel('Amplitude');
ylabel('Frequency');
%saveas(gcf,'MauritiusFig_Amp','png');
toc