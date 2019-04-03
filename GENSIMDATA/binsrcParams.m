clear;
%% Extract parameters of sources in frequency bin X (Mauritius Poster)
% Load the frequency bin edges from the search parameter file for bin X.
load('~/Research/PulsarTiming/SimDATA/Mauritius/searchParams_GWBsimDataSKA/searchParams_GWBsimDataSKA2.mat');
%%
% Load the source parameters across the entire frequency range
load('~/Research/PulsarTiming/SimDATA/Mauritius/GWBsimDataSKA/GWBsimDataSKASrlz1Nrlz9.mat',...
     'omega', 'timingResiduals_tmp', 'yr','snr_chr');
% load the estimated source 
load('~/Research/PulsarTiming/SimDATA/Mauritius/Mauritius_results_mat/3_GWBsimDataSKASrlz1Nrlz9.mat',...
    'bestRealLoc');
%%
% Find the sources with frequencies in bin 4
bin5srcOmgIndx = find(omega >=xmaxmin(3,2) & omega <=xmaxmin(3,1));
% Get their frequencies in Hz
bin5srcOmega = omega(bin5srcOmgIndx);
bin5srcF = (bin5srcOmega/(2*pi))/(24*3600*365);
% Get their SNR
bin5srcSNR = snr_chr(bin5srcOmgIndx);

%% Estimated source
%path_to_simulationData = '~/Research/PulsarTiming/SimDATA/Mauritius/GWBsimDataSKA/GWBsimDataSKASrlz1Nrlz9.mat';
path_to_estimatedData = '~/Research/PulsarTiming/SimDATA/Mauritius/Mauritius_results_mat/3_GWBsimDataSKASrlz1Nrlz9.mat';
path_to_pulsar_catalog = 'survey_ska.mat';
estFreq = bestRealLoc(3)/(2*pi*365*24*3600);
estSNR = convertAmp2snr(path_to_estimatedData,path_to_pulsar_catalog);
%%
% Plot the FFT of the timing residuals of the sources in bin 4
timingResFFt = fft(timingResiduals_tmp');
timingResPdg = timingResFFt(1:(floor(130/2)+1),:);
timingResPdg = abs(timingResPdg);
freqVec = (0:(floor(130/2)))*(1/((yr(end)-yr(1))*365*24*3600));
figure(1)
plot(freqVec,timingResPdg);
hold on
plot(bin5srcF,6e-6,'.')
%%
% Reproduce the Mauritius poster figure (Source frequency vs SNR) for bin 4
figure(2)
plot(bin5srcSNR,bin5srcF,'o',estSNR,estFreq,'r*');
xlabel('Network Signal to Noise Ratio');
ylabel('Source Frequency (Hz)');
legend('Injected source','Estimated source');
title('Realization 9 bin 10')
% Spacing of source frequencies relative to Fourier frequency spacing
%figure(3)
[bin5srcFsort,bin5srcFsortIndx] = sort(bin5srcF,'descend');
figure(4)
plot(bin5srcSNR(bin5srcFsortIndx),[0,diff(sort(bin5srcFsort))/(freqVec(2)-freqVec(1))],'o');
xlabel('Network Signal to Noise Ratio');
ylabel(' Frequency spacing\times Data duration');