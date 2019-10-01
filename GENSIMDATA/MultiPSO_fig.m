clear;
tic
%% Extract parameters of sources in frequency bin X (Mauritius Poster)
% Load the frequency bin edges from the search parameter file for bin X.
simParamsDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test3/searchParams_Nyquist';
simDataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test6/Combo2';
estDataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test6/Combo2/Combo2/Results';
inputFileName = 'GWBsimDataSKASrlz1Nrlz3.mat';
% Load the source parameters across the entire frequency range
load([simDataDir,filesep,inputFileName],'omega',...
    'timingResiduals_tmp', 'yr','snr_chr','simParams');

%% setting fig axis
ybin_up = [];
ybin_low = [];
x = [];
y = [];

sy = [];
sx = [];

%% reading the files
inParamsList = dir([simParamsDir,filesep,'*.mat']);
inDataList = dir([estDataDir,filesep,'*Nrlz3*.mat']);
nFile = dir([estDataDir,filesep,'3_*.hdf5']); % count how many iterations are used.
num_ite = length(nFile);
inParamNames = {};
inDataNames = {};
N = length(inParamsList);% number of files
% get the name and sort it
for n = 1:N
    inParamNames = [inParamNames inParamsList(n).name];
end

% need to be improved, temporary use.
for m = 1:2*num_ite % 2* means 2 bands are used.
    inDataNames = [inDataNames inDataList(m).name]; 
end

inParamNames = sort_nat(inParamNames);
inDataNames = sort_nat(inDataNames);

%% data pre-processing
for i = 1:N
    % load bands and estimated data
    load([simParamsDir,filesep,char(inParamNames(i))]);
    ybin_up = [ybin_up searchParams.angular_velocity(1)];% saving the band boundary
    ybin_low = [ybin_low searchParams.angular_velocity(2)];
    for j = 1:2*num_ite % need to be imporved.
        load([estDataDir,filesep,char(inDataNames(j))],'bestRealLoc');
        path_to_estimatedData = [estDataDir,filesep,char(inDataNames(j))];
        
        % Find the sources with frequencies in bin 4
        binsrcOmgIndx = find(omega >= searchParams.angular_velocity(2) & ...
            omega <= searchParams.angular_velocity(1));
        
        
        if isempty(binsrcOmgIndx)
            disp(["There's no signal injected in band",num2str(i)]);
            continue
        end
        % Get their frequencies in Hz
        binsrcOmega = omega(binsrcOmgIndx);
        y = [y binsrcOmega];
        %binsrcF = (binsrcOmega/(2*pi))/(24*3600*365);
        % Get their SNR
        binsrcSNR = snr_chr(binsrcOmgIndx);
        x = [x binsrcSNR];
        %% Estimated source
        %path_to_simulationData = '~/Research/PulsarTiming/SimDATA/Mauritius/GWBsimDataSKA/GWBsimDataSKASrlz1Nrlz9.mat';
        %     path_to_pulsar_catalog = 'survey_ska.mat';
        phiI = bestRealLoc(8:1007);% esimated pulsar phases
        estFreq = bestRealLoc(3)/(2*pi*365*24*3600);% convert unit from yr-1 to Hz
        [sourceParams]=ColSrcParams(path_to_estimatedData);
        %     [pulsarParams]=ColPsrParams(path_to_pulsar_catalog);
        %disp(sourceParams.Amp)
        [estSNR,estTimRes] = Amp2Snr(sourceParams,simParams,phiI,yr);
        sy = [sy estFreq];
        sx = [sx estSNR];
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
binSNR = 0:1:450;
ybin_up = ybin_up/(2*pi*365*24*3600);
ybin_up = repmat(ybin_up,length(binSNR),1);% stack itself vertically to broadcast to the dimension of x
ybin_low = ybin_low/(2*pi*365*24*3600);
ybin_low = repmat(ybin_low,length(binSNR),1);
stage = 1:1:num_ite; % number of subtractions.

%% plot the entire map
figname = 'b2b3';
figure(1)
%yyaxis right
%loglog(x,y,'o',sx,sy,'kd','MarkerSize',10);
plot(x,y,'o',sx,sy,'s')
% disp(sy)
% plot the grid
hold on
for k=1:N
    %semilogy(binSNR,ybin_up(:,j),'b-');
    plot(binSNR,ybin_up(:,k),'b-');
    %semilogy(binSNR,ybin_low(:,j),'b--');
    plot(binSNR,ybin_low(:,k),'b--');
end
hold off

xlabel('SNR');
ylabel('Frequency');
legend('True','Estimated','Location','northeast');
title([figname,' combo']);
saveas(gcf,figname,'png');
savefig(figname);
%save('estTimRes01.mat','estTimRes');

figure(2)
plot(stage,sx(1:10),'b',stage,sx(11:20),'r');
title([figname,' SNR vs. Stage']);
xlabel('Stage');
ylabel('SNR');
legend('band 2','band 3');
saveas(gcf,[figname,'_stage'],'png');
savefig([figname,'_stage']);

toc