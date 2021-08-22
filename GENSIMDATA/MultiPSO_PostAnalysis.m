% Multiple simulation data version for MultiPSO_fig.m
clear;
tic
%% Extract parameters of sources in frequency bin X (Mauritius Poster)
% Load the frequency bin edges from the search parameter file for bin X.
simParamsDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/searchParams/Band_opt';
simDataDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/simData/Band_opt_diff';
estDataDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/results_diff_opt_iMBLT';
inputFileName = 'GWBsimDataSKASrlz*Nrlz1'; % change here to switch between multiple and single realizations
% Load the simulated source parameters.
simDataList = dir([simDataDir,filesep,inputFileName,'.mat']);
simDataList = sort_nat({simDataList.name});
rlzs = length(simDataList);

for rlz = 1:rlzs
    load([simDataDir,filesep,simDataList{rlz}],'omega','alpha','delta',...
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
    
    %%%%%%%%%%%%%%%%%%% DON'T FORGET %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [~,simFileName,~] = fileparts(simDataList{rlz});
    nFile = dir([estDataDir,filesep,'1_',simFileName,'*.mat']); % count how many iterations are used. For initial PSO est.
%     nFile = dir([estDataDir,filesep,simFileName,'*.mat']); % For MBLT and other tests with irregular filename.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% reading the files
    %     inParamsList = dir([simParamsDir,filesep,'searchParams','*.mat']);
    inParamsList = dir([simParamsDir,filesep,simFileName,filesep,'searchParams','*.mat']); % corresponds to each realization's own optimal searchParams
    inDataList = dir([estDataDir,filesep,'*',simFileName,'*.mat']);
    
    num_ite = length(nFile);
    
    inParamNames = sort_nat({inParamsList.name});
    exp = 'searchParams_Nyquist\d.mat'; % regular expressions for desire file names
    inParamNames = regexp(inParamNames,exp,'match');
    inParamNames = inParamNames(~cellfun(@isempty,inParamNames)); % get rid of empty cells
    N = length(inParamNames);% number of bands
    %     N = 1; % For plotting only 1 band
    inDataNames = sort_nat({inDataList.name});
    
    %% data pre-processing
    bx = zeros(N,3);
    by = zeros(N,3);
    binsize = zeros(N,1);
    %     esize =zeros(N,1);
    etyband = 0; % band doesn't have sources inside.
    
    for i = 1:N
        % load bands and estimated data
        %         load([simParamsDir,filesep,inParamNames{i}{1}]); % for initial bands
        load([simParamsDir,filesep,simFileName,filesep,inParamNames{i}{1}]); % for optimal band
        ybin_up = [ybin_up searchParams.angular_velocity(1)];% saving the band boundary
        ybin_low = [ybin_low searchParams.angular_velocity(2)];% Find the sources with frequencies in specific band
        Indx = find(omega >= searchParams.angular_velocity(2) & ...
            omega <= searchParams.angular_velocity(1));
        
        
        if isempty(Indx)
            disp(["There's no signal injected in band ",num2str(i)]);
            etyband = i;
            continue
        end
        % Get their frequencies in Hz
        binsrcOmega = omega(Indx);
        y = [y binsrcOmega]; % total y for all bands.
        binsize(i) = length(binsrcOmega);
        by(i,1:binsize(i)) = binsrcOmega;
        %binsrcF = (binsrcOmega/(2*pi))/(24*3600*365);
        % Get their SNR
        binsrcSNR = snr_chr(Indx);
        x = [x binsrcSNR]; % total x for all bands
        bx(i,1:binsize(i)) = binsrcSNR;
        % Get the sky location
        sra = [sra alpha(Indx)]; % simulated RA
        sdec = [sdec delta(Indx)]; % simulated Dec
        for j = 1:num_ite
            load([estDataDir,filesep,char(inDataNames(j + num_ite * (i-1)))],'bestRealLoc');
            disp(['File: ',char(inDataNames(j + num_ite * (i-1))),' loaded']);
            path_to_estimatedData = [estDataDir,filesep,char(inDataNames(j + num_ite * (i-1)))];
            
            %% Estimated source
            %path_to_simulationData = '~/Research/PulsarTiming/SimDATA/Mauritius/GWBsimDataSKA/GWBsimDataSKASrlz1Nrlz9.mat';
            %     path_to_pulsar_catalog = 'survey_ska.mat';
            estFreq = bestRealLoc(3)/(2*pi*365*24*3600);% convert unit from yr-1 to Hz
            [sourceParams]=ColSrcParams(path_to_estimatedData, simParams.Np);
            [estSNR,estTimRes] = Amp2Snr(sourceParams,simParams,yr);
            sy = [sy estFreq]; % estimated frequency
            sx = [sx estSNR]; % estimated SNR
            %             esize(i) = length(sx);
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
    
    scalx = sx .* simParams.sd(1)/(100*10^(-9)); % rescale SNR into 100 ns.
    
    y = y/(2*pi*365*24*3600);
    by = by/(2*pi*365*24*3600); % change from rad/yr to Hz
    uplim = max(max(x),max(scalx))+30;
    binSNR = 0:1:uplim;
    binSNR_log = logspace(-5,3,length(binSNR));
    ybin_up = ybin_up/(2*pi*365*24*3600);
    ybin_up = repmat(ybin_up,length(binSNR),1);% stack itself vertically to broadcast to the dimension of x
    ybin_low = ybin_low/(2*pi*365*24*3600);
    ybin_low = repmat(ybin_low,length(binSNR),1);
    stage = 1:1:num_ite;
    sortsy = sort(sy); % sort sy in ascending order.
    %% Noise processing
    %     noisedir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/noise/results';
    %     disp('Processing noise-only data');
    %     [avgnoise] = noiseprocess(noisedir,simNoiseDir,num_ite,N);
    %simNoiseDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/noise';
    %noisefile = [simNoiseDir,filesep,'noise.mat'];
    %load(noisefile);
    
    % Save Estimated Sources Info
    save([estDataDir,filesep,simFileName,filesep,'EstSrc'],'sx','sy','ra','dec');
    
    %% plot settings
    %     prefix = [estDataDir,filesep,'fig',filesep,simDataFileNames{lp}];
    prefix = [estDataDir,filesep,'fig',filesep,simFileName]; % for a single realization
    mkdir(prefix);
    figname = [simFileName,'-iMBLT'];
    
    %% Plot
    figure(1)
    % yyaxis right
    % loglog(x,y,'o',sx,sy,'kd','MarkerSize',10);
    %     b1srcid = find(sy <= ybin_up(1,1)); % find sources in band 1.
    %     b2srcid = find(sy > ybin_up(2,1));
    %     plot(x,y,'ob',sx(b1srcid),sy(b1srcid),'sr',sx(b2srcid),sy(b2srcid),'dm');
    %     catalog sources according to band edges, compatible with "Union"
    %     results.
    % plot(x,y,'ob',sx(1:num_ite),sy(1:num_ite),'sr',sx(num_ite + 1:end),sy(num_ite + 1:end),'dm'); % catalog sources according to files loading sequence, compatible with number of sources in band 1 is the same as that in band 2.
    plot(x,y,'ob',sx,sy,'rs') % do not separete estimated sources into different band
    %     semilogx(x,y,'o',sx,sy,'s');
    % disp(sy)
    
    hold on
    
    % plot grid
    for k=1:N
        %                 semilogx(binSNR_log,ybin_up(:,k),'b-');
        plot(binSNR,ybin_up(:,k),'b-');
        %                 semilogx(binSNR_log,ybin_low(:,k),'b--');
        plot(binSNR,ybin_low(:,k),'b--');
    end
    hold off
    
    xlabel('SNR');
    xlim([0 uplim]);
    ylabel('Frequency');
    % legend('True','Estimated B1','Estimated B2','Location','northeast');
    legend('True','Estimated','Location','best');
    title(figname);
    saveas(gcf,[prefix,filesep,figname],'png');
    savefig([prefix,filesep,figname]);
    % save('estTimRes01.mat','estTimRes');
    
    for m = 1:N
        figure(10+m)
        % add label to each est. and true sources
        plot(bx(m,1:binsize(m)),by(m,1:binsize(m)),'ob',sx((m-1)*num_ite+1:m*num_ite),sy((m-1)*num_ite+1:m*num_ite),'sr')
        hold on
        plot(binSNR,ybin_up(:,m),'b-');
        plot(binSNR,ybin_low(:,m),'b--');
        hold off
        xlabel('SNR');
        xlim([0 uplim]);
        ylabel('Frequency');
        legend('True','Estimated','Location','northeast');
        title([figname,'-Band',num2str(m)]);
        [sorx,id_true] = sort(bx(m,1:binsize(m)),'descend');
        text(sorx+1,by(m,id_true),num2str((1:numel(bx(m,1:binsize(m))))'),'Color','#0072BD')
        text(sx((m-1)*num_ite+1:m*num_ite)+1,sy((m-1)*num_ite+1:m*num_ite),num2str((1:num_ite)'),'Color','#D95319')
        saveas(gcf,[prefix,filesep,figname,'Band',num2str(m)],'png');
        savefig([prefix,filesep,figname,'Band',num2str(m)]);
    end
    
    %  % SNR vs. Stage figure
    %     figure(2)
    %     Legend = {N,1};
    %     hold on
    %     for i = 1:N - etyband
    %         plot(stage,sx(num_ite*(i-1)+1:num_ite*i));
    % %         semilogy(stage,sx(num_ite*(i-1)+1:num_ite*i));
    %         Legend{i} = ['Band ', num2str(i)+etyband];
    %         if i == (N - etyband)
    %             plot(stage, avgnoise,'--k');
    % %             semilogy(stage,avgnoise,'--k');
    %             Legend{i+1} = 'Noise';
    %         end
    %     end
    %     hold off
    %     legend(Legend);
    %     title([figname,' SNR vs. Stage']);
    %     xlabel('Stage');
    %     xlim([1 num_ite]);
    %     ylabel('SNR');
    %     saveas(gcf,[prefix,filesep,figname,' SNR-Stage'],'png');
    %     savefig([prefix,filesep,figname,' SNR-Stage']);
    
    figure(3)
    plot(sra,sdec,'ob',ra,dec,'sr');
    title(['Sky location for ',figname]);
    legend('Simulated source','Estimated source');
    xlabel('RA');
    ylabel('Dec');
    saveas(gcf,[prefix,filesep,figname,' skyloc'],'png')
    savefig([prefix,filesep,figname,' skyloc']);
    
    figure(4)
    plot(sra,y,'ob',ra,sy,'sr');
    title(['RA vs. Freq for ',figname]);
    xlabel('RA');
    ylabel('Freq.');
    legend('True','Estimated')
    saveas(gcf,[prefix,filesep,figname,' RA'],'png')
    savefig([prefix,filesep,figname,' RA']);
    
    
    figure(5)
    plot(sdec,y,'ob',dec,sy,'sr');
    title(['DEC vs. Freq for ',figname]);
    xlabel('DEC');
    ylabel('Freq.');
    legend('True','Estimated')
    saveas(gcf,[prefix,filesep,figname,' DEC'],'png')
    savefig([prefix,filesep,figname,' DEC']);
    
    % Serial SNR plots
    figure
    plot(sx,'o')
    title('Serial SNR')
    xlabel('Serial Number')
    ylabel('SNR')
    saveas(gcf,[prefix,filesep,figname,'Serial_SNR'],'png')
    savefig([prefix,filesep,figname,'Serial_SNR'])
    
    % consecutive SNR diff
    figure
    plot(diff(sx),'o')
    hold on
    plot([0,length(sx)],[0,0],'r') % plot diff 0 line
    hold off
    title('Serial SNR Diff')
    xlabel('Serial Number')
    ylabel('SNR Diff')
    saveas(gcf,[prefix,filesep,figname,'Serial_SNR_Diff'],'png')
    savefig([prefix,filesep,figname,'Serial_SNR_Diff'])
    %% set up cutoff
    SNRcut = 50;
    Idx = find(x >= SNRcut);% set SNR cutoff for simulated sources
    Sx = x(Idx);
    Sy = y(Idx);
    sra = sra(Idx);
    sdec = sdec(Idx);
    
    Sidx = find(sx >= SNRcut);% SNR cutoff for est. sources
    sx = sx(Sidx);
    % sx = sx .* simParams.sd(1)/(100*10^(-9)); % rescale SNR into 100 ns.
    sy = sy(Sidx);
    ra = ra(Sidx);
    dec = dec(Sidx);
    
    figure(6)
    plot(sra,Sy,'ob',ra,sy,'sr');
    title(['RA vs. Freq for ',figname]);
    xlabel('RA');
    ylabel('Freq.');
    legend('True','Estimated')
    saveas(gcf,[prefix,filesep,figname,' RA-SNRcutoff ',num2str(SNRcut)],'png')
    savefig([prefix,filesep,figname,' RA-SNRcutoff ',num2str(SNRcut)]);
    
    
    figure(7)
    plot(sdec,Sy,'ob',dec,sy,'sr');
    title(['DEC vs. Freq for ',figname]);
    xlabel('DEC');
    ylabel('Freq.');
    legend('True','Estimated')
    saveas(gcf,[prefix,filesep,figname,' DEC-SNRcutoff ',num2str(SNRcut)],'png')
    savefig([prefix,filesep,figname,' DEC-SNRcutoff ',num2str(SNRcut)]);
    
    figure(8)
    plot(sra,sdec,'ob',ra,dec,'sr');
    title(['Sky location for ',figname]);
    legend('Simulated source','Estimated source');
    xlabel('RA');
    ylabel('Dec');
    saveas(gcf,[prefix,filesep,figname,' skyloc-SNRcutoff',num2str(SNRcut)],'png')
    savefig([prefix,filesep,figname,' skyloc']);
    
    figure(9)
    % yyaxis right
    % loglog(x,y,'o',sx,sy,'kd','MarkerSize',10);
    plot(Sx,Sy,'o',sx,sy,'s')
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
    saveas(gcf,[prefix,filesep,figname,' SNRcutoff ',num2str(SNRcut)],'png');
    savefig([prefix,filesep,figname,' SNRcutoff ',num2str(SNRcut)]);
    
    %% freq-only set up
    yb = Sy(2*ybin_up(1,1)/3 < Sy & Sy <= ybin_up(1,1));% simulated sources above SNR cutoff and blow the frequency boundary which separates the bands.
    yb2 = Sy(ybin_up(1,1)/3 < Sy & Sy <= 2*ybin_up(1,1));
    yb3 = Sy(Sy <= ybin_up(1,1)/3 );
    szyb = length(yb);
    szyb2 = length(yb2);
    szyb3 = length(yb3);
    cstb = zeros(szyb,1);
    cstb2 = zeros(szyb2,1);
    cstb3 = zeros(szyb3,1);
    yu = Sy(Sy > ybin_up(1,1));
    szyu = length(yu);
    cstu = zeros(szyu,1);
    
    syb = sy(2*ybin_up(1,1)/3 < sy & sy <= ybin_up(1,1)); % est. sources...
    syb2 = sy(ybin_up(1,1)/3 < sy & sy <= 2*ybin_up(1,1)/3);
    syb3 = sy(sy <= ybin_up(1,1)/3);
    szsyb = length(syb);
    szsyb2 = length(syb2);
    szsyb3 = length(syb3);
    scstb = zeros(szsyb,1);
    scstb2 = zeros(szsyb2,1);
    scstb3 = zeros(szsyb3,1);
    syu = sy(sy > ybin_up(1,1));
    szsyu = length(syu);
    scstu = zeros(szsyu,1);
    
    
    %% Freq-only plotting
    figure(10)
    subplot(2,2,3)
    plot(yb,cstb,'ob',syb,scstb,'sr');
    xlabel('Frequency');
    ylabel('Constant');
    title(' Band 1 upper');
    legend('Ture','Estimated');
    
    subplot(2,2,4)
    plot(yu,cstu,'ob',syu,scstu,'sr')
    xlabel('Frequency');
    ylabel('Constant');
    title(' Band 2');
    legend('Ture','Estimated');
    
    subplot(2,2,2)
    plot(yb2,cstb2,'ob',syb2,scstb2,'sr')
    xlabel('Frequency');
    ylabel('Constant');
    title(' Band 1 Medium');
    legend('Ture','Estimated');
    
    subplot(2,2,1)
    plot(yb3,cstb3,'ob',syb3,scstb3,'sr')
    xlabel('Frequency');
    ylabel('Constant');
    title(' Band 1 lower');
    legend('Ture','Estimated');
    
    sgtitle([figname,' Freq-only SNRcutoff ',num2str(SNRcut)]);
    saveas(gcf,[prefix,filesep,figname,' Freq-only SNRcutoff ',num2str(SNRcut)],'png')
    savefig([prefix,filesep,figname,' Freq-only SNRcutoff ',num2str(SNRcut)]);
    
    %% Frequency histogram
    dy = diff(sortsy);
    figure
    histogram(dy)
    saveas(gcf,[prefix,filesep,figname,'Diff-freq-histogram'],'png')
    savefig([prefix,filesep,figname,'Diff-freq-histogram']);
    figure
    plot(dy,sortsy(2:end),'ob')
    xlabel('Frequency differences')
    ylabel('Frequency')
    title('Dff. Frequency vs. Frequency')
    saveas(gcf,[prefix,filesep,figname,'Diff-Frequency'],'png')
    savefig([prefix,filesep,figname,'Diff-Frequency']);
    
    % close all figures
    close all
end

toc

%END
