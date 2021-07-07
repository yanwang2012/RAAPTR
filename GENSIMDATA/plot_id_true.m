% script to plot identified sources vs. True sources
% Identified sources are carried out by cross correlating confirmed sources
% with true sources.
% (1) 2D skymap
% (2) 3D skymap

% Author: QYQ
% 1/7/2021
clear;
%% Load data
simDataDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/simData/Band_opt_diff';
DataDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/results_diff_opt_xMBLT';
Filename = 'GWBsimDataSKASrlz*Nrlz1';
%%
SNR_threshold = 20;
tSNR_cut = 5;
psr_t = 0.6;

confirmFilename = ['Confirmed_Src_Est_SNR',num2str(SNR_threshold)];
% Use all true sources
% idtyFilename = ['Identified_Src_SNR',num2str(SNR_threshold)];
% matchedFilename = ['Matched_Sources_Est_SNR',num2str(SNR_threshold)];
% Use filtered true sources
idtyFilename = ['Identified_Src_SNR',num2str(SNR_threshold),'tSNR_',num2str(tSNR_cut),'_psrT_',num2str(psr_t)];
matchedFilename = ['Matched_Sources_Est_SNR',num2str(SNR_threshold),'tSNR_',num2str(tSNR_cut),'_psrT_',num2str(psr_t)];
ext = '.mat';

%% Files
simFile = dir([simDataDir,filesep,Filename,ext]);
Nrlzs = length(simFile);
simFileNames = sort_nat({simFile.name});

for rlz = 1:Nrlzs
    [~,simFileName,~] = fileparts(simFileNames{rlz});
    simFile = [simDataDir,filesep,simFileName,ext];
    cfFile = [DataDir,filesep,simFileName,filesep,confirmFilename,ext];
    idFile = [DataDir,filesep,simFileName,filesep,idtyFilename,ext];
    id2true = [DataDir,filesep,simFileName,filesep,matchedFilename,ext];
    load(simFile);
    load(cfFile);
    load(id2true);
    load(idFile);
    
    %% Sky location
    simRA_nm = []; % not matched true sources
    simDec_nm = [];
    simSNR_nm = [];
    
    cnfrmRA = [];
    idRA = [];
    cnfrmDec = [];
    idDec = [];
    cnfrmSNR = [];
    idSNR = [];
    Nband = length(NidsrcBand);
    confirm_src = CnfrmSrc_SNR;
    cnfrmFreq = [];
    idFreq = [];
    % NcnfrmsrcBand = length(confirm_src);
    % idBandSrc = zeros(Nband,1);
    
    % get # of confirmed sources in each band
    % for band = 1:Nband
    % idBandSrc(band) = sum(~cellfun('isempty',confirm_src(band,:)));
    % end
    
    for b = 1:Nband
        %     idx = setdiff(1:length(SrcAlpha{b}),id_max_cnfrm(:,b)); % get the
        %     index of not matched true sources corresponds to confirmed sources
        idx = setdiff(1:length(SrcAlpha{b}),id_max_idty(:,b)); % .......... corresponds to identified sources.
        for i = 1:NcnfrmsrcBand(b)
            cnfrmRA = [cnfrmRA confirm_src{b,i}.alpha];
            cnfrmDec = [cnfrmDec confirm_src{b,i}.delta];
            [SNR_tmp,~] = Amp2Snr(confirm_src{b,i},simParams,yr);
            cnfrmSNR = [cnfrmSNR SNR_tmp];
            cnfrmFreq = [cnfrmFreq confirm_src{b,i}.omega/(2*pi*24*365*3600)];
        end
        
        for j = 1: NidsrcBand(b)
            idRA = [idRA id_src{b,j}.alpha];
            idDec = [idDec id_src{b,j}.delta];
            [SNR_tmp,~] = Amp2Snr(id_src{b,j},simParams,yr);
            idSNR = [idSNR SNR_tmp];
            idFreq = [idFreq id_src{b,j}.omega/(2*pi*24*365*3600)];
        end
        simRA_nm = [simRA_nm SrcAlpha{b}(idx)];
        simDec_nm = [simDec_nm SrcDelta{b}(idx)];
        simSNR_nm = [simSNR_nm SrcSNR{b}(idx)];
    end
    % Use all true sources
    %     save([DataDir,filesep,simFileName,filesep,'simSrc_nm_sky_Est_SNR',num2str(SNR_threshold)],'simRA_nm','simDec_nm','simSNR_nm');
    %     save([DataDir,filesep,simFileName,filesep,'matSrc_sky_Est_SNR',num2str(SNR_threshold)],'matched_alpha_cnfrm','matched_dec_cnfrm','matched_snr_cnfrm', ...,
    %         'matched_freq_cnfrm','matched_alpha','matched_dec','matched_snr','matched_freq');
    %     save([DataDir,filesep,simFileName,filesep,'cnfrmSrc_sky_Est_SNR',num2str(SNR_threshold)],'cnfrmRA','cnfrmDec','cnfrmSNR','cnfrmFreq');
    %     save([DataDir,filesep,simFileName,filesep,'idSrc_sky_Est_SNR',num2str(SNR_threshold)],'idRA','idDec','idSNR','idFreq');
    
    % Use filtered true sources
    save([DataDir,filesep,simFileName,filesep,'simSrc_nm_sky_Est_SNR',num2str(SNR_threshold),'tSNR_',num2str(tSNR_cut),'_psrT_',num2str(psr_t),ext],'simRA_nm','simDec_nm','simSNR_nm');
    save([DataDir,filesep,simFileName,filesep,'matSrc_sky_Est_SNR',num2str(SNR_threshold),'tSNR_',num2str(tSNR_cut),'_psrT_',num2str(psr_t),ext],'matched_alpha_cnfrm','matched_dec_cnfrm','matched_snr_cnfrm', ...,
        'matched_freq_cnfrm','matched_alpha','matched_dec','matched_snr','matched_freq');
    save([DataDir,filesep,simFileName,filesep,'cnfrmSrc_sky_Est_SNR',num2str(SNR_threshold),'tSNR_',num2str(tSNR_cut),'_psrT_',num2str(psr_t),ext],'cnfrmRA','cnfrmDec','cnfrmSNR','cnfrmFreq');
    save([DataDir,filesep,simFileName,filesep,'idSrc_sky_Est_SNR',num2str(SNR_threshold),'tSNR_',num2str(tSNR_cut),'_psrT_',num2str(psr_t),ext],'idRA','idDec','idSNR','idFreq');
    
    %% plot
    load('/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/GENSIMDATA/Acond for SKA/CondMap.mat'); % load skymap condition number
    
    figure
    % plot condition numbers
    ax1 = axes;
    imagesc(ax1,condMap.RA, condMap.Dec, condMap.cond);
    
    xlabel(ax1,'\alpha')
    ylabel(ax1,'\delta')
    title(ax1,'Sky Location')
    
    % plot(simRA,simDec,'o','MarkerEdgeColor','##D4D8D9') % plot true sources
    ax2 = axes;
    scatter(ax2,simRA_nm,simDec_nm,'o','MarkerEdgeColor','#D4D8D9') % use scatter to make marker size change
    
    hold on
    
    % plot(idRA,idDec,'rs') % plot identified sources
    % scatter(idRA,idDec,idSNR,'rs') % with different size
    scatter(ax2,idRA,idDec,[],'^')
    scatter(ax2,cnfrmRA,cnfrmDec,[],cnfrmSNR, 'filled') % with different color
    
    % pause
    
    % plot(matched_alpha,matched_dec,'ob') % plot truly matched true sources
    % scatter(matched_alpha,matched_dec,matched_snr,'ob')
    scatter(ax2,matched_alpha_cnfrm,matched_dec_cnfrm,[],matched_snr_cnfrm,'s')
    
    % pause
    % connect confirmed sources with true sources
    % for j = 1:length(idRA)
    %     plot(ax2,[idRA(j),matched_alpha(j)],[idDec(j),matched_dec(j)],'Color','m') % connect identified and matched sources
    %     %     pause
    % end
    
    % connect reported sources with true sources
    for j = 1:length(cnfrmRA)
        plot(ax2,[cnfrmRA(j),matched_alpha_cnfrm(j)],[cnfrmDec(j),matched_dec_cnfrm(j)],'Color','m') % connect identified and matched sources
        %     pause
    end
    
    xlim(ax2,[0 2*pi])
    
    %     Link two axes together
    linkaxes([ax1,ax2])
    
    %     Hide the top axis
    ax2.Visible = 'off';
    ax2.XTick = [];
    ax2.YTick = [];
    set([ax1,ax2],'Position',[.17 .11 .685 .815]);
    
    
    colormap(ax1)
    % cb1 = colorbar(ax1,'Location','northoutside');
    cb1 = colorbar(ax1,'Position',[.05 .11 .04 .815]);
    cb1.Label.String = 'Condition Number';
    cb1.Label.FontSize = 14;
    
    cmap = colormap(ax2,'autumn');
    colormap(ax2,flipud(cmap)); % flip 'autumn' upside down
    %     cb2 = colorbar(ax2,'Location','southoutside');
    cb2 = colorbar(ax2,'Position',[.90 .11 .04 .815]);
    cb2.Label.String = 'SNR';
    cb2.Label.FontSize = 14;
    %     ylabel(ax2,'\delta');
    %     xlabel(ax2,'\alpha');
    legend(ax2,{'True Srcs', 'Identified Srcs', 'Confirmed Srcs','Matched True Srcs','Matched & Confrm.'},'Location','best')
    
    saveas(gcf,[DataDir,filesep,'fig',filesep,simFileName,filesep,'SkyLocationC_Est_SNR',num2str(SNR_threshold),'tSNR_',num2str(tSNR_cut),'_psrT_',num2str(psr_t),'.png'])
    savefig([DataDir,filesep,'fig',filesep,simFileName,filesep,'SkyLocationC_Est_SNR',num2str(SNR_threshold),'tSNR_',num2str(tSNR_cut),'_psrT_',num2str(psr_t),'.fig'])
    close all
end

%END