% script to plot confirmed sources vs. True sources
% Confirmed sources are carried out by cross correlating reported sources
% with true sources.
% (1) 2D skymap
% (2) 3D skymap

% Author: QYQ
% 1/7/2021
clear;
%% Load data
simDataDir = '/Users/yiqianqian/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/simData/Band_opt_diff';
idDataDir = '/Users/yiqianqian/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/results_diff_opt_iMBLT';
Filename = 'GWBsimDataSKASrlz*Nrlz1';
%% IMPORTANT:MAKE SURE THIS IS CORRECT
SNR_Threshold = 7;
NMTC_t = 0.65;
tSNR = 5; % True sources SNR threshold
%%
confirmFilename = ['Confirmed_Src_SNR',num2str(SNR_Threshold),'tSNR_',num2str(tSNR),'_NMTC',num2str(NMTC_t)];
repFileName = ['RepSrc_SNR',num2str(SNR_Threshold)];
matchedFileName = ['Matched_Sources_SNR',num2str(SNR_Threshold),'tSNR_',num2str(tSNR),'_NMTC',num2str(NMTC_t)];
ext = '.mat';

%% Files
simFile = dir([simDataDir,filesep,Filename,ext]);
Nrlzs = length(simFile);
simFileNames = sort_nat({simFile.name});

for rlz = 1:Nrlzs
    [~,simFileName,~] = fileparts(simFileNames{rlz});
    simFile = [simDataDir,filesep,simFileName,ext];
    cfFile = [idDataDir,filesep,simFileName,filesep,confirmFilename,ext];
    reportFile = [idDataDir,filesep,simFileName,filesep,repFileName,ext];
    cnfrm2true = [idDataDir,filesep,simFileName,filesep,matchedFileName,ext];
    load(simFile);
    load(cfFile);
    load(cnfrm2true);
    load(reportFile);
    
    %% Sky location
    simRA_nm = []; % not matched true sources
    simDec_nm = [];
    
    cnfrmRA = [];
    repRA = [];
    cnfrmDec = [];
    repDec = [];
    cnfrmSNR = [];
    repSNR = [];
    cnfrmFreq = [];
    repFreq = [];
    Nband = length(NcnfrmsrcBand);
    % idBandSrc = zeros(Nband,1);
    
    % get # of confirmed sources in each band
    % for band = 1:Nband
    % idBandSrc(band) = sum(~cellfun('isempty',confirm_src(band,:)));
    % end
    
    for b = 1:Nband
        N = NcnfrmsrcBand(b);
        %     idx = setdiff(1:length(SrcAlpha{b}),id_max_cnfrm(:,b)); % get the
        %     index of not matched true sources corresponds to confirmed sources
        idx = setdiff(1:length(SrcAlpha{b}),id_max(:,b)); % .......... corresponds to reported sources.
        for i = 1:N
            cnfrmRA = [cnfrmRA confirm_src{b,i}.alpha];
            cnfrmDec = [cnfrmDec confirm_src{b,i}.delta];
            [idSNR_tmp,~] = Amp2Snr(confirm_src{b,i},simParams,yr);
            cnfrmSNR = [cnfrmSNR idSNR_tmp];
            cnfrmFreq = [cnfrmFreq confirm_src{b,i}.omega/(2*pi*24*365*3600)];
        end
        
        for j = 1: NrepsrcBand(b)
            repRA = [repRA RepSrc_SNR{b,j}.alpha];
            repDec = [repDec RepSrc_SNR{b,j}.delta];
            [repSNR_tmp,~] = Amp2Snr(RepSrc_SNR{b,j},simParams,yr);
            repSNR = [repSNR repSNR_tmp];
            repFreq = [repFreq RepSrc_SNR{b,j}.omega/(24*365*3600*2*pi)];
        end
        simRA_nm = [simRA_nm SrcAlpha{b}(idx)];
        simDec_nm = [simDec_nm SrcDelta{b}(idx)];
    end
    % Use all true sources
    %     save([idDataDir,filesep,simFileName,filesep,'simSrc_nm_sky',num2str(SNR_Threshold)],'simRA_nm','simDec_nm');
    %     save([idDataDir,filesep,simFileName,filesep,'matSrc_sky',num2str(SNR_Threshold)],'matched_alpha_rep','matched_dec_rep','matched_snr_rep',...,
    %         'matched_freq_rep','matched_alpha','matched_dec','matched_snr','matched_freq');
    %     save([idDataDir,filesep,simFileName,filesep,'cnfrmSrc_sky',num2str(SNR_Threshold)],'cnfrmRA','cnfrmDec','cnfrmSNR','cnfrmFreq');
    %     save([idDataDir,filesep,simFileName,filesep,'repSrc_sky',num2str(SNR_Threshold)],'repRA','repDec','repSNR','repFreq');
    
    % Use filtered true sources
    save([idDataDir,filesep,simFileName,filesep,'simSrc_nm_sky',num2str(SNR_Threshold),'tSNR_',num2str(tSNR),'_NMTC',num2str(NMTC_t),'.mat'],'simRA_nm','simDec_nm');
    save([idDataDir,filesep,simFileName,filesep,'matSrc_sky',num2str(SNR_Threshold),'tSNR_',num2str(tSNR),'_NMTC',num2str(NMTC_t),'.mat'],'matched_alpha_rep','matched_dec_rep','matched_snr_rep',...,
        'matched_freq_rep','matched_alpha','matched_dec','matched_snr','matched_freq');
    save([idDataDir,filesep,simFileName,filesep,'cnfrmSrc_sky',num2str(SNR_Threshold),'tSNR_',num2str(tSNR),'_NMTC',num2str(NMTC_t),'.mat'],'cnfrmRA','cnfrmDec','cnfrmSNR','cnfrmFreq');
    save([idDataDir,filesep,simFileName,filesep,'repSrc_sky',num2str(SNR_Threshold),'tSNR_',num2str(tSNR),'_NMTC',num2str(NMTC_t),'.mat'],'repRA','repDec','repSNR','repFreq');
    
    %% plot
    load('/Users/yiqianqian/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/GENSIMDATA/Acond for SKA/CondMap.mat'); % load skymap condition number
    
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
    scatter(ax2,repRA,repDec,[],'^')
    scatter(ax2,cnfrmRA,cnfrmDec,[],cnfrmSNR, 'filled') % with different color
    
    % pause
    
    % plot(matched_alpha,matched_dec,'ob') % plot truly matched true sources
    % scatter(matched_alpha,matched_dec,matched_snr,'ob')
    scatter(ax2,matched_alpha,matched_dec,[],matched_snr,'s')
    
    % pause
    % connect confirmed sources with true sources
    % for j = 1:length(idRA)
    %     plot(ax2,[idRA(j),matched_alpha(j)],[idDec(j),matched_dec(j)],'Color','m') % connect identified and matched sources
    %     %     pause
    % end
    
    % connect reported sources with true sources
    for j = 1:length(cnfrmRA)
        plot(ax2,[cnfrmRA(j),matched_alpha(j)],[cnfrmDec(j),matched_dec(j)],'Color','m') % connect identified and matched sources
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
    legend(ax2,{'True Srcs', 'Reported Srcs', 'Confirmed Srcs','Matched True Srcs','Matched & Report.'},'Location','best')
    
    saveas(gcf,[idDataDir,filesep,'fig',filesep,simFileName,filesep,'SkyLocationC_SNR',num2str(SNR_Threshold),'_NMTC',num2str(NMTC_t),'.png'])
    savefig([idDataDir,filesep,'fig',filesep,simFileName,filesep,'SkyLocationC_SNR',num2str(SNR_Threshold),'_NMTC',num2str(NMTC_t),'.fig'])
    close all
end

%% 3D sphere plot
% [x,y,z] = sphere;
% [cfr_x,cfr_y,cfr_z] = sph2cart(cnfrmRA,cnfrmDec,1);
% [sim_x,sim_y,sim_z] = sph2cart(simRA_nm,simDec_nm,1);
% [rep_x,rep_y,rep_z] = sph2cart(repRA,repDec,1);
% [mat_x_rep,mat_y_rep,mat_z_rep] = sph2cart(matched_alpha_rep,matched_dec_rep,1);
% 
% figure
% surf(x,y,z,'LineStyle','--','FaceColor','w','EdgeColor','k')
% 
% hold on
% 
% scatter3(cfr_x,cfr_y,cfr_z,40,cnfrmSNR,'filled');
% scatter3(sim_x,sim_y,sim_z,'o','MarkerEdgeColor','#888E8F')
% scatter3(rep_x,rep_y,rep_z,'^')
% scatter3(mat_x_rep,mat_y_rep,mat_z_rep,[],matched_snr_rep,'s');
% 
% center = zeros(1,3); % center of sphere
% for src = 1:length(repRA)
%     [v] = GreatCircle(center,[rep_x(src),rep_y(src),rep_z(src)],[mat_x_rep(src),mat_y_rep(src),mat_z_rep(src)],1);
%     plot3(v(1,:),v(2,:),v(3,:),'m')
%     % pause
% end
% hold off
% axis equal
% 
% title('3D Skymap')
% legend({'Grid','Confirmed Source','True Source','Reported Source','Matched True','Report to Matched'},'location','best')
% saveas(gcf,[idDataDir,filesep,'fig',filesep,'Skymap3D.png'])
% savefig([idDataDir,filesep,'fig',filesep,'Skymap3D'])

%% 3D sphere plot with Condition number
% load([idDataDir,filesep,'Acond.mat']); % load skymap condition numberload([idDataDir,filesep,'Acond.mat']); % load skymap condition number
%
% [x,y,z] = sphere(300); % generate a 300 faces unit sphere
% % convert sources coord. to spherical
% [id_x,id_y,id_z] = sph2cart(idRA,idDec,1);
% [sim_x,sim_y,sim_z] = sph2cart(simRA_nm,simDec_nm,1);
% [matched_x,matched_y,matched_z] = sph2cart(matched_alpha,matched_dec,1);
%
% figure
%
% ax1 = axes;
% surf(ax1,x,y,z,Acond,'EdgeColor','Interp','FaceColor','Interp');
% axis equal
%
% ax2 = axes;
% scatter3(ax2,id_x,id_y,id_z,40,idSNR,'filled');
% hold on
% scatter3(ax2,sim_x,sim_y,sim_z,'o','MarkerEdgeColor','#888E8F')
% scatter3(ax2,matched_x,matched_y,matched_z,[],matched_snr,'s');
%
% % connect identified sources with matched true sources with great circle on
% % sphere
% center = zeros(1,3); % center of sphere
% for src = 1:length(idRA)
%     [v] = GreatCircle(center,[id_x(src),id_y(src),id_z(src)],[matched_x(src),matched_y(src),matched_z(src)],1);
%     plot3(ax2,v(1,:),v(2,:),v(3,:),'m')
%     % pause
% end
%
% axis equal
%
% % Link two axes together
% hLink = linkprop([ax1,ax2],{'XLim','YLim','ZLim','CameraUpVector','CameraPosition','CameraTarget'});
%
% % Hide the top axes
% ax2.Visible = 'off';
% ax2.XTick = [];
% ax2.YTick = [];
%
% % Give each one its colormap
% colormap(ax1)
% colormap(ax2,'spring')
%
% % get everthin lined up
% cb1 = colorbar(ax1,'Position',[0.1 0.1 0.05 0.815]); % four-elements vector to specify Position [left bottom width height]
% cb2 = colorbar(ax2,'Position',[0.81 0.1 0.05 0.815]);
% cb1.Label.String = 'Condition Number';
% cb2.Label.String = 'SNR';
% cb1.Label.FontSize = 14;
% cb2.Label.FontSize = 14;
%
% legend(ax2,{'Identified Sources','True Sources','Matched True Source','Matched & Identi.'},'Location','southeast')
% setappdata(gcf,'StoreTheLink',hLink); % store the link so that they can rotate and zoom synchronically
%
% saveas(gcf,[idDataDir,filesep,'fig',filesep,'SkyLocation3D.png'])
% savefig([idDataDir,filesep,'fig',filesep,'SkyLocation3D'])

%END