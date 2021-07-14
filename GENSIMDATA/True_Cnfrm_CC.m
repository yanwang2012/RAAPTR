% Cross-Correlation Coefficients Matrix
% True sources vs. Confirmed sources which are caried out by using two
% different sets of results.
% Update to multi version

% Author: QYQ
% 06/24/2021

clear;
tic

%% Dir settings
searchParamsDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/searchParams/Whole';
simdataDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/simData/Band_opt_diff';
cnfrmdataDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/results_diff_opt_iMBLT';
Filename = 'GWBsimDataSKASrlz*Nrlz1';

SNR_threshold = 20;
cnfrmFilename = ['Confirmed_Src_Est_SNR',num2str(SNR_threshold)];
ext = '.mat';

%% Files
simFile = dir([simdataDir,filesep,Filename,ext]);
Nrlzs = length(simFile);
simFileNames = sort_nat({simFile.name});

%% Main
for rlz = 1:Nrlzs
    [~,simFileName,~] = fileparts(simFileNames{rlz});
    paraFile = dir([searchParamsDir,filesep,simFileName,filesep,'searchParams*',ext]);
    simFile = [simdataDir,filesep,simFileName,ext];
    % idFolder = 'id-Union2-xMBLT-vs-Union2-xMBLT-iMBLT';
    cnfrmFile = [cnfrmdataDir,filesep,simFileName,filesep,cnfrmFilename,ext];
    paraFilename = sort_nat({paraFile.name});
    exp = 'searchParams\d.mat'; % regular expressions for desire file names
    paraFilename = regexp(paraFilename,exp,'match');
    paraFilename = paraFilename(~cellfun(@isempty,paraFilename)); % get rid of empty cells
    Nband = length(paraFilename);
    
    
    load(simFile);
    load(cnfrmFile);
    confirm_src = CnfrmSrc_SNR;
    
    %% Seperate sources into different bands
    % Ntsrc = length(alpha); % Number of true sources.
    SrcSNR = {};
    SrcAlpha = {};
    SrcAmp = {};
    SrcDelta = {};
    SrcIota = {};
    SrcOmega = {};
    SrcPhi0 = {};
    SrcThetaN = {};
    
    for i = 1:Nband
        load([searchParamsDir,filesep,simFileName,filesep,paraFilename{i}{1}]);
        Indx = find(omega >= searchParams.angular_velocity(2) & ...
            omega <= searchParams.angular_velocity(1));
        
        SrcSNR{i} = snr_chr(Indx);
        SrcAlpha{i} = alpha(Indx);
        SrcDelta{i} = delta(Indx);
        SrcAmp{i} = Amp(Indx);
        SrcIota{i} = iota(Indx);
        SrcOmega{i} = omega(Indx);
        SrcPhi0{i} = phi0(Indx);
        SrcThetaN{i} = thetaN(Indx);
        
    end
    
    %% Sort sources in different bands
    
    for j = 1:Nband
        [~,id] = sort(SrcSNR{j},'descend'); % sort true sources according to SNR value
        SrcSNR{j} = SrcSNR{j}(id);
        SrcAlpha{j} = SrcAlpha{j}(id);
        SrcDelta{j} = SrcDelta{j}(id);
        SrcAmp{j} = SrcAmp{j}(id);
        SrcIota{j} = SrcIota{j}(id);
        SrcOmega{j} = SrcOmega{j}(id);
        SrcPhi0{j} = SrcPhi0{j}(id);
        SrcThetaN{j} = SrcThetaN{j}(id);
    end
    simSrc = struct('SrcSNR',SrcSNR,'SrcAlpha',SrcAlpha,'SrcDelta',SrcDelta,'SrcAmp',SrcAmp,...
        'SrcIota',SrcIota,'SrcOmega',SrcOmega,'SrcPhi0',SrcPhi0,'SrcThetaN',SrcThetaN); % Simulated sources parameters
    save([cnfrmdataDir,filesep,simFileName,filesep,'simSrc_Est'],'simSrc');
    % create new sets which excludes sources below certain criteria
    snr_cut = 5; % excludes true sources below snr 3
    RSrcSNR = {};
    RSrcAlpha = {};
    RSrcAmp = {};
    RSrcDelta = {};
    RSrcIota = {};
    RSrcOmega = {};
    RSrcPhi0 = {};
    RSrcThetaN = {};
    
    % apply snr cut off
    for j = 1:Nband
        % Assign SNR cut to simulated sources
        id_snr = SrcSNR{j} >= snr_cut;
        RSrcSNR{j} = SrcSNR{j}(id_snr);
        RSrcAlpha{j} = SrcAlpha{j}(id_snr);
        RSrcDelta{j} = SrcDelta{j}(id_snr);
        RSrcAmp{j} = SrcAmp{j}(id_snr);
        RSrcIota{j} = SrcIota{j}(id_snr);
        RSrcOmega{j} = SrcOmega{j}(id_snr);
        RSrcPhi0{j} = SrcPhi0{j}(id_snr);
        RSrcThetaN{j} = SrcThetaN{j}(id_snr);
    end
    RsimSrc = struct('SrcSNR',RSrcSNR,'SrcAlpha',RSrcAlpha,'SrcDelta',RSrcDelta,'SrcAmp',RSrcAmp,...
        'SrcIota',RSrcIota,'SrcOmega',RSrcOmega,'SrcPhi0',RSrcPhi0,'SrcThetaN',RSrcThetaN); % Simulated sources parameters but exclude weaker sources
    save([cnfrmdataDir,filesep,simFileName,filesep,'RsimSrc_Est'],'RsimSrc');
   
    %% Get identified sources info
    % idsrcBand1 = sum(~cellfun('isempty',idsrc(1,:))); % number of sources in a band.
    % idsrcBand2 = sum(~cellfun('isempty',idsrc(2,:)));
    % idsrcBand = struct('Band1',idsrcBand1,'Band2',idsrcBand2);
    % idsrcBand = length(confirm_src);
    
    %% Cross-Corelation
    % Normalized MTC
    psr_t = 0.9; % NMTC value threshold per-pulsar
    [rho,rho_max,id_max,estSNR] = NMTC(Nband,NcnfrmsrcBand,RsimSrc,confirm_src,simParams,yr,psr_t);
    
    
    % Minimum distance Maximum CC.
    % [rho,rho_max,dif_freq_max,dif_ra_max,dif_dec_max,id_max,estSNR] =
    % MinDMaxC(Nband,NestsrcBand,SrcAlpha,SrcDelta,SrcOmega,SrcPhi0,SrcIota,SrcThetaN,SrcAmp,EstSrc,simParams,yr);
    %     save([cnfrmdataDir,filesep,simFileName,filesep,'NMTC_Est_SNR',num2str(SNR_threshold)],'rho','rho_max');
    save([cnfrmdataDir,filesep,simFileName,filesep,'NMTC_Est_SNR',num2str(SNR_threshold),'tSNR_',num2str(snr_cut),'_psrT_',num2str(psr_t),ext],'rho','rho_max');
    
    %% Eliminating spurious sources
    t = 0.70; % NMTC threshold used to identify sources.
    id_src = {}; % identified sources.
    r = {}; % rows
    c = {}; % columns
    id_src_alpha = [];
    id_src_dec = [];
    id_src_snr = [];
    id_src_freq = [];
    for b = 1:Nband
        [r{b},c{b},~] = find(rho_max{b} > t); % r is the row of rho, c is the column of rho.
        % in rho, rows correspond to EstSrc2, columns correspond to EstSrc1.
        % select the identified sources from est. sources.
        for rr = 1:length(r{b})
            id_src{b,rr} = confirm_src{r{b}(rr)};
            id_src_alpha = [id_src_alpha confirm_src{r{b}(rr)}.alpha];
            id_src_dec = [id_src_dec confirm_src{r{b}(rr)}.delta];
            [cnfrm_snr,~] = Amp2Snr(confirm_src{r{b}(rr)},simParams,yr);
            id_src_snr = [id_src_snr cnfrm_snr];
            id_src_freq = [id_src_freq confirm_src{b,r{b}(rr)}.omega/(2*pi*365*24*3600)]; % convert from rad/yr to Hz
        end
    end
    
    NidsrcBand = zeros(Nband,1);
    for idb = 1:Nband
        NidsrcBand(idb) = sum(~cellfun('isempty',id_src(idb,:))); % # of identified sources in each band
    end
%     save([cnfrmdataDir,filesep,simFileName,filesep,'Identified_Src_SNR',num2str(SNR_threshold)],'id_src','NidsrcBand',...
%         'id_src_alpha','id_src_dec','id_src_freq','id_src_snr');
    save([cnfrmdataDir,filesep,simFileName,filesep,'Identified_Src_SNR',num2str(SNR_threshold),'tSNR_',num2str(snr_cut),'_psrT_',num2str(psr_t),ext],'id_src','NidsrcBand',...
        'id_src_alpha','id_src_dec','id_src_freq','id_src_snr','snr_cut','psr_t');
    
    %% Save matched true sources
    % save sky locations
    matched_alpha = []; % right ascension
    matched_alpha_cnfrm = [];
    matched_dec = []; % declination
    matched_dec_cnfrm = [];
    matched_snr = [];
    matched_snr_cnfrm = [];
    matched_freq = [];
    matched_freq_cnfrm = [];
    id_max_idty = zeros(size(id_max));
    cnfrm_src_nm = []; % confirm sources which don't find a true match
    
    for band = 1:Nband
        % matching true sources to confirmed sources
        matched_alpha_cnfrm = [matched_alpha_cnfrm SrcAlpha{band}(id_max(id_max(:,band) ~= 0, band))]; % exclude 0 elements
        matched_dec_cnfrm = [matched_dec_cnfrm SrcDelta{band}(id_max(id_max(:,band) ~= 0, band))];
        matched_snr_cnfrm = [matched_snr_cnfrm SrcSNR{band}(id_max(id_max(:,band) ~= 0, band))];
        matched_freq_cnfrm = [matched_freq_cnfrm SrcOmega{band}(id_max(id_max(:,band) ~= 0, band)) / (2*pi*24*365*3600)];
        % matching true sources to identified source
        matched_alpha = [matched_alpha SrcAlpha{band}(id_max(r{band},band))];
        matched_dec = [matched_dec SrcDelta{band}(id_max(r{band},band))];
        matched_snr = [matched_snr SrcSNR{band}(id_max(r{band},band))];
        matched_freq = [matched_freq SrcOmega{band}(id_max(r{band},band))/(2*pi*365*24*3600)]; % convert rad/yr to Hz
        id_max_idty(r{band},band) = id_max(r{band},band);
        cnfrm_src_nm = [cnfrm_src_nm confirm_src{id_max(:,band) == 0}]; 
    end
    
    %     save([cnfrmdataDir,filesep,simFileName,filesep,'Matched_Sources_Est_SNR',num2str(SNR_threshold),ext],'id_max','matched_alpha','matched_dec','matched_snr',...
    %         'matched_freq','SrcAlpha','SrcSNR','SrcDelta','id_max_idty','matched_alpha_cnfrm','matched_dec_cnfrm','matched_snr_cnfrm','matched_freq_cnfrm');
    save([cnfrmdataDir,filesep,simFileName,filesep,'Matched_Sources_Est_SNR',num2str(SNR_threshold),'tSNR_',num2str(snr_cut),'_psrT_',num2str(psr_t),ext],'id_max','matched_alpha','matched_dec','matched_snr',...
        'matched_freq','SrcAlpha','SrcSNR','SrcDelta','id_max_idty','matched_alpha_cnfrm','matched_dec_cnfrm','matched_snr_cnfrm','matched_freq_cnfrm', 'cnfrm_src_nm');
    %% Plotting
    metric = 'NMTC';
    methods = 'Confimred vs True';
    prefix = [cnfrmdataDir,filesep,'fig',filesep,simFileName,filesep,metric,'-',methods,'_SNR',num2str(SNR_threshold),'tSNR_',num2str(snr_cut),'_psrT_',num2str(psr_t)];
    mkdir(prefix);
    
    figname1 = metric;
    for fig = 1:Nband
        figure
        imagesc(rho_max{fig});
        a = colorbar;
        xlabel('True Source')
        ylabel('Confirmed Source')
        ylabel(a,'Corss-Correlation Coefficients')
        title(['Band ',num2str(fig)])
        saveas(gcf,[prefix,filesep,figname1,'Band ',num2str(fig)],'png');
        savefig([prefix,filesep,figname1,'Band ',num2str(fig)]);
    end
    
    
    figname2 = [metric,'-SNR'];
    for fig2 = 1:Nband
        N = NcnfrmsrcBand(fig2);
        figure
        plot(estSNR(fig2,1:N),max(rho_max{fig2},[],2),'ob')
        xlabel('Estimated SNR')
        ylabel(metric)
        title(['Band ',num2str(fig2)])
        saveas(gcf,[prefix,filesep,figname2,'Band ',num2str(fig2)],'png');
        savefig([prefix,filesep,figname2,'Band ',num2str(fig2)]);
    end
    
    close all
end

toc
%END