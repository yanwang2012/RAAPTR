% Cross-Correlation Coefficients analysis for two sets of est. sources.
% without spliting into different bands.

% Author: QYQ
% 06/22/2021

clear;
tic

%% Dir settings
simdataDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/simData/Band_opt_diff';
estSrc1Dir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/results_diff_opt_xMBLT';
estsrc1 = 'xMBLT';
estSrc2Dir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/results_diff_opt_iMBLT';
estsrc2 = 'iMBLT';
Filename = 'GWBsimDataSKASrlz*Nrlz1';
ext = '.mat';

%% Files
estSrc1File = dir([estSrc1Dir,filesep,'*',Filename,'*',ext]);
estSrc2File = dir([estSrc2Dir,filesep,'*',Filename,'*',ext]);
simFile = dir([simdataDir,filesep,Filename,'*',ext]);
estSrc2Filename = sort_nat({estSrc2File.name});
estSrc1Filename = sort_nat({estSrc1File.name});
simFilename = sort_nat({simFile.name});

Nrlzs = length(simFile);

for rlz = 1:Nrlzs
    load([simdataDir,filesep,simFilename{rlz}]);
    %% Get estimated sources info
    EstSrc2 = {};
    EstSrc1 = {};
    [~,baseName] = fileparts(simFilename{rlz});
    exp1 = ['\d+_',baseName,'(?=_|\.mat)_?\d{0,2}\.mat']; % for initial band selection
    estSrc1File_tmp = regexp(estSrc1Filename,exp1,'match');
    estSrc1File_tmp = estSrc1File_tmp(~cellfun(@isempty, estSrc1File_tmp));
    exp2 = ['\d+_',baseName,'_\d{0,2}iMBLT\.mat']; % for iMBLT results
    estSrc2File_tmp = regexp(estSrc2Filename,exp2,'match');
    estSrc2File_tmp = estSrc2File_tmp(~cellfun(@isempty, estSrc2File_tmp));
    Nite1 = length(estSrc1File_tmp);
    Nite2 = length(estSrc2File_tmp);
    for k = 1:Nite1
        path_to_estimatedDataestSrc1 = [estSrc1Dir,filesep,estSrc1File_tmp{k}{1}];
        EstSrc1{k} = ColSrcParams(path_to_estimatedDataestSrc1, simParams.Np);
    end
    
    for j = 1:Nite2
        path_to_estimatedDataestSrc2 = [estSrc2Dir,filesep,estSrc2File_tmp{j}{1}];
        EstSrc2{j} = ColSrcParams(path_to_estimatedDataestSrc2, simParams.Np);
    end
    % Cross-Corelation
    [gamma,rho,id_max,estSNR1,estSNR2] = ESNMTCW(EstSrc1,EstSrc2,simParams,yr,0.90);
    
    %% Eliminating spurious sources
    t = 0.70; % NMTC threshold used to identify sources.
    confirm_src = {}; % confirmed sources
    cnfrm_src_snr = [];
    
    [r,c,~] = find(gamma > t); % r is the row of rho, c is the column of gamma.
    % in gamma, rows correspond to EstSrc1, columns correspond to EstSrc2.
    % select the identified sources from est. sources.
    for rr = 1:length(r)
        confirm_src{rr} = EstSrc2{c(rr)};
        cnfrm_src_snr = [cnfrm_src_snr estSNR2(c(rr))];
    end
    
    % select confirmed sources using SNR threshold
    snr_trs = 20; % set SNR threshold
    logits = cnfrm_src_snr > snr_trs;
    CnfrmSrc_SNR_tmp = cell(size(confirm_src)); % reported sources
    CnfrmSrc_SNR = cell(size(confirm_src));
    CnfrmSrc_SNR_tmp(logits) = confirm_src(logits);
    idx = ~cellfun('isempty',CnfrmSrc_SNR_tmp);
    NcnfrmsrcBand = sum(idx,2);
    % remove the empty cells in each row
    
    CnfrmSrc_SNR(1:NcnfrmsrcBand) = CnfrmSrc_SNR_tmp(idx);
    CnfrmSrc_SNR = CnfrmSrc_SNR(~cellfun('isempty',CnfrmSrc_SNR)); % remove trailing blank cells
    
    mkdir([estSrc2Dir,filesep,baseName]);
    methods = 'xMBLT'; % use 'Est' for one band + xMBLT, while 'xMBLT' for xMBLT + iMBLT
    save([estSrc2Dir,filesep,baseName,filesep,'Confirmed_Src_xMBLT_SNR',num2str(snr_trs)],'NcnfrmsrcBand','CnfrmSrc_SNR','snr_trs');
    
end

toc

%END