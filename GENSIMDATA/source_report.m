% Report estimated sources using different scheme
% Updated to handle multiple realizations

% Author: QYQ
% 05/27/2021

clear;
tic

%% Dir settings
searchParamsDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/searchParams/Whole';
simdataDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/simData/Band_opt_diff';
estdataDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/results_diff_one';
Filename = 'GWBsimDataSKASrlz*Nrlz1';
ext = '.mat';

%% Files
simFile = dir([simdataDir,filesep,Filename,ext]);
Nrlzs = length(simFile);
simFileNames = sort_nat({simFile.name});

%% Main
for rlz = 1:Nrlzs
    [~,simFileName,~] = fileparts(simFileNames{rlz}); 
    paraFile = dir([searchParamsDir,filesep,simFileName,filesep,'searchParams','*.mat']);
    estFile = dir([estdataDir,filesep,'*',simFileName,'*',ext]);
    Nestsrc = length(estFile);
    
    paraFilename = sort_nat({paraFile.name});
    exp = 'searchParams\d.mat'; % regular expressions for desire file names
    paraFilename = regexp(paraFilename,exp,'match');
    paraFilename = paraFilename(~cellfun(@isempty,paraFilename)); % get rid of empty cells
    Nband = length(paraFilename);
    
    estFilename = sort_nat({estFile.name});
    load([simdataDir,filesep,simFileName,ext]);
    
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
    folderName = [estdataDir,filesep,simFileName];
    if ~exist(folderName,'dir')
        mkdir(folderName)
    end
    save([folderName,filesep,'simSrc'],'simSrc');
    
    %% Get estimated sources info
    NestsrcBand = zeros(Nband,1);
    for nb = 1:Nband
        NestsrcBand(nb) = Nestsrc/Nband;
    end
    % NestsrcBand = struct('Band1',estsrcBand1,'Band2',estsrcBand2,'N',EN);
    EstSrc = {};
    EstSNR = zeros(Nband,Nestsrc/Nband);
    for band = 1:Nband
        estsrcBand = NestsrcBand(band);
        for k = 1:estsrcBand
            path_to_estimatedData = [estdataDir,filesep,char(estFilename((band - 1) * NestsrcBand(band) + k))];
            EstSrc{band,k} = ColSrcParams(path_to_estimatedData,simParams.Np);
            EstSNR(band,k) = Amp2Snr(EstSrc{band,k},simParams,yr);
        end
    end
    
    save([folderName,filesep,'estSrc'],'EstSrc','EstSNR');
    
    % report sources using SNR scheme
    snr_trs = 20; % set SNR threshold
    logits = EstSNR > snr_trs;
    RepSrc_SNR_tmp = cell(size(EstSrc)); % reported sources
    RepSrc_SNR = cell(size(EstSrc));
    RepSrc_SNR_tmp(logits) = EstSrc(logits);
    idx = ~cellfun('isempty',RepSrc_SNR_tmp);
    NrepsrcBand = sum(idx,2);
    % remove the empty cells in each row
    for band = 1:Nband
        RepSrc_SNR(band,1:NrepsrcBand(band)) = RepSrc_SNR_tmp(band,idx(band,:));
    end
    save([folderName,filesep,'RepSrc_SNR',num2str(snr_trs)],'RepSrc_SNR','NrepsrcBand');
end

toc
% END