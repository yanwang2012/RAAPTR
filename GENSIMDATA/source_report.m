% Report estimated sources using different scheme

% Author: QYQ
% 04/13/2021

clear;
tic

%% Dir settings
searchParamsDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/searchParams/Band_opt/New';
simdataDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/Band_opt/simData';
estdataDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/Band_opt/iMBLT';
Filename = 'GWBsimDataSKASrlz1Nrlz3';
ext = '.mat';

%% Files
paraFile = dir([searchParamsDir,filesep,'searchParams','*.mat']);
simFile = [simdataDir,filesep,Filename,ext];
estFile = dir([estdataDir,filesep,'*',Filename,'*',ext]);
Nestsrc = length(estFile);

paraFilename = sort_nat({paraFile.name});
exp = 'searchParams_Nyquist\d.mat'; % regular expressions for desire file names
paraFilename = regexp(paraFilename,exp,'match');
paraFilename = paraFilename(~cellfun(@isempty,paraFilename)); % get rid of empty cells
Nband = length(paraFilename);

estFilename = sort_nat({estFile.name});
load(simFile);

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
    load([searchParamsDir,filesep,char(paraFilename{i})]);
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
save([estdataDir,filesep,'simSrc'],'simSrc');

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

save([estdataDir,filesep,'estSrc'],'EstSrc','EstSNR');

% report sources using SNR scheme
snr_trs = 15; % set SNR threshold
logits = EstSNR > snr_trs;
EstSrc_SNR = cell(size(EstSrc));
EstSrc_SNR(logits) = EstSrc(logits);
NrepsrcBand = zeros(Nband,1);
for idb = 1:Nband
    NrepsrcBand(idb) = sum(~cellfun('isempty',EstSrc_SNR(idb,:))); % # of identified sources in each band
end
save([estdataDir,filesep,'EstSrc_SNR'],'EstSrc_SNR','NrepsrcBand');

toc
         