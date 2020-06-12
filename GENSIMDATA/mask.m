% "Mask" automate script

%% add extra noise to timing residual
clear;
simDataDir = '/Users/qianyiqian/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Mask/sd_500/subtract';
simFileName = 'GWBsimDataSKASrlz1Nrlz3';

load([simDataDir,filesep,simFileName,'_sub','.mat'])
N = simParams.N;
stage = 2; % stages of Mask
etr_sd = 4.0 * 10^(-7); % standard deviation of extra noise
simParams.sd = sqrt(simParams.sd.^2 + etr_sd^2);
etr_noise = etr_sd * rand(1,N);
timingResiduals = timingResiduals_tmp + etr_noise;

outputDir = ['/Users/qianyiqian/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Mask/sd_',num2str(etr_sd*10^9)];
mkdir(outputDir)
newFile = [simFileName,'_Mask',num2str(stage)];
copyfile([simDataDir,filesep,simFileName,'_sub','.mat'],[outputDir,filesep,newFile,'.mat']);
m = matfile([outputDir,filesep,newFile],'Writable',true);
m.simParams = simParams;
m.timingResiduals = timingResiduals;

%% subtract est. sources from ori. sim data
clear;
OriSimDataDir = '/Users/qianyiqian/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Mask/';
simDataDir = '/Users/qianyiqian/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Mask/sd_500/';
simFileName = 'GWBsimDataSKASrlz1Nrlz3';
estDataDir = '/Users/qianyiqian/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Mask/sd_500/results';
load([simDataDir,filesep,simFileName]);
estFiles = dir([estDataDir,filesep,'*',simFileName,'*','.mat']);

nFiles = length(estFiles);
snrhold = 50; % SNR threshold.
estFileNames = sort_nat({estFiles.name});
ResCell = {}; % store estimated timing residuals
SNRarray = []; % store estimated SNR

% collect est. sources info
for i = 1:nFiles
    [srcParams] = ColSrcParams([estDataDir,filesep,char(estFileNames(i))]);
    [SNR,timRes_tmp] = Amp2Snr(srcParams,simParams,yr);
    ResCell = [ResCell timRes_tmp];
    SNRarray = [SNRarray SNR];
end

index = find(SNRarray >= snrhold);
Ni = length(index);
timRes = zeros(simParams.Np,simParams.N);

for idx = 1:Ni
    x = index(idx);
    timRes = timRes + ResCell{x};
end

% subtract from orig. file
output = '/Users/qianyiqian/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Mask/sd_500/subtract';
mkdir(output);
newFile = [output,filesep,simFileName,'_sub','.mat'];
copyfile([OriSimDataDir,filesep,simFileName,'.mat'],newFile);
m = matfile(newFile,'Writable',true);

m.timingResiduals = m.timingResiduals - timRes;

%END