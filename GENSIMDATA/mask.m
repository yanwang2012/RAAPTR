% "Mask" automate script

%% add extra noise to timing residual
clear;
simDataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Mask/sd_300_whole/subtract';
simFileName = 'GWBsimDataSKASrlz1Nrlz3';

load([simDataDir,filesep,simFileName,'_sub','_whole','.mat'])
N = simParams.N;
Np = simParams.Np;
stage = 2; % stages of Mask
etr_sd = zeros(Np,1);
etr_noise = zeros(Np,N);
for i = 1:Np
    etr_sd(i) = 2.0 * 10^(-7); % standard deviation of extra noise
    etr_noise(i,:) = etr_sd(i) * randn(1,N);
    timingResiduals(i,:) = timingResiduals(i,:) + etr_noise(i,:);
end
simParams.sd = sqrt(simParams.sd.^2 + etr_sd.^2);

outputDir = ['~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Mask/sd_',num2str(etr_sd(1)*10^9),'_whole'];
mkdir(outputDir)
newFile = [simFileName,'_Mask',num2str(stage)];
% copyfile([simDataDir,filesep,simFileName,'_sub','.mat'],[outputDir,filesep,newFile,'.mat']);
copyfile([simDataDir,filesep,simFileName,'_sub','_whole','.mat'],[outputDir,filesep,newFile,'.mat']);

m = matfile([outputDir,filesep,newFile],'Writable',true);
m.simParams = simParams;
m.timingResiduals = timingResiduals;

%% subtract est. sources from precedent subtrac data
clear;

PreSimDataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Mask/sd_489.9/subtract';
simDataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Mask/sd_300_whole';
simFileName = 'GWBsimDataSKASrlz1Nrlz3';
estDataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Mask/sd_300_whole/results_whole';
stage = 2;
load([simDataDir,filesep,simFileName,'_Mask',num2str(stage),'.mat'],'simParams','yr');
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
output = [simDataDir,filesep,'subtract'];
mkdir(output);
newFile = [output,filesep,simFileName,'_sub','_whole','.mat'];
copyfile([PreSimDataDir,filesep,simFileName,'_sub','_whole','.mat'],newFile); % can't substute '_sub' with '*', the later will copy the whole folder.
% copyfile([PreSimDataDir,filesep,simFileName,'.mat'],newFile); % for stage 1 specific.
m = matfile(newFile,'Writable',true);

m.timingResiduals = m.timingResiduals - timRes;

%END