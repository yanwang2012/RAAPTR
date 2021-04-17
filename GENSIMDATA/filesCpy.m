%% Copy input files for different band use
clear;
simDataDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/simData/diff_srlz_cos';
outputDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/simData/diff_srlz_cos';
% outputDir = simDataDir;
mkdir(outputDir);
inFiles = dir([simDataDir,filesep,'GWBsimDataSKASrlz*Nrlz1.mat']);
inFileNames = sort_nat({inFiles.name});
NumBand = 2;
N = length(inFiles);
for i = 1:N
    for j = 1:NumBand
        copyfile([simDataDir,filesep,inFileNames{i}],[outputDir,filesep,num2str(j) '_' inFileNames{i}])
    end
end