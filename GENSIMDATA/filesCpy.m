%% Copy input files for different band use
clear;
simDataDir = '/Users/qyq/Research/PulsarTiming/YuYang_data';
% outputDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands/supperNarrow_iMBLT1_after_20/Results_20/2_iMBLT/results/2iMBLT_after/HDF5';
outputDir = simDataDir;
mkdir(outputDir);
inFiles = dir([simDataDir,filesep,'*.mat']);
NumBand = 5;
N = length(inFiles);
for i = 1:N
    for j = 1:NumBand
        copyfile([simDataDir,filesep,inFiles(i).name],[outputDir,filesep,num2str(j) '_' inFiles(i).name])
    end
end