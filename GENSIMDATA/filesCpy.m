%% Copy input files for different band use
clear;
simDataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands/supperNarrow_iMBLT1_after_20/Results_20/2_iMBLT/results/2iMBLT_after/results/3_iMBLT/results/3iMBLT_after/results/4_iMBLT/results/4iMBLT_after/results/5_iMBLT/results/5iMBLT_after/results/6_iMBLT/results/6iMBLT_after/results/7_iMBLT/results/7iMBLT_after/HDF5';
% outputDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands/supperNarrow_iMBLT1_after_20/Results_20/2_iMBLT/results/2iMBLT_after/HDF5';
outputDir = simDataDir;
mkdir(outputDir);
inFiles = dir([simDataDir,filesep,'*.hdf5']);
NumBand = 2;
N = length(inFiles);
for i = 1:N
    for j = 1:NumBand
        copyfile([simDataDir,filesep,inFiles(i).name],[outputDir,filesep,num2str(j) '_' inFiles(i).name])
    end
end