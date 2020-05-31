% convert .mat file to .hdf5 file
clear;
DataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands/supperNarrow_iMBLT1_after_20/Results_20/2_iMBLT/results/2iMBLT_after/results/3_iMBLT/results/3iMBLT_after/results/4_iMBLT/results/4iMBLT_after/results/5_iMBLT/results/5iMBLT_after/results/6_iMBLT/results/6iMBLT_after/results/7_iMBLT/results/7iMBLT_after/results/8_iMBLT';
outDir = [DataDir,filesep,'HDF5'];
inFileList = dir([DataDir,filesep,'*GWB*.mat']); % input file directory
mkdir(outDir); % output directory
for lpc = 1:length(inFileList) 

inFile = [DataDir,filesep,inFileList(lpc).name]; 

[~,inFileName,~]  = fileparts(inFile); 

outFile = [outDir,filesep,inFileName,'.hdf5']; 

mpavinfile2hdf5(inFile,outFile); 
end