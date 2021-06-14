% convert .mat file to .hdf5 file
clear;
DataDir = '/work2/05884/qyqstc/stampede2/MultiPSO/Yuyang/1000Psr/results_xMBLT2/1_iMBLT/results/1iMBLT_after/results/2_iMBLT/results/2iMBLT_after/results/3_iMBLT/results/3iMBLT_after';
outDir = [DataDir,filesep,'HDF5'];
inFileList = dir([DataDir,filesep,'*GWB*.mat']); % input file directory
mkdir(outDir); % output directory
for lpc = 1:length(inFileList) 

inFile = [DataDir,filesep,inFileList(lpc).name]; 

[~,inFileName,~]  = fileparts(inFile); 

outFile = [outDir,filesep,inFileName,'.hdf5']; 

mpavinfile2hdf5(inFile,outFile); 
end
