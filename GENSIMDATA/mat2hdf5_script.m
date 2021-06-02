% convert .mat file to .hdf5 file
clear;
DataDir = '/work2/05884/qyqstc/stampede2/MultiPSO/Final/realizations/2bands/results_diff_opt_xMBLT2/1_iMBLT/results/1iMBLT_after';
outDir = [DataDir,filesep,'HDF5'];
inFileList = dir([DataDir,filesep,'*GWB*.mat']); % input file directory
mkdir(outDir); % output directory
for lpc = 1:length(inFileList) 

inFile = [DataDir,filesep,inFileList(lpc).name]; 

[~,inFileName,~]  = fileparts(inFile); 

outFile = [outDir,filesep,inFileName,'.hdf5']; 

mpavinfile2hdf5(inFile,outFile); 
end
