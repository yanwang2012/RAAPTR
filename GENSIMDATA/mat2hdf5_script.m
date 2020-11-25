% convert .mat file to .hdf5 file
clear;
DataDir = '/work/05884/qyqstc/lonestar/MultiPSO/Task8/BANDEDGE/2bands/superNarrow/Union2_xMBLT/results/Union2_xMBLT2/results/1_iMBLT/results';
outDir = [DataDir,filesep,'HDF5'];
inFileList = dir([DataDir,filesep,'*GWB*.mat']); % input file directory
mkdir(outDir); % output directory
for lpc = 1:length(inFileList) 

inFile = [DataDir,filesep,inFileList(lpc).name]; 

[~,inFileName,~]  = fileparts(inFile); 

outFile = [outDir,filesep,inFileName,'.hdf5']; 

mpavinfile2hdf5(inFile,outFile); 
end
