% convert .mat file to .hdf5 file
clear;
DataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test10/FullBand/test/MBLT2';
outDir = [DataDir,filesep,'HDF5'];
inFileList = dir([DataDir,filesep,'*band*.mat']); % input file directory
mkdir(outDir); % output directory
for lpc = 1:length(inFileList) 

inFile = [DataDir,filesep,inFileList(lpc).name]; 

[~,inFileName,~]  = fileparts(inFile); 

outFile = [outDir,filesep,inFileName,'.hdf5']; 

mpavinfile2hdf5(inFile,outFile); 

end