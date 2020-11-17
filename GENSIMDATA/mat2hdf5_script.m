% convert .mat file to .hdf5 file
clear;
DataDir = '/Users/qyq/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands/SuperNarrow/Results_supNar_rand1/GWBsimDataSKASrlz1Nrlz3_xMBLT/results/Union/GWBsimDataSKASrlz1Nrlz3_Union_xMBLT3';
outDir = [DataDir,filesep,'HDF5'];
inFileList = dir([DataDir,filesep,'*GWB*.mat']); % input file directory
mkdir(outDir); % output directory
for lpc = 1:length(inFileList) 

inFile = [DataDir,filesep,inFileList(lpc).name]; 

[~,inFileName,~]  = fileparts(inFile); 

outFile = [outDir,filesep,inFileName,'.hdf5']; 

mpavinfile2hdf5(inFile,outFile); 
end
