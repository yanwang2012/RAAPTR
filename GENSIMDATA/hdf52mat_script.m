% A script convert .hdf5 file to .mat file
clear;
inputFileDir = '/Users/qyq/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands/SuperNarrow/Results_supNar_rand1/GWBsimDataSKASrlz1Nrlz3_xMBLT/results/Union/GWBsimDataSKASrlz1Nrlz3_Union_xMBLT3/results';
% outputFileDir = [inputFileDir,filesep,'MAT'];
inFileList = dir([inputFileDir,filesep,'*.hdf5']); % input file directory
% mkdir(outputFileDir)% output directory
for lpc = 1:length(inFileList) 

inFile = [inputFileDir,filesep,inFileList(lpc).name]; 

[~,inFileName,~]  = fileparts(inFile); 

outFile = [inputFileDir,filesep,inFileName,'.mat']; 

mpavoutfile2mat(inFile,outFile); 

end 
