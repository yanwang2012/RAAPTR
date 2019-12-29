% A script convert .hdf5 file to .mat file
clear;
inputFileDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test10/FullBand/test/MBLT2/Results';
outputFileDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test10/FullBand/test/MBLT2/Results';
inFileList = dir([inputFileDir,filesep,'*.hdf5']); % input file directory
mkdir(outputFileDir)% output directory
for lpc = 1:length(inFileList) 

inFile = [inputFileDir,filesep,inFileList(lpc).name]; 

[~,inFileName,~]  = fileparts(inFile); 

outFile = [outputFileDir,filesep,inFileName,'.mat']; 

mpavoutfile2mat(inFile,outFile); 

end 