% A script convert .hdf5 file to .mat file
clear;
inputFileDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/results_diff_one';
% outputFileDir = [inputFileDir,filesep,'MAT'];
inFileList = dir([inputFileDir,filesep,'*.hdf5']); % input file directory
% mkdir(outputFileDir)% output directory
for lpc = 1:length(inFileList) 

inFile = [inputFileDir,filesep,inFileList(lpc).name]; 

[~,inFileName,~]  = fileparts(inFile); 

outFile = [inputFileDir,filesep,inFileName,'.mat']; 

mpavoutfile2mat(inFile,outFile); 

end 
