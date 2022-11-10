% A script convert .hdf5 file to .mat file
clear;
inputFileDir = '/Users/yiqianqian/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/RAAPTR/MxAvPhaseC/TESTH5IO/results';
% outputFileDir = [inputFileDir,filesep,'MAT'];
inFileList = dir([inputFileDir,filesep,'NANOGrav12.5yv4_results*.hdf5']); % input file directory
% mkdir(outputFileDir)% output directory
for lpc = 1:length(inFileList) 

inFile = [inputFileDir,filesep,inFileList(lpc).name]; 

[~,inFileName,~]  = fileparts(inFile); 

outFile = [inputFileDir,filesep,inFileName,'.mat']; 

mpavoutfile2mat(inFile,outFile); 

end 
