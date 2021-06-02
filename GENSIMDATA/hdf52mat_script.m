% A script convert .hdf5 file to .mat file
clear;
inputFileDir = '/work2/05884/qyqstc/stampede2/MultiPSO/Final/realizations/2bands/results_diff_opt_xMBLT2';
% outputFileDir = [inputFileDir,filesep,'MAT'];
inFileList = dir([inputFileDir,filesep,'*.hdf5']); % input file directory
% mkdir(outputFileDir)% output directory
for lpc = 1:length(inFileList) 

inFile = [inputFileDir,filesep,inFileList(lpc).name]; 

[~,inFileName,~]  = fileparts(inFile); 

outFile = [inputFileDir,filesep,inFileName,'.mat']; 

mpavoutfile2mat(inFile,outFile); 

end 
