% A script convert .hdf5 file to .mat file
clear;
inputFileDir = '/work/05884/qyqstc/lonestar/MultiPSO/Task8/BANDEDGE/2bands/SupNar_xMBLT_iMBLT20/GWBsimDataSKASrlz1Nrlz3_xMBLT2/results/1_iMBLT/results/1iMBLT_after/results/2_iMBLT/results/2iMBLT_after/results/3_iMBLT/results/3iMBLT_after/results/4_iMBLT/results/4iMBLT_after/results/5_iMBLT/results/5iMBLT_after/results/6_iMBLT/results/6iMBLT_after/results/7_iMBLT/results/7iMBLT_after/results/8_iMBLT/results/8iMBLT_after/results/9_iMBLT/results/9iMBLT_after/results/10_iMBLT/results/10iMBLT_after/results/11_iMBLT/results/11iMBLT_after/results/12_iMBLT/results/12iMBLT_after/results/13_iMBLT/results/13iMBLT_after/results/14_iMBLT/results/14iMBLT_after/results/15_iMBLT/results/15iMBLT_after/results/16_iMBLT/results/16iMBLT_after/results/17_iMBLT/results/17iMBLT_after/results/18_iMBLT/results/18iMBLT_after/results/19_iMBLT/results/19iMBLT_after/results';
% outputFileDir = [inputFileDir,filesep,'MAT'];
inFileList = dir([inputFileDir,filesep,'*.hdf5']); % input file directory
% mkdir(outputFileDir)% output directory
for lpc = 1:length(inFileList) 

inFile = [inputFileDir,filesep,inFileList(lpc).name]; 

[~,inFileName,~]  = fileparts(inFile); 

outFile = [inputFileDir,filesep,inFileName,'.mat']; 

mpavoutfile2mat(inFile,outFile); 

end 
