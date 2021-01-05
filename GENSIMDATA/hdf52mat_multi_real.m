% A script convert .hdf5 file to .mat file
% Extended for multiple realizations
% Author: QYQ
% 1/5/2021

clear;
inputFileDir = '/Users/qyq/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands';
folderName = dir([inputFileDir,filesep,'*xMBLT']);
folderName = sort_nat({folderName.name});
Nreal = 50; % number of realizations

for r = 1:Nreal
    % outputFileDir = [inputFileDir,filesep,'MAT'];
    inFileList = dir([inputFileDir,filesep,folderName{r},'*.hdf5']); % input file directory
    % mkdir(outputFileDir)% output directory
    for lpc = 1:length(inFileList)
        
        inFile = [inputFileDir,filesep,inFileList(lpc).name];
        
        [~,inFileName,~]  = fileparts(inFile);
        
        outFile = [inputFileDir,filesep,inFileName,'.mat'];
        
        mpavoutfile2mat(inFile,outFile);
        
    end
end