% convert .mat file to .hdf5 file
% extended to handle multiple realizations
% Author: QYQ
% 1/5/2021

clear;
DataDir = '/Users/qyq/Research/PulsarTiming/YuYang_data/simData';
folderName = dir([DataDir,filesep,'*xMBLT']);
folderName = sort_nat({folderName.name});
Nreal = length(folderName); % # of realizations

for r = 1:Nreal
    outDir = [DataDir,filesep,folderName{r},filesep,'HDF5'];
    inFileList = dir([DataDir,filesep,folderName{r},filesep,'*GWB*.mat']); % input file directory
    mkdir(outDir); % output directory
    for lpc = 1:length(inFileList)
        
        inFile = [DataDir,filesep,folderName{r},filesep,inFileList(lpc).name];
        
        [~,inFileName,~]  = fileparts(inFile);
        
        outFile = [outDir,filesep,inFileName,'.hdf5'];
        
        mpavinfile2hdf5(inFile,outFile);
    end
end