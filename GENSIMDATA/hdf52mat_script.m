% A script convert .hdf5 file to .mat file
inFileList = dir('~/Research/PulsarTiming/SimDATA/Mauritius/Mauritius_uni_results/*.hdf5'); % input file directory
mkdir('~/Research/PulsarTiming/SimDATA/Mauritius/Mauritius_uni_results_mat')% output directory
for lpc = 1:length(inFileList) 

inFile = ['~/Research/PulsarTiming/SimDATA/Mauritius/Mauritius_uni_results',...
          filesep,inFileList(lpc).name]; 

[~,inFileName,~]  = fileparts(inFile); 

outFile = ['~/Research/PulsarTiming/SimDATA/Mauritius/Mauritius_uni_results_mat',...
          filesep,inFileName,'.mat']; 

mpavoutfile2mat(inFile,outFile); 

end 