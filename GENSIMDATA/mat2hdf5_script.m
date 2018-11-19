% convert .mat file to .hdf5 file

inFileList = dir('/Users/qianyiqian/Research/DATA/searchParams_GWBsimDataSKA/*.mat'); % input file directory
mkdir('/Users/qianyiqian/Research/DATA/searchParams_GWBsimDataSKA')% output directory
for lpc = 1:length(inFileList) 

inFile = ['/Users/qianyiqian/Research/DATA/searchParams_GWBsimDataSKA',filesep,inFileList(lpc).name]; 

[~,inFileName,~]  = fileparts(inFile); 

outFile = ['/Users/qianyiqian/Research/DATA/searchParams_GWBsimDataSKA_HDF5',filesep,inFileName,'.hdf5']; 

mpavinfile2hdf5(inFile,outFile); 

end 