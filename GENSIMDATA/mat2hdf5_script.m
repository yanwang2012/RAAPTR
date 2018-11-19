% convert .mat file to .hdf5 file

inFileList = dir('/Users/qianyiqian/Research/DATA/TESTDATA/*.mat'); % input file directory
mkdir('/Users/qianyiqian/Research/DATA/TESTDATA_HDF5')% output directory
for lpc = 1:length(inFileList) 

inFile = ['/Users/qianyiqian/Research/DATA/TESTDATA',filesep,inFileList(lpc).name]; 

[~,inFileName,~]  = fileparts(inFile); 

outFile = ['/Users/qianyiqian/Research/DATA/TESTDATA_HDF5',filesep,inFileName,'.hdf5']; 

mpavinfile2hdf5(inFile,outFile); 

end 