% convert .mat file to .hdf5 file

inFileList = dir('TESTDATA/*.mat'); 
mkdir('TESTDATA_HDF5')
for lpc = 1:length(inFileList) 

inFile = ['TESTDATA',filesep,inFileList(lpc).name]; 

[~,inFileName,~]  = fileparts(inFile); 

outFile = ['TESTDATA_HDF5',filesep,inFileName,'.hdf5']; 

mpavinfile2hdf5(inFile,outFile); 

end 