% A script convert searchParams file from .mat to hdf5

inFileList = dir('/Users/qianyiqian/Research/DATA/searchParams_GWBsimDataSKA/*.mat');
for lpc = 1:length(inFileList)
    inFile = ['/Users/qianyiqian/Research/DATA/searchParams_GWBsimDataSKA',filesep,inFileList(lpc).name];
    [~,inFileName,~] = fileparts(inFile);
    outFileName = ['/Users/qianyiqian/Research/DATA/searchParams_GWBsimDataSKA',filesep,inFileName,'.hdf5'];
    inFileInfo = load(inFile);
    fldName = 'xmaxmin';
    %Create HDF5 file
    fid = H5F.create(outFileName);
    H5F.close(fid);
    nDim = length(size((inFileInfo.(fldName))'));
    baseSz = ones(1, nDim);
    h5create(outFileName,['/',fldName],...
        max([size((inFileInfo.(fldName))'); baseSz]));
    h5write(outFileName,['/',fldName],(inFileInfo.(fldName))');
end