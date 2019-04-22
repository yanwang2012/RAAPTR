% A script convert searchParams file from .mat to hdf5
clear;
inFileList = dir('~/Research/PulsarTiming/SimDATA/Mauritius/uniform search params/*.mat');
mkdir('~/Research/PulsarTiming/SimDATA/Mauritius/uniform search params HDF5');
for lpc = 1:length(inFileList)
    inFile = ['~/Research/PulsarTiming/SimDATA/Mauritius/uniform search params',filesep,inFileList(lpc).name];
    [~,inFileName,~] = fileparts(inFile);
    outFileName = ['~/Research/PulsarTiming/SimDATA/Mauritius/uniform search params HDF5',filesep,inFileName,'.hdf5'];
    inFileInfo = load(inFile);
    struName = 'searchParams';
    fldName = 'xmaxmin';
    bandsIn = 'Bands_Info';
    %Create HDF5 file
    fid = H5F.create(outFileName);
    H5F.close(fid);
    nDim = length(size((inFileInfo.(struName).alpha)'));
    baseSz = ones(1, nDim);
    
    h5create(outFileName,['/',bandsIn],2);
    h5write(outFileName,['/',bandsIn],inFileInfo.(struName).band_num,1,1);
    h5writeatt(outFileName,['/',bandsIn],'bandNumber',inFileInfo.(struName).band_num);
    h5write(outFileName,['/',bandsIn],inFileInfo.NumBands,2,1);
    h5writeatt(outFileName,['/',bandsIn],'NumBands',inFileInfo.NumBands);
    
%     h5create(outFileName,['/',fldName],...
%         max([size((inFileInfo.(struName).alpha)'); baseSz]));
    h5create(outFileName,['/',fldName],[7,2])
    h5write(outFileName,['/',fldName],(inFileInfo.(struName).alpha)',[1,1],[1,2]);
    h5writeatt(outFileName,['/',fldName],'alpha',(inFileInfo.(struName).alpha)');
    
    h5write(outFileName,['/',fldName],(inFileInfo.(struName).delta)',[2,1],[1,2]);
    h5writeatt(outFileName,['/',fldName],'delta',(inFileInfo.(struName).delta)');
    
    h5write(outFileName,['/',fldName],(inFileInfo.(struName).angular_velocity)',[3,1],[1,2]);
    h5writeatt(outFileName,['/',fldName],'angular velocity',(inFileInfo.(struName).angular_velocity)');
    
    h5write(outFileName,['/',fldName],(inFileInfo.(struName).phi0)',[4,1],[1,2]);
    h5writeatt(outFileName,['/',fldName],'phi0',(inFileInfo.(struName).phi0)');
    
    h5write(outFileName,['/',fldName],(inFileInfo.(struName).amplitude)',[5,1],[1,2]);
    h5writeatt(outFileName,['/',fldName],'amplitude',(inFileInfo.(struName).amplitude)');
    
    h5write(outFileName,['/',fldName],(inFileInfo.(struName).inclination)',[6,1],[1,2]);
    h5writeatt(outFileName,['/',fldName],'inclination',(inFileInfo.(struName).inclination)');
    
    h5write(outFileName,['/',fldName],(inFileInfo.FreqRange)',[7,1],[1,2]);
    h5writeatt(outFileName,['/',fldName],'Frequency range',(inFileInfo.FreqRange)');
end