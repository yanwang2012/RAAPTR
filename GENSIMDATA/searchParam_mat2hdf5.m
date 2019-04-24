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
    h5create(outFileName,['/',fldName],[2,7])% fuck it hdf5 [colum row] but wright horizentally
    h5write(outFileName,['/',fldName],(inFileInfo.(struName).alpha),[1,1],[2,1]);
    h5writeatt(outFileName,['/',fldName],'alpha',(inFileInfo.(struName).alpha)');
    
    h5write(outFileName,['/',fldName],(inFileInfo.(struName).delta),[1,2],[2,1]);
    h5writeatt(outFileName,['/',fldName],'delta',(inFileInfo.(struName).delta)');
    
    h5write(outFileName,['/',fldName],(inFileInfo.(struName).angular_velocity),[1,3],[2,1]);
    h5writeatt(outFileName,['/',fldName],'angular velocity',(inFileInfo.(struName).angular_velocity)');
    
    h5write(outFileName,['/',fldName],(inFileInfo.(struName).phi0),[1,4],[2,1]);
    h5writeatt(outFileName,['/',fldName],'phi0',(inFileInfo.(struName).phi0)');
    
    h5write(outFileName,['/',fldName],(inFileInfo.(struName).amplitude),[1,5],[2,1]);
    h5writeatt(outFileName,['/',fldName],'amplitude',(inFileInfo.(struName).amplitude)');
    
    h5write(outFileName,['/',fldName],(inFileInfo.(struName).inclination),[1,6],[2,1]);
    h5writeatt(outFileName,['/',fldName],'inclination',(inFileInfo.(struName).inclination)');
    
    h5write(outFileName,['/',fldName],(inFileInfo.(struName).polarization),[1,7],[2,1]);
    h5writeatt(outFileName,['/',fldName],'polarization',(inFileInfo.(struName).polarization)');
    
    %h5write(outFileName,['/',fldName],(inFileInfo.FreqRange),[1,7],[2,1]);
    h5writeatt(outFileName,['/',fldName],'Frequency range',(inFileInfo.FreqRange)');
end