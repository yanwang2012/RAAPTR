% A script convert searchParams file from .mat to .hdf5 from different
% directoies.
tic
clear;
FolderDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/searchParams/bandTrials';
FolderName = dir(FolderDir);
FolderName = FolderName(~startsWith({FolderName.name},'.')); % ignore folders start with '.'
FolderName = sort_nat({FolderName.name}); % sort the folder names naturally
N = length(FolderName);

for i = 1:N
    OutDir = [FolderDir,filesep,char(FolderName(i)),filesep,'HDF5'];
    inFileList = dir([FolderDir,filesep,char(FolderName(i)),filesep,'*.mat']);
    mkdir(OutDir);
    for lpc = 1:length(inFileList)
        inFile = [FolderDir,filesep,char(FolderName(i)),filesep,inFileList(lpc).name];
        [~,inFileName,~] = fileparts(inFile);
        outFileName = [OutDir,filesep,inFileName,'.hdf5'];
        inFileInfo = load(inFile);
        struName = 'searchParams';
        fldName = 'xmaxmin';
        bandsIn = 'Bands_Info';
        %Create HDF5 file
        fid = H5F.create(outFileName);
        H5F.close(fid);
        nDim = length(size((inFileInfo.(struName).alpha)'));
        baseSz = ones(1, nDim);
        
        h5create(outFileName,['/',bandsIn],[2,1]); %[col, row]
        h5write(outFileName,['/',bandsIn],inFileInfo.(struName).band_num,[1,1],[1,1])%[start point],[size of data].
        h5writeatt(outFileName,['/',bandsIn],'bandNumber',inFileInfo.(struName).band_num);
        h5write(outFileName,['/',bandsIn],inFileInfo.NumBands,[2,1],[1,1]);
        h5writeatt(outFileName,['/',bandsIn],'NumBands',inFileInfo.NumBands);
        
        %     h5create(outFileName,['/',fldName],...
        %         max([size((inFileInfo.(struName).alpha)'); baseSz]));
        h5create(outFileName,['/',fldName],[2,7])% [2,7] 2 colums 7 rows
        h5write(outFileName,['/',fldName],(inFileInfo.(struName).alpha),[1,1],[2,1]); % [2,1] 2 rows 1 colum but writing horizontally into that dataset.
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
end


toc
% END