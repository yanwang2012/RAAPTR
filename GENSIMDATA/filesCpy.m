%% Copy input files for different band use
simDataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test10/FullBand/WLSRC-band1/HDF5';
inFiles = dir([simDataDir,filesep,'*rm*.hdf5']);
NumBand = 5;
N = length(inFiles);
for i = 1:N
    for j = 1:NumBand
        copyfile([simDataDir,filesep,inFiles(i).name],[simDataDir,filesep,num2str(j) '_' inFiles(i).name])
    end
end