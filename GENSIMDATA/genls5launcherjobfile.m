clear;
simDataDir = '~/Research/PulsarTiming/SimDATA/Mauritius/GWBsimDataSKA_HDF5'; % local
%simDataDir = '/work/05884/qyqstc/lonestar/GWBsimDataSKA_HDF5'; % on ls5
search_paramDataDir = '~/Research/PulsarTiming/SimDATA/Mauritius/uniform search params HDF5';
inFileList = dir([simDataDir,filesep,'*.hdf5']);
inFileList_param = dir([search_paramDataDir,filesep,'*.hdf5']);
nParamFiles = length(inFileList_param);
nFiles = length(inFileList);
fid = fopen('Maricius_UNI_JOBFILE.txt','w');

%% single source
for lpc = 1:nFiles
    fprintf(fid,'/work/05884/qyqstc/lonestar/PULSARTIMING/RAAPTR/MxAvPhaseC/perfeval_spmd.out ');
    % fprintf(fid,'/work/05884/qyqstc/lonestar/PULSARTIMING/RAAPTR/MxAvPhaseC/searchParams_SimDATASKA.hdf5 ');
    fprintf(fid,'/work/05884/qyqstc/lonestar/searchParams_HDF5/%s ', inFileList(lpc).name);
    fprintf(fid,'/work/05884/qyqstc/lonestar/TESTDATA_HDF5/GWBsimDataSKASrlz1Nrlz1.hdf5 ');
    % fprintf(fid,'/Users/qianyiqian/Research/PulsarTiming/SimDATA/TESTDATA_HDF5/%s ',inFileList(lpc).name);
    % fprintf(fid,'/Users/qianyiqian/Research/PulsarTiming/SimDATA/Out/%s ',inFileList(lpc).name);
    fprintf(fid,'/work/05884/qyqstc/lonestar/Multi_Results/%s ', inFileList(lpc).name);
    fprintf(fid,'maxPhase');
    fprintf(fid,'\n');
end
fclose(fid);

%% multiple source
for lpc = 1:nFiles
    for ppc = 1:nParamFiles
        bandNum = h5readatt([search_paramDataDir,'/',inFileList_param(ppc).name],...
            ['/','Bands_Info'],'bandNumber');
        fprintf(fid,'/work/05884/qyqstc/lonestar/PULSARTIMING/RAAPTR/MxAvPhaseC/perfeval_spmd.out ');
        fprintf(fid,'/work/05884/qyqstc/lonestar/uniform_searchParams_HDF5/%s ',inFileList_param(ppc).name);
        fprintf(fid,'/work/05884/qyqstc/lonestar/GWBsimDataSKA_HDF5/%s ',inFileList(lpc).name);
        fprintf(fid,'/work/05884/qyqstc/lonestar/Mauritius_uni_results/%d_%s ',bandNum,inFileList(lpc).name);
        fprintf(fid,'avPhase');
        fprintf(fid,'\n');
    end
end
fclose(fid);