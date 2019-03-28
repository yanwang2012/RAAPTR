simDataDir = '~/Research/PulsarTiming/SimDATA/Maricius results/GWBsimDataSKA_HDF5'; % local
%simDataDir = '/work/05884/qyqstc/lonestar/GWBsimDataSKA_HDF5'; % on ls5
search_paramDataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/searchParams_HDF5';
inFileList = dir([simDataDir,filesep,'*.hdf5']);
inFileList_param = dir([search_paramDataDir,filesep,'*.hdf5']);
nParamFiles = length(inFileList_param);
nFiles = length(inFileList);
fid = fopen('Maricius_JOBFILE.txt','w');

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
        fprintf(fid,'/work/05884/qyqstc/lonestar/PULSARTIMING/RAAPTR/MxAvPhaseC/perfeval_spmd.out ');
        fprintf(fid,'/work/05884/qyqstc/lonestar/searchParams_HDF5/%s ',inFileList_param(ppc).name);
        fprintf(fid,'/work/05884/qyqstc/lonestar/GWBsimDataSKA_HDF5/%s ',inFileList(lpc).name);
        fprintf(fid,'/work/05884/qyqstc/lonestar/Maricius_results/%d_%s ',ppc,inFileList(lpc).name);
        fprintf(fid,'AvPhase');
        fprintf(fid,'\n');
    end
end
fclose(fid);