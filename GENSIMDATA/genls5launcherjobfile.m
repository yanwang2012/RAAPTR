inFileList = dir(['/Users/qianyiqian/Research/RAAPTR/GENSIMDATA/TESTDATA_HDF5',filesep,'*.hdf5']);
nFiles = length(inFileList);
fid = fopen('GWBsimDataSKA_JOBFILE.txt','w');
for lpc = 1:nFiles
    fprintf(fid,'/work/05884/qyqstc/lonestar/PULSARTIMING/RAAPTR/MxAvPhaseC/perfeval_spmd.out ');
    fprintf(fid,'/work/05884/qyqstc/lonestar/PULSARTIMING/RAAPTR/MxAvPhaseC/searchParams_simDataSKA.hdf5 ');
    fprintf(fid,'/work/05884/qyqstc/lonestar/TESTDATA_HDF5/%s ', inFileList(lpc).name);
    fprintf(fid, '/work/05884/qyqstc/lonestar/RESULTS/%s ',inFileList(lpc).name);
    fprintf(fid,'MaxPhase');
    fprintf(fid,'\n');
end
fclose(fid);