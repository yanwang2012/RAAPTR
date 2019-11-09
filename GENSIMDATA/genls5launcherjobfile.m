clear;
simDataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test10/FullBand/WLSRC-band1/HDF5'; % local
%simDataDir = '/work/05884/qyqstc/lonestar/GWBsimDataSKA_HDF5'; % on ls5
search_paramDataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test10/searchParams/HDF5';
inFileList = dir([simDataDir,filesep,'GWB*rm*.hdf5']);% original file
cpyFileList = dir([simDataDir,filesep,'*_GWB*_rm*.hdf5']);% file copies
inFileList_param = dir([search_paramDataDir,filesep,'*.hdf5']);
nParamFiles = length(inFileList_param);
cpyFiles = length(cpyFileList); % including bands
nFiles = length(inFileList);
fid = fopen('~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test10/JOBFILERMBAND1.txt','w');

%% single source
% for lpc = 1:nFiles
%     fprintf(fid,'/work/05884/qyqstc/lonestar/PULSARTIMING/RAAPTR/MxAvPhaseC/perfeval_spmd.out ');
%     fprintf(fid,'/work/05884/qyqstc/lonestar/PULSARTIMING/RAAPTR/MxAvPhaseC/searchParams_SimDATASKA.hdf5 ');
%     % fprintf(fid,'/work/05884/qyqstc/lonestar/searchParams_HDF5/%s ', inFileList(lpc).name);
%     % fprintf(fid,'/work/05884/qyqstc/lonestar/TESTDATA_HDF5/GWBsimDataSKASrlz1Nrlz1.hdf5 ');
%     fprintf(fid,'/Users/qianyiqian/Research/PulsarTiming/SimDATA/TESTDATA_HDF5/%s ',inFileList(lpc).name);
%     % fprintf(fid,'/Users/qianyiqian/Research/PulsarTiming/SimDATA/Out/%s ',inFileList(lpc).name);
%     fprintf(fid,'/work/05884/qyqstc/lonestar/Multi_Results/%s ', inFileList(lpc).name);
%     fprintf(fid,'maxPhase');
%     fprintf(fid,'\n');
% end
% fclose(fid);

%% multiple source
% sort the input file lists' name naturally
inFileName = {};
excut = '/work/05884/qyqstc/lonestar/RAAPTR/MxAvPhaseC/Multi_PSO.out ';
searchParamDir = '/work/05884/qyqstc/lonestar/MultiPSO/Task7/searchParams/';
dataDir = '/work/05884/qyqstc/lonestar/MultiPSO/Task7/FullBand/WLSRC-band1/';
resultsDir = '/work/05884/qyqstc/lonestar/MultiPSO/Task7/FullBand/WLSRC-band1/Results/';
for k = 1:cpyFiles
    inFileName = [inFileName cpyFileList(k).name];
end
inFileName = sort_nat(inFileName);
for lpc = 1:nParamFiles
    for ppc = 1:nFiles
        bandNum = h5readatt([search_paramDataDir,'/',inFileList_param(lpc).name],...
            ['/','Bands_Info'],'bandNumber');
        fprintf(fid,excut);
        fprintf(fid,searchParamDir + "%s ",inFileList_param(lpc).name); % double quotes create string, "plus" combine strings.
        fprintf(fid,dataDir + "%s ",char(inFileName((lpc-1)*nFiles+ppc)));
        fprintf(fid,resultsDir + "%s ",char(inFileName((lpc-1)*nFiles+ppc)));
        fprintf(fid,'avPhase ');
        fprintf(fid,'%d',10); % number of iteration
        fprintf(fid,'\n'); % if this line is the last line delete the blank line, it will affect Launcher counting the number of jobs.
    end
end
fclose(fid);