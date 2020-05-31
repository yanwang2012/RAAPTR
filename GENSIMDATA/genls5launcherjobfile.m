clear;
simDataDir = '/Users/qianyiqian/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/HDF5'; % local
%simDataDir = '/work/05884/qyqstc/lonestar/GWBsimDataSKA_HDF5'; % on ls5
search_paramDataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/searchParams/2bands/HDF5';
inFileList = dir([simDataDir,filesep,'GWB*','.hdf5']);% original file
cpyFileList = dir([simDataDir,filesep,'*_GWB*','.hdf5']);% total files including duplicate files for different bands.
inFileList_param = dir([search_paramDataDir,filesep,'*','.hdf5']);
nParamFiles = length(inFileList_param);
cpyFiles = length(cpyFileList); % Multiple bands
if isempty(inFileList) == 1
    nFiles = length(cpyFileList);% Multiple band
    inFileName = sort_nat({cpyFileList.name});

else
    nFiles = length(inFileList); % Single band
    inFileName = sort_nat({inFileList.name}); 
end
fid = fopen('~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/2bandRealizations.txt','w');% 'a' for appending data at the end of the file.

%% LS5 config
excut = '/work/05884/qyqstc/lonestar/RAAPTR/MxAvPhaseC/Multi_PSO.out ';
searchParamDir = '/work/05884/qyqstc/lonestar/MultiPSO/Task8/searchParams/2bands/HDF5/';
dataDir = '/work/05884/qyqstc/lonestar/MultiPSO/Final/realizations/2bands/';
resultsDir = [dataDir,'results/'];
ite = 20;

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
% inFileName = {};
% for k = 1:cpyFiles
%     inFileName = [inFileName cpyFileList(k).name];
% end

bandFiles = nFiles/nParamFiles;

for lpc = 1:nParamFiles
    for ppc = 1:bandFiles
        bandNum = h5readatt([search_paramDataDir,'/',inFileList_param(lpc).name],...
            ['/','Bands_Info'],'bandNumber');
        fprintf(fid,excut);
        fprintf(fid,searchParamDir + "%s ",inFileList_param(lpc).name); % double quotes create string, "plus" combine strings.
        fprintf(fid,dataDir + "%s ",char(inFileName((lpc-1) * bandFiles + ppc)));
        fprintf(fid,resultsDir + "%s ",char(inFileName((lpc-1) * bandFiles + ppc)));
        fprintf(fid,'avPhase ');
        fprintf(fid,'%d',ite); % number of iteration
        fprintf(fid,'\n'); % if this line is the last line delete the blank line, it will affect Launcher counting the number of jobs.
    end
end

fclose(fid);

% EOF