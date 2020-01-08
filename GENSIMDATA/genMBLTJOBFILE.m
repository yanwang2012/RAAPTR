clear;
simDataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/Realizations/GWBsimDataSKASrlz1Nrlz5/HDF5'; % local
%simDataDir = '/work/05884/qyqstc/lonestar/GWBsimDataSKA_HDF5'; % on ls5
search_paramDataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test10/searchParams/HDF5';
inFileList = dir([simDataDir,filesep,'GWB*.hdf5']);% original file
inFileList_param = dir([search_paramDataDir,filesep,'*.hdf5']);
nParamFiles = length(inFileList_param);
nFiles = length(inFileList);
fid = fopen('~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/MBLTJOBFILE_REA2.txt','a'); % 'a' for appending data at the end of the file.
excut = '/work/05884/qyqstc/lonestar/RAAPTR/MxAvPhaseC/Multi_PSO.out ';
searchParamDir = '/work/05884/qyqstc/lonestar/MultiPSO/Task7/searchParams/';
dataDir = '/work/05884/qyqstc/lonestar/MultiPSO/Task8/Realizations/GWBsimDataSKASrlz1Nrlz5/';
resultsDir = '/work/05884/qyqstc/lonestar/MultiPSO/Task8/Realizations/GWBsimDataSKASrlz1Nrlz5/Results/';

inFileName = {};
for k = 1:nFiles
    inFileName = [inFileName inFileList(k).name];
end
inFileName = sort_nat(inFileName);

for ppc = 1:nFiles
    bandNum = h5readatt([search_paramDataDir,'/',inFileList_param(ppc).name],...
        ['/','Bands_Info'],'bandNumber');
    fprintf(fid,excut);
    fprintf(fid,searchParamDir + "%s ",inFileList_param(ppc).name); % double quotes create string, "plus" combine strings.
    fprintf(fid,dataDir + "%s ",char(inFileName(ppc)));
    fprintf(fid,resultsDir + "%s ",char(inFileName(ppc)));
    fprintf(fid,'avPhase ');
    fprintf(fid,'%d',10); % number of iteration
    fprintf(fid,'\n'); % if this line is the last line delete the blank line, it will affect Launcher counting the number of jobs.
end