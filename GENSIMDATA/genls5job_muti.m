% Generate Lonestar5 job file
% Multiple realizations version
% Author: QYQ
% 1/5/2021
%% Setting up Path
clear;
DataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands'; % local
%simDataDir = '/work/05884/qyqstc/lonestar/GWBsimDataSKA_HDF5'; % on ls5
search_paramDataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/searchParams/2bands/HDF5';
folderName = dir([DataDir,filesep,'*xMBLT']);
folderName = sort_nat({folderName.name});
inFileList_para = dir([search_paramDataDir,filesep,'*','.hdf5']);
ParaNames = sort_nat({inFileList_para.name});
exp = 'searchParams\d.hdf5';
ParaNames = regexp(ParaNames,exp,'match');
ParaNames = ParaNames(~cellfun(@isempty,ParaNames));
nParamFiles = length(ParaNames);

fid = fopen('~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/Realizations_xMBLT.txt','w');% 'a' for appending data at the end of the file.

Nreal = 50; % # of realizations

%% LS5 config
excut = '/work/05884/qyqstc/lonestar/RAAPTR/MxAvPhaseC/Multi_PSO.out ';
searchParamDir = '/work/05884/qyqstc/lonestar/MultiPSO/Task8/searchParams/2bands/HDF5/';
dataDir = '/work/05884/qyqstc/lonestar/MultiPSO/Final/realizations/2bands/xMBLT/';
resultsDir = [dataDir,'results/'];
ite = 20;



%% Write file
for r = 1:Nreal
    inFileList = dir([DataDir,filesep,folderName{r},filesep,'GWB*','.hdf5']);% original file
    cpyFileList = dir([DataDir,filesep,folderName{r},filesep,'HDF5',filesep,'*GWB*','.hdf5']);% total files including duplicate files for different bands.
    cpyFiles = length(cpyFileList); % Multiple bands
    
    if isempty(inFileList) == 1
        nFiles = length(cpyFileList);% Multiple band
        inFileName = sort_nat({cpyFileList.name});
        
    else
        nFiles = length(inFileList); % Single band
        inFileName = sort_nat({inFileList.name});
    end
    
    bandFiles = nFiles/nParamFiles;
    
    for lpc = 1:nParamFiles
        for ppc = 1:bandFiles
            bandNum = h5readatt([search_paramDataDir,'/',char(ParaNames{lpc})],...
                ['/','Bands_Info'],'bandNumber');
            fprintf(fid,excut);
            fprintf(fid,searchParamDir + "%s ",char(ParaNames{lpc})); % double quotes create string, "plus" combine strings.
            fprintf(fid,dataDir + "%s ",[folderName{r},filesep,char(inFileName((lpc-1) * bandFiles + ppc))]);
            fprintf(fid,resultsDir + "%s ",[folderName{r},filesep,char(inFileName((lpc-1) * bandFiles + ppc))]);
            fprintf(fid,'avPhase ');
            fprintf(fid,'%d',ite); % number of iteration
            fprintf(fid,'\n'); % if this line is the last line delete the blank line, it will affect Launcher counting the number of jobs.
        end
    end
    
end

fclose(fid);

% END