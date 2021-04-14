% Generate Lonestar5 job file
% Multiple realizations version
% Author: QYQ
% 1/5/2021
%% Setting up Path
clear;
simDataDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/simData/diff_srlz';
DataDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/simData/diff_srlz'; % local
%simDataDir = '/work/05884/qyqstc/lonestar/GWBsimDataSKA_HDF5'; % on ls5
search_paramDataDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/searchParams/2bands/superNarrow';
% folderName = dir([search_paramDataDir,filesep,'*GWB*']);
% folderName = sort_nat({folderName.name});

NsimFiles = dir([simDataDir,filesep,'GWBsimDataSKASrlz*.mat']);
Nreal = length(NsimFiles); % # of realizations

%% LS5 config
excut = '/work/05884/qyqstc/lonestar/RAAPTR/MxAvPhaseC/Multi_PSO.out ';
searchParamDir = '/work/05884/qyqstc/lonestar/MultiPSO/Task8/searchParams/2bands/superNarrow/'; % end with /
dataDir = '/work/05884/qyqstc/lonestar/MultiPSO/Final/realizations/2bands/simData/diff_srlz/'; % end with /
resultsDir = [dataDir,'results/'];
ite = 20;


fid = fopen('/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/JOBFILE_Diff_Srlz.txt','w');% 'a' for appending data at the end of the file.

%% Write file
for r = 1:4 %Nreal
    %     inFileList_para =
    %     dir([search_paramDataDir,filesep,folderName{r},filesep,'*.hdf5']); %
    %     uncomment if each realization has its own searchParams
    inFileList_para = dir([search_paramDataDir,filesep,'*.hdf5']);
    ParaNames = sort_nat({inFileList_para.name});
    exp = 'searchParams\d.hdf5';
    ParaNames = regexp(ParaNames,exp,'match');
    ParaNames = ParaNames(~cellfun(@isempty,ParaNames));
    nParamFiles = length(ParaNames);
    inFileList = dir([DataDir,filesep,'GWB*','.hdf5']);% original file
    cpyFileList = dir([DataDir,filesep,'HDF5',filesep,'*GWB*','.hdf5']);% total files including duplicate files for different bands.
    cpyFiles = length(cpyFileList); % Multiple bands
    
    if isempty(inFileList) == 1
        nFiles = length(cpyFileList);% Multiple band
        inFileName = sort_nat({cpyFileList.name});
        
    else
        nFiles = length(inFileList); % Single band
        inFileName = sort_nat({inFileList.name});
    end
    
    %     bandFiles = nFiles/nParamFiles;
    
    for lpc = 1:nParamFiles
        
        %             bandNum = h5readatt([search_paramDataDir,'/',char(ParaNames{lpc})],...
        %                 ['/','Bands_Info'],'bandNumber');
        fprintf(fid,excut);
        %         fprintf(fid,'%s%s/%s ',searchParamDir, folderName{r},
        %         ParaNames{lpc}{1}); % uncomment if each realization has different
        %         searchParams.
        fprintf(fid,'%s/%s ',searchParamDir, ParaNames{lpc}{1});
        % select files for corresponding realizations
        %         exp = ['\d+_',folderName{r},'(?=_|\.hdf5)_?\d{0,2}.hdf5'];
%         exp = ['\d+_','(?=_|\.hdf5)_?\d{0,2}.hdf5']; % for same srlz but diffrent nrlz.
        exp = ['\d+_GWBsimDataSKASrlz',num2str(r),'.*\.hdf5'];
        FilenameList = regexp(inFileName,exp,'match');
        FilenameList = FilenameList(~cellfun(@isempty,FilenameList));
        fprintf(fid,dataDir + "%s ",FilenameList{lpc}{1});
        fprintf(fid,resultsDir + "%s ",FilenameList{lpc}{1});
        fprintf(fid,'avPhase ');
        fprintf(fid,'%d',ite); % number of iteration
        fprintf(fid,'\n'); % if this line is the last line delete the blank line, it will affect Launcher counting the number of jobs.
        
    end
    
end

fclose(fid);

% END