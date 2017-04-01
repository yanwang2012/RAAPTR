%% Reconfigure data file directory for run restarts
%Assumed data repository structure: recnfDir has *.mat input data files.
%Results from runs are stored in folder with resultsDirPattern* wildcard
%matched directories. Results are in the form of *.mat files with exactly
%the same name as the corresponding input files. 
recnfDir = 'simDataSKA_loc7';
resultsDirPattern = 'results_';

%Directory on which to rerun
newDir = [recnfDir,'_restart_1'];

%Number of workers to use
nWorkers = 8;

%List all results directories
resultsDirs = dir([recnfDir,filesep,resultsDirPattern,'*']);
nResultsDirs = length(resultsDirs);

%Record all input files already processed
nResultsFiles = zeros(nResultsDirs,1);
for lp = 1:nResultsDirs
    nResultsFiles(lp) = length(dir([recnfDir,filesep,resultsDirs(lp).name, filesep, '*.mat']));
end
namesResultFiles = cell(sum(nResultsFiles),1);
countFiles = 1;
for lp = 1:nResultsDirs
    fileList = dir([recnfDir,filesep,resultsDirs(lp).name, filesep, '*.mat']);
    for lp2 = 1:nResultsFiles(lp)
        namesResultFiles{countFiles} = fileList(lp2).name;
        countFiles = countFiles+1;
    end
end
countFiles = countFiles -1;

%Get list of all input files
inputFileList = dir([recnfDir,filesep,'*.mat']);
nInputFiles = length(inputFileList);
namesInputFiles = cell(nInputFiles, 1);
for lp = 1:nInputFiles
    namesInputFiles{lp} = inputFileList(lp).name;
end

%Remove names of result files from list of input files
indxKeepFiles = ones(nInputFiles,1);
for lp = 1:nInputFiles
    inputFileName = namesInputFiles{lp};
    for lp2 = 1:countFiles
        if strcmp(inputFileName,namesResultFiles{lp2})
            indxKeepFiles(lp) = 0;
            break;
        end
    end
end

%Create restart directory and copy input files
mkdir(newDir);
copyFileCount = 0;
for lp = 1:nInputFiles
    if(indxKeepFiles(lp))
        copyfile([recnfDir,filesep,namesInputFiles{lp}],[newDir,filesep,namesInputFiles{lp}]);
        copyFileCount = copyFileCount +1;
    end
end

%Number of files per worker
nfw = nfilesperworker(nWorkers, copyFileCount);
fprintf(1,'%d %d\n',nfw');