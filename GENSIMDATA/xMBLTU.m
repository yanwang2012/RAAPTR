% xMBLT Implementation for united est. sources

% Author: QYQ 09/16/2020

clear;
tic
%% Set up

simParamsDir = '/Users/qyq/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/searchParams/2bands/superNarrow';
simParamsName = 'searchParams';
inParamsList = dir([simParamsDir,filesep,simParamsName,'*.mat']);
simDataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands';
estDataDir = '/Users/qyq/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands/SuperNarrow/Results_supNar_rand1/GWBsimDataSKASrlz1Nrlz3_xMBLT/results/Union2';
inputFileName = 'GWBsimDataSKASrlz1Nrlz3';

inParamNames = sort_nat({inParamsList.name});
exp = 'searchParams\d.mat'; % regular expressions for desire file names
inParamNames = regexp(inParamNames,exp,'match');
inParamNames = inParamNames(~cellfun(@isempty,inParamNames)); % get rid of empty cells
Npara = length(inParamNames);

% repartition united sources according to searchParams file
outNames = {};
outfiles1 = dir([estDataDir,filesep,'1_*.mat']);
outNames{1} = sort_nat({outfiles1.name});
outfiles2 = dir([estDataDir,filesep,'2_*.mat']);
outNames{2} = sort_nat({outfiles2.name});
outfilesAll = dir([estDataDir,filesep,'*',inputFileName,'*.mat']);
Allfilesname = sort_nat({outfilesAll.name});
AllSrc = length(Allfilesname);
Nfiles = zeros(Npara,1);% number of output files
Nfiles(1) = length(outNames{1});
Nfiles(2) = length(outNames{2});

% check if band n sources above band n edge, if so change filenames to
% (n+1)_INPUTFILENAME*.mat. * Note *: Now only works for 2 bands case.
for i = 1:Npara
    load([simParamsDir,filesep,char(inParamNames{i})],'searchParams');
    sim_angular_velocity = searchParams.angular_velocity;
    for j = 1:Nfiles(i)
        load([estDataDir,filesep,char(outNames{i}(j))],'bestRealLoc');
        est_angular_velocity = bestRealLoc(3);
        switch i
            case 1
                if est_angular_velocity > sim_angular_velocity(i)
                    movefile([estDataDir,filesep,char(outNames{i}(j))],[estDataDir,filesep,num2str(swap(i)),'_',char(outNames{i}(j))]);
                end
            case 2
                if est_angular_velocity <= sim_angular_velocity(i)
                    movefile([estDataDir,filesep,char(outNames{i}(j))],[estDataDir,filesep,num2str(swap(i)),'_',char(outNames{i}(j))]);
                end
        end
    end
end

%% xMBLT
% Load the simulated source parameters.
load([simDataDir,filesep,inputFileName,'.mat']);
estTimRes = zeros(simParams.Np,simParams.N);

[file,Index]=rassign(estDataDir,Allfilesname,AllSrc,Nfiles(1),simParams,yr);
Filename = 'xMBLT';
OutputDir = [estDataDir,filesep,Filename];
mkdir(OutputDir);
for i = 1:Npara
    outputfiles = dir([estDataDir,filesep,num2str(swap(i)),'_*',inputFileName,'*.mat']);
    NestSrc = length(outputfiles);
    outputfilenames = sort_nat({outputfiles.name});
%     nFile = dir([estDataDir,filesep,num2str(swap(i)),'_*',inputFileName,'*.mat']); % count how many sources in each band.
%         num_ite = length(nFile);
    for j = 1:NestSrc
        if ~isempty(file) && strcmp(char(outNames{swap(i)}(j)),char(file(1))) == 1
            continue
        else
            %                 disp("j is:"+j);
            path_estData = [estDataDir,filesep,char(outputfilenames(j))];
            disp(['File loaded: ',char(outputfilenames(j))]);
            [srcParams]=ColSrcParams(path_estData);
            [~,estTimRes_tmp] = Amp2Snr(srcParams,simParams,yr);
            estTimRes = estTimRes + estTimRes_tmp;
        end
    end
    newFile = strcat(OutputDir,filesep,num2str(i),'_',inputFileName,'.mat');
    copyfile([simDataDir,filesep,inputFileName,'.mat'],newFile);
    m = matfile(newFile,'Writable',true);
    m.timingResiduals = timingResiduals - estTimRes;
    estTimRes = zeros(simParams.Np,simParams.N);
end


toc
% END