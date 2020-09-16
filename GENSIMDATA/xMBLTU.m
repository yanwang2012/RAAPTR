% xMBLT Implementation for united est. sources

% Author: QYQ 09/16/2020

clear;
tic
%% Set up
simParamsDir = '/Users/qyq/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/searchParams/2bands/superNarrow';
simParamsName = 'searchParamsRand';
inParamsList = dir([simParamsDir,filesep,simParamsName,'*.mat']);
simDataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands';
estDataDir = '/Users/qyq/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands/Results_supNar_rand1/GWBsimDataSKASrlz1Nrlz3_xMBLT/results/Union';
inputFileName = 'GWBsimDataSKASrlz1Nrlz3';
Npara = length(inParamsList);

% Load the simulated source parameters.
load([simDataDir,filesep,inputFileName,'.mat']);
estTimRes = zeros(simParams.Np,simParams.N);

%% xMBLT
Filename = 'GWBsimDataSKASrlz1Nrlz3_xMBLT_Union';
OutputDir = [estDataDir,filesep,Filename];
mkdir(OutputDir);
for i = 1:Npara
    outputfiles = dir([estDataDir,filesep,num2str(swap(i)),'_',inputFileName,'*.mat']);
    NestSrc = length(outputfiles);
    outputfilenames = sort_nat({outputfiles.name});
    [file,Index]=rassign(estDataDir,outputfiles,NestSrc,simParams,yr);
    nFile = dir([estDataDir,filesep,num2str(swap(i)),'_',inputFileName,'*.mat']); % count how many sources in each band.
%     num_ite = length(nFile);
    for j = 1:NestSrc
            if ismember(j,Index) == 1
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