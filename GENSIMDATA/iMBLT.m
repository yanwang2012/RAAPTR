% MBLT in ONE band
% QYQ 02/01/2019
clear;
tic
%% Set up
simParamsDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/searchParams/2bands';
simParamsName = 'searchParams';
inParamsList = dir([simParamsDir,filesep,simParamsName,'*.mat']);
simDataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands/ONE/band1';
estDataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands/ONE/band1/Results10';
inputFileName = 'GWBsimDataSKASrlz1Nrlz3';
bandNum = 1;
outputfiles = dir([estDataDir,filesep,num2str(bandNum),'_',inputFileName,'*.mat']);
Npara = length(inParamsList);
NestSrc = length(outputfiles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DON'T FORGET TO CHECK THE NAME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nFile = dir([estDataDir,filesep,inputFileName,'*.mat']); % count how many iterations are used.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_ite = length(nFile);
% Load the simulated source parameters.
load([simDataDir,filesep,inputFileName,'.mat']);
estTimRes = zeros(simParams.Np,simParams.N);
ResCell = {}; % a cell of timing residuals to store all the estimated residuals.
SNRarray = [];

%% MBLT
[file,Index]=rassign(estDataDir,outputfiles,NestSrc,simParams,yr);
% disp(["File needs to be skipped: ",file]);
outputFilename = 'GWBsimDataSKASrlz1Nrlz3';
OutputDir = [simDataDir,filesep,outputFilename];
mkdir(OutputDir);

for j = 1:NestSrc
    
    if j == Index
        continue
    else
        %                 disp("j is:"+j);
        path_estData = [estDataDir,filesep,outputfiles(j).name];
        [srcParams]=ColSrcParams(path_estData);
        [SNR,estTimRes_tmp] = Amp2Snr(srcParams,simParams,yr);
       ResCell = [ResCell estTimRes_tmp];
       SNRarray = [SNRarray SNR];
    end
    
end

for i = 1:NestSrc-1
    estTimRes = estTimRes + ResCell{i+1};
end

newFile = strcat(OutputDir,filesep,inputFileName,'band ',num2str(bandNum),'Only','.mat');
copyfile([simDataDir,filesep,inputFileName,'.mat'],newFile);
m = matfile(newFile,'Writable',true);
m.timingResiduals = timingResiduals - estTimRes;



toc
% END