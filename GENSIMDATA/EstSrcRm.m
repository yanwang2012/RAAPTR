%% Remove estimated sources
% Given a set of output files, remove all the estimated sources from
% simulation data.
% QYQ 03/11/2019

clear;
tic
%% Set up
simDataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands';
estDataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands/supperNarrow_iMBLT1/Results_5';
inputFileName = 'GWBsimDataSKASrlz1Nrlz3';
ext = '.mat';
outputfiles = dir([estDataDir,filesep,'*',inputFileName,'*',ext]);
% Npara = length(inParamsList);
NestSrc = length(outputfiles);

% Load the simulated source parameters.
load([simDataDir,filesep,inputFileName,'.mat']);
estTimRes = zeros(simParams.Np,simParams.N);
ResCell = {}; % a cell of timing residuals to store all the estimated residuals.
SNRarray = [];

%% MBLT
[file,Index]=rassign(estDataDir,outputfiles,NestSrc,simParams,yr);
% disp(["File needs to be skipped: ",file]);
foldername = 'supperNarrow_iMBLT1_after';
OutputDir = [simDataDir,filesep,foldername];
outputfilenames = sort_nat({outputfiles.name});

mkdir(OutputDir);

for j = 1:NestSrc
    
    if j == Index
        continue
    else
        %                 disp("j is:"+j);
        path_estData = [estDataDir,filesep,char(outputfilenames(j))];
        [srcParams]=ColSrcParams(path_estData);
        [SNR,estTimRes_tmp] = Amp2Snr(srcParams,simParams,yr);
       ResCell = [ResCell estTimRes_tmp];
       SNRarray = [SNRarray SNR];
    end
    
end

% Accumulate the timing residulas for different sources.
for nsrc = 2:NestSrc
    estTimRes = estTimRes + ResCell{nsrc};
end

newFile = strcat(OutputDir,filesep,inputFileName,'_After_iMBLT_1','.mat');
copyfile([simDataDir,filesep,inputFileName,'.mat'],newFile);
m = matfile(newFile,'Writable',true);
m.timingResiduals = timingResiduals - estTimRes;



toc
% END