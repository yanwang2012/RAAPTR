%% Remove estimated sources
% Given a set of output files, remove all the estimated sources from
% simulation data.
% QYQ 03/11/2019

clear;
tic
%% Set up
simDataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands/supperNarrow_iMBLT1_after_20/Results_20/2_iMBLT/results/2iMBLT_after/results/3_iMBLT/results/3iMBLT_after/results/4_iMBLT/results/4iMBLT_after/results/5_iMBLT/results/5iMBLT_after/results/6_iMBLT/results/6iMBLT_after/';
estDataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands/supperNarrow_iMBLT1_after_20/Results_20/2_iMBLT/results/2iMBLT_after/results/3_iMBLT/results/3iMBLT_after/results/4_iMBLT/results/4iMBLT_after/results/5_iMBLT/results/5iMBLT_after/results/6_iMBLT/results/6iMBLT_after/results/7_iMBLT/results';
inputFileName = 'GWBsimDataSKASrlz1Nrlz3';
ext = '.mat';
outputfiles = dir([estDataDir,filesep,'*',inputFileName,'*',ext]);
% Npara = length(inParamsList);
NestSrc = length(outputfiles);
stage = 7; % stage number of iMBLT

% Load the simulated source parameters.
if stage == 1
    simDataDir = '/Users/qianyiqian/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/';
    load([simDataDir,filesep,inputFileName,ext])
else
    load([simDataDir,filesep,inputFileName,'_after_',num2str(stage-1),'iMBLT',ext]);% subtraction in iMBLT needs to base on the new inputdata generated from last subtraction, not from the original data.
end
estTimRes = zeros(simParams.Np,simParams.N);
ResCell = {}; % a cell of timing residuals to store all the estimated residuals.
SNRarray = [];

%% MBLT
[file,Index]=rassign(estDataDir,outputfiles,NestSrc,simParams,yr);
% disp(["File needs to be skipped: ",file]);
foldername = [num2str(stage),'iMBLT_after'];
OutputDir = [estDataDir,filesep,foldername];
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
for nsrc = 1:NestSrc
    estTimRes = estTimRes + ResCell{nsrc};
end

newFile = strcat(OutputDir,filesep,inputFileName,'_after_',num2str(stage),'iMBLT',ext);

if stage == 1
    copyfile([simDataDir,filesep,inputFileName,ext],newFile);
else
    copyfile([simDataDir,filesep,inputFileName,'_after_',num2str(stage-1),'iMBLT',ext],newFile);
end
m = matfile(newFile,'Writable',true);
m.timingResiduals = timingResiduals - estTimRes;



toc
% END