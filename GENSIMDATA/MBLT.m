% MBLT Implementation for Multisource
% QYQ 23/11/2019
clear;
tic
%% Set up
simParamsDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test10/searchParams';
simParamsName = 'searchParams_Nyquist';
inParamsList = dir([simParamsDir,filesep,simParamsName,'*.mat']);
simDataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test10/FullBand/test';
estDataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test10/FullBand/test/MBLT2/Input';
outputfiles = dir([estDataDir,filesep,'*.mat']);
inputFileName = 'GWBsimDataSKASrlz1Nrlz3';
Npara = length(inParamsList);
NestSrc = length(outputfiles);
nFile = dir([estDataDir,filesep,'1_',inputFileName,'*.mat']); % count how many iterations are used.
num_ite = length(nFile);
% Load the simulated source parameters.
load([simDataDir,filesep,inputFileName,'.mat']);
estTimRes = zeros(simParams.Np,simParams.N);

%% MBLT
[file,Index]=rassign(estDataDir,outputfiles,NestSrc,simParams,yr);
% disp(["File needs to be skipped: ",file]);
OutputDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test10/FullBand/test/MBLT2';
for i = 1:Npara
    for j = 1:NestSrc
        if j<= (i-1)*num_ite || j > i*num_ite
            if j == Index
                continue
            else
                disp("j is:"+j);
                path_estData = [estDataDir,filesep,outputfiles(j).name];
                [srcParams]=ColSrcParams(path_estData);
                [~,estTimRes_tmp] = Amp2Snr(srcParams,simParams,yr);
                estTimRes = estTimRes + estTimRes_tmp;
            end
        end
    end
    newFile = strcat(OutputDir,filesep,inputFileName,'band ',num2str(i),'.mat');
    copyfile([simDataDir,filesep,inputFileName,'.mat'],newFile);
    m = matfile(newFile,'Writable',true);
    m.timingResiduals = timingResiduals - estTimRes;
    estTimRes = zeros(simParams.Np,simParams.N);
end


toc
% END