% MBLT in ONE band
% QYQ 02/01/2019
clear;
tic
%% Set up
simDataDir = '/work/05884/qyqstc/lonestar/MultiPSO/Task8/BANDEDGE/2bands/SupNar_xMBLT_iMBLT20/GWBsimDataSKASrlz1Nrlz3_xMBLT2';
iMBLTDataDir = '/work/05884/qyqstc/lonestar/MultiPSO/Task8/BANDEDGE/2bands/SupNar_xMBLT_iMBLT20/GWBsimDataSKASrlz1Nrlz3_xMBLT2';
estDataDir = '/work/05884/qyqstc/lonestar/MultiPSO/Task8/BANDEDGE/2bands/SupNar_xMBLT_iMBLT20/GWBsimDataSKASrlz1Nrlz3_xMBLT2/results';
inputFileName = 'GWBsimDataSKASrlz1Nrlz3';
bandNum = 2; % number of band
stage = 1; % iMBLT stages
for band = 1:bandNum
    if bandNum > 1
        outputfiles = dir([estDataDir,filesep,num2str(band),'_',inputFileName,'*.mat']);
    else
        outputfiles = dir([estDataDir,filesep,'*',inputFileName,'*.mat']);
    end
    % Npara = length(inParamsList);
    NestSrc = length(outputfiles);
    
    % Load the simulated source parameters.
    load([simDataDir,filesep,num2str(band),'_',inputFileName,'.mat']);
    % load([iMBLTDataDir,filesep,inputFileName,'_after_',num2str(stage-1),'iMBLT','.mat'])
    estTimRes = zeros(simParams.Np,simParams.N);
    ResCell = {}; % a cell of timing residuals to store all the estimated residuals.
    SNRarray = [];
    
    %% iMBLT
    [file,Index]=rassign(estDataDir,outputfiles,NestSrc,simParams,yr);
    % disp(["File needs to be skipped: ",file]);
    foldername = [num2str(stage),'_iMBLT'];
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
    if NestSrc >= 6
        for nsrc = 2:6%NestSrc remove first 5 src.
            estTimRes = estTimRes + ResCell{nsrc};
        end
    elseif 1 < NestSrc && NestSrc < 6
        for nsrc = 2:NestSrc % when output files less than 5, remove all the rest sources lower than the target.
            estTimRes = estTimRes + ResCell{nsrc};
        end
    else
        disp("No need to subtract more. Go estimation!")
    end
    
    % remove sources from input data
    newFile = strcat(OutputDir,filesep,num2str(band),'_',inputFileName,'_',num2str(stage),'iMBLT','.mat');
    copyfile([simDataDir,filesep,num2str(band),'_',inputFileName,'.mat'],newFile);
    % copyfile([iMBLTDataDir,filesep,inputFileName,'_after_',num2str(stage-1),'iMBLT','.mat'],newFile); % can't use * to match files
    
    m = matfile(newFile,'Writable',true); % .mat file needs to be V7.3, if can't load, use command save('filename','-v7.3') convert
    m.timingResiduals = timingResiduals - estTimRes;
    
end

toc
% END
