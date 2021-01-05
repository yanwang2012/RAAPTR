% MBLT in ONE band
% QYQ 02/01/2019
clear;
tic
%% Set up
simDataDir = '/work/05884/qyqstc/lonestar/MultiPSO/Task8/BANDEDGE/2bands/superNarrow/simData/Union2_xMBLT_iMBLT';
iMBLTDataDir = '/work/05884/qyqstc/lonestar/MultiPSO/Task8/BANDEDGE/2bands/superNarrow/Union2_xMBLT/results/Union2_xMBLT2/results/1_iMBLT/results/1iMBLT_after/results/2_iMBLT/results/2iMBLT_after/results/3_iMBLT/results/3iMBLT_after/results/4_iMBLT/results/4iMBLT_after/results/5_iMBLT/results/5iMBLT_after/results/6_iMBLT/results/6iMBLT_after/results/7_iMBLT/results/7iMBLT_after/results/8_iMBLT/results/8iMBLT_after/results/9_iMBLT/results/9iMBLT_after/results/10_iMBLT/results/10iMBLT_after/results/11_iMBLT/results/11iMBLT_after/results/12_iMBLT/results/12iMBLT_after/results/13_iMBLT/results/13iMBLT_after/results/14_iMBLT/results/14iMBLT_after/results/15_iMBLT/results/15iMBLT_after/results/16_iMBLT/results/16iMBLT_after/results/17_iMBLT/results/17iMBLT_after/results/18_iMBLT/results/18iMBLT_after';
estDataDir = '/work/05884/qyqstc/lonestar/MultiPSO/Task8/BANDEDGE/2bands/superNarrow/Union2_xMBLT/results/Union2_xMBLT2/results/1_iMBLT/results/1iMBLT_after/results/2_iMBLT/results/2iMBLT_after/results/3_iMBLT/results/3iMBLT_after/results/4_iMBLT/results/4iMBLT_after/results/5_iMBLT/results/5iMBLT_after/results/6_iMBLT/results/6iMBLT_after/results/7_iMBLT/results/7iMBLT_after/results/8_iMBLT/results/8iMBLT_after/results/9_iMBLT/results/9iMBLT_after/results/10_iMBLT/results/10iMBLT_after/results/11_iMBLT/results/11iMBLT_after/results/12_iMBLT/results/12iMBLT_after/results/13_iMBLT/results/13iMBLT_after/results/14_iMBLT/results/14iMBLT_after/results/15_iMBLT/results/15iMBLT_after/results/16_iMBLT/results/16iMBLT_after/results/17_iMBLT/results/17iMBLT_after/results/18_iMBLT/results/18iMBLT_after/results';
inputFileName = 'GWBsimDataSKASrlz1Nrlz3';
bandNum = 2; % number of band
stage = 19; % iMBLT stages
for band = 1:bandNum
    if bandNum > 1
        outputfiles = dir([estDataDir,filesep,num2str(band),'_',inputFileName,'*.mat']);
    else
        outputfiles = dir([estDataDir,filesep,'*',inputFileName,'*.mat']);
    end
    % Npara = length(inParamsList);
    NestSrc = length(outputfiles);
    
    % Load the simulated source parameters.
    % load([simDataDir,filesep,num2str(band),'_',inputFileName,'.mat']); % for first iMBLT stage.
    load([iMBLTDataDir,filesep,num2str(band),'_',inputFileName,'_after_',num2str(stage-1),'iMBLT','.mat'])
    estTimRes = zeros(simParams.Np,simParams.N);
    ResCell = {}; % a cell of timing residuals to store all the estimated residuals.
    SNRarray = [];
    
    %% iMBLT
    % [file,Index]=rassign(estDataDir,outputfiles,NestSrc,simParams,yr);
    % disp(["File needs to be skipped: ",file]);
    foldername = [num2str(stage),'_iMBLT'];
    OutputDir = [estDataDir,filesep,foldername];
    outputfilenames = sort_nat({outputfiles.name});
    
    mkdir(OutputDir);
    
    for j = 1:NestSrc
        
        %         if ismember(j,Index) == 1
        %             continue
        %         else
        %                 disp("j is:"+j);
        path_estData = [estDataDir,filesep,char(outputfilenames(j))];
        [srcParams]=ColSrcParams(path_estData);
        [SNR,estTimRes_tmp] = Amp2Snr(srcParams,simParams,yr);
        ResCell = [ResCell estTimRes_tmp];
        SNRarray = [SNRarray SNR];
        %         end
        
    end
    
    % Accumulate the timing residulas for different sources.
%     nskp = length(Index); % number of skipped files
    
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
    % copyfile([simDataDir,filesep,num2str(band),'_',inputFileName,'.mat'],newFile); % for first iMBLT stage
    copyfile([iMBLTDataDir,filesep,num2str(band),'_',inputFileName,'_after_',num2str(stage-1),'iMBLT','.mat'],newFile); % can't use * to match files
    
    m = matfile(newFile,'Writable',true); % .mat file needs to be V7.3, if can't load, use command save('filename','-v7.3') convert
    m.timingResiduals = timingResiduals - estTimRes;
    
end

toc
% END
