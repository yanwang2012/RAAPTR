%% Remove estimated sources
% Given a set of output files, remove all the estimated sources from
% simulation data.
% QYQ 03/11/2019

clear;
tic
%% Set up
simDataDir = '/work/05884/qyqstc/lonestar/MultiPSO/Task8/BANDEDGE/2bands/superNarrow/Union2_xMBLT/results/Union2_xMBLT2';
estDataDir = '/work/05884/qyqstc/lonestar/MultiPSO/Task8/BANDEDGE/2bands/superNarrow/Union2_xMBLT/results/Union2_xMBLT2/results/1_iMBLT/results';
inputFileName = 'GWBsimDataSKASrlz1Nrlz3';
ext = '.mat';
stage = 1; % stage number of iMBLT
bandNum = 2;

% Load the simulated source parameters.
for i = 1:bandNum
    outputfiles = dir([estDataDir,filesep,num2str(i),'_',inputFileName,'*',ext]);
    NestSrc = length(outputfiles);
    if stage == 1
        simDataDir = '/work/05884/qyqstc/lonestar/MultiPSO/Task8/BANDEDGE/2bands/superNarrow/Union2_xMBLT/results/Union2_xMBLT2';
        load([simDataDir,filesep,num2str(i),'_',inputFileName,ext])
    else
        load([simDataDir,filesep,num2str(i),'_',inputFileName,'_after_',num2str(stage-1),'iMBLT',ext]);% subtraction in iMBLT needs to base on the new inputdata generated from last subtraction, not from the original data.
    end
    estTimRes = zeros(simParams.Np,simParams.N);
    ResCell = {}; % a cell of timing residuals to store all the estimated residuals.
    SNRarray = [];
    
    %% MBLT
    % [file,Index]=rassign(estDataDir,outputfiles,NestSrc,simParams,yr);
    % disp(["File needs to be skipped: ",file]);
    foldername = [num2str(stage),'iMBLT_after'];
    OutputDir = [estDataDir,filesep,foldername];
    outputfilenames = sort_nat({outputfiles.name});
    
    mkdir(OutputDir);
    
    for j = 1:NestSrc
        
        %         if j == Index
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
    for nsrc = 1:NestSrc
        estTimRes = estTimRes + ResCell{nsrc};
    end
    
    newFile = strcat(OutputDir,filesep,num2str(i),'_',inputFileName,'_after_',num2str(stage),'iMBLT',ext);
    
    if stage == 1
        copyfile([simDataDir,filesep,num2str(i),'_',inputFileName,ext],newFile);
    else
        copyfile([simDataDir,filesep,num2str(i),'_',inputFileName,'_after_',num2str(stage-1),'iMBLT',ext],newFile);
    end
    m = matfile(newFile,'Writable',true);
    m.timingResiduals = timingResiduals - estTimRes;
end


toc
% END
