% MBLT in ONE band
% Updated to multiple version
% 05/31/2021 QYQ

clear;
tic
%% Set up
simDataDir = '/work2/05884/qyqstc/stampede2/MultiPSO/Final/realizations/2bands/results_diff_opt_xMBLT2/1_iMBLT/results/1iMBLT_after/results/2_iMBLT/results/2iMBLT_after/results/3_iMBLT/results/3iMBLT_after/results/4_iMBLT/results/4iMBLT_after/results/5_iMBLT/results/5iMBLT_after'; % for 1st stage
% searchParamsDir = '/work2/05884/qyqstc/stampede2/MultiPSO/Final/realizations/2bands/searchParams/Band_opt';
searchParamsDir = '/work2/05884/qyqstc/stampede2/MultiPSO/Yuyang/searchParams/Band_opt';
iMBLTDataDir = '/work2/05884/qyqstc/stampede2/MultiPSO/Yuyang/1000Psr/results_xMBLT2/1_iMBLT/results/1iMBLT_after/results/2_iMBLT/results/2iMBLT_after';
estDataDir = '/work2/05884/qyqstc/stampede2/MultiPSO/Yuyang/1000Psr/results_xMBLT2/1_iMBLT/results/1iMBLT_after/results/2_iMBLT/results/2iMBLT_after/results';
% srlzDir = '/work2/05884/qyqstc/stampede2/MultiPSO/Final/realizations/2bands/simData/Band_opt_diff';
srlzDir = '/work2/05884/qyqstc/stampede2/MultiPSO/Yuyang/1000Psr';
inputFileName = 'GWBsimDataSKASrlz*Nrlz1';
NsimFiles = dir([srlzDir,filesep,inputFileName,'.mat']);
simFileNames = sort_nat({NsimFiles.name});
Nrlzs = length(NsimFiles);
spName = 'searchParams_Nyquist*';
stage = 3; % iMBLT stages

for rlz = 1:Nrlzs
    [~,simFileBaseName,~] = fileparts(simFileNames{rlz});
    NspFiles = dir([searchParamsDir,filesep,simFileBaseName,filesep,spName]);
    bandNum = length(NspFiles); % number of band
    
    for band = 1:bandNum
        if bandNum > 1
            outputfiles = dir([estDataDir,filesep,num2str(band),'_',simFileBaseName,'*.mat']);
        else
            outputfiles = dir([estDataDir,filesep,'*',simFileBaseName,'*.mat']);
        end
        % Npara = length(inParamsList);
        NestSrc = length(outputfiles);
        
        % Load the simulated source parameters.
        if stage == 1
            load([simDataDir,filesep,num2str(band),'_',simFileBaseName,'.mat']); % for first iMBLT stage.
        else
            load([iMBLTDataDir,filesep,num2str(band),'_',simFileBaseName,'_after_',num2str(stage-1),'iMBLT','.mat'])
        end
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
            path_estData = [estDataDir,filesep,outputfilenames{j}];
            [srcParams]=ColSrcParams(path_estData,simParams.Np);
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
        newFile = strcat(OutputDir,filesep,num2str(band),'_',simFileBaseName,'_',num2str(stage),'iMBLT','.mat');
        if stage == 1
            copyfile([simDataDir,filesep,num2str(band),'_',simFileBaseName,'.mat'],newFile); % for first iMBLT stage
        else
            copyfile([iMBLTDataDir,filesep,num2str(band),'_',simFileBaseName,'_after_',num2str(stage-1),'iMBLT','.mat'],newFile); % can't use * to match files
        end
        m = matfile(newFile,'Writable',true); % .mat file needs to be V7.3, if can't load, use command save('filename','-v7.3') convert
        m.timingResiduals = timingResiduals - estTimRes;
        
    end
end


toc
% END
