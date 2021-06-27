%% Remove estimated sources
% Given a set of output files, remove all the estimated sources from
% simulation data.
% QYQ 03/11/2019

clear;
tic
%% Set up
simDataDir = '/work2/05884/qyqstc/stampede2/MultiPSO/Yuyang/1000Psr/results_xMBLT2/1_iMBLT/results/1iMBLT_after/results/2_iMBLT/results/2iMBLT_after';
% searchParamsDir = '/work2/05884/qyqstc/stampede2/MultiPSO/Final/realizations/2bands/searchParams/Band_opt';
searchParamsDir = '/work2/05884/qyqstc/stampede2/MultiPSO/Yuyang/searchParams/Band_opt'
estDataDir = '/work2/05884/qyqstc/stampede2/MultiPSO/Yuyang/1000Psr/results_xMBLT2/1_iMBLT/results/1iMBLT_after/results/2_iMBLT/results/2iMBLT_after/results/3_iMBLT/results';
% srlzDir = '/work2/05884/qyqstc/stampede2/MultiPSO/Final/realizations/2bands/simData/Band_opt_diff';
srlzDir = '/work2/05884/qyqstc/stampede2/MultiPSO/Yuyang/1000Psr'
inputFileName = 'GWBsimDataSKASrlz*Nrlz1';
ext = '.mat';

NsimFiles = dir([srlzDir,filesep,inputFileName,ext]);
simFileNames = sort_nat({NsimFiles.name});
exp = 'GWBsimDataSKASrlz\dNrlz1.*';
simFileNames = regexp(simFileNames,exp,'match');
simFileNames = simFileNames(~cellfun(@isempty,simFileNames));
Nrlzs = length(simFileNames);
stage = 3; % stage number of iMBLT
spName = 'searchParams_Nyquist*';

for rlz = 1:Nrlzs
    [~,simFileBaseName,~] = fileparts(simFileNames{rlz}{1});
    NspFiles = dir([searchParamsDir,filesep,simFileBaseName,filesep,spName]);
    bandNum = length(NspFiles); % number of band
    
    % Load the simulated source parameters.
    for i = 1:bandNum
        outputfiles = dir([estDataDir,filesep,num2str(i),'_',simFileBaseName,'*',ext]);
        NestSrc = length(outputfiles);
        if stage == 1
            simDataDir = '/work2/05884/qyqstc/stampede2/MultiPSO/Yuyang/1000Psr/xMBLT2';
            load([simDataDir,filesep,num2str(i),'_',simFileBaseName,ext])
        else
            load([simDataDir,filesep,num2str(i),'_',simFileBaseName,'_after_',num2str(stage-1),'iMBLT',ext]);% subtraction in iMBLT needs to base on the new inputdata generated from last subtraction, not from the original data.
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
            [srcParams]=ColSrcParams(path_estData, simParams.Np);
            [SNR,estTimRes_tmp] = Amp2Snr(srcParams,simParams,yr);
            ResCell = [ResCell estTimRes_tmp];
            SNRarray = [SNRarray SNR];
            %         end
        end
        
        % Accumulate the timing residulas for different sources.
        for nsrc = 1:NestSrc
            estTimRes = estTimRes + ResCell{nsrc};
        end
        
        newFile = strcat(OutputDir,filesep,num2str(i),'_',simFileBaseName,'_after_',num2str(stage),'iMBLT',ext);
        
        if stage == 1
            copyfile([simDataDir,filesep,num2str(i),'_',simFileBaseName,ext],newFile);
        else
            copyfile([simDataDir,filesep,num2str(i),'_',simFileBaseName,'_after_',num2str(stage-1),'iMBLT',ext],newFile);
        end
        m = matfile(newFile,'Writable',true);
        m.timingResiduals = timingResiduals - estTimRes;
    end
end

toc
% END
