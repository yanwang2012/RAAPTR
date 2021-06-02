%% Remove estimated sources
% Given a set of output files, remove all the estimated sources from
% simulation data.
% QYQ 03/11/2019

clear;
tic
%% Set up
simDataDir = '/work/05884/qyqstc/lonestar/MultiPSO/Task8/BANDEDGE/2bands/superNarrow/Union2_xMBLT/results/Union2_xMBLT2/results/1_iMBLT/results/1iMBLT_after/results/2_iMBLT/results/2iMBLT_after/results/3_iMBLT/results/3iMBLT_after/results/4_iMBLT/results/4iMBLT_after/results/5_iMBLT/results/5iMBLT_after/results/6_iMBLT/results/6iMBLT_after/results/7_iMBLT/results/7iMBLT_after/results/8_iMBLT/results/8iMBLT_after/results/9_iMBLT/results/9iMBLT_after/results/10_iMBLT/results/10iMBLT_after/results/11_iMBLT/results/11iMBLT_after/results/12_iMBLT/results/12iMBLT_after/results/13_iMBLT/results/13iMBLT_after/results/14_iMBLT/results/14iMBLT_after/results/15_iMBLT/results/15iMBLT_after/results/16_iMBLT/results/16iMBLT_after/results/17_iMBLT/results/17iMBLT_after/results/18_iMBLT/results/18iMBLT_after';
searchParamsDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/searchParams/Band_opt';
estDataDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/results_diff_opt_xMBLT2/1_iMBLT/results';
srlzDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/simData/Band_opt_diff';
inputFileName = 'GWBsimDataSKASrlz*Nrlz1';
ext = '.mat';

NsimFiles = dir([srlzDir,filesep,inputFileName,ext]);
simFileNames = sort_nat({NsimFiles.name});
exp = 'GWBsimDataSKASrlz\dNrlz1.*';
simFileNames = regexp(simFileNames,exp,'match');
simFileNames = simFileNames(~cellfun(@isempty,simFileNames));
Nrlzs = length(simFileNames);
stage = 1; % stage number of iMBLT
spName = 'searchParams_Nyquist*';

for rlz = 1:Nrlzs
    [~,simFileBaseName,~] = fileparts(simFileNames{rlz});
    NspFiles = dir([searchParamsDir,filesep,simFileBaseName,filesep,spName]);
    bandNum = length(NspFiles); % number of band
    
    % Load the simulated source parameters.
    for i = 1:bandNum
        outputfiles = dir([estDataDir,filesep,num2str(i),'_',simFileBaseName,'*',ext]);
        NestSrc = length(outputfiles);
        if stage == 1
            simDataDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/simData/Band_opt_diff/Band_opt_diff_xMBLT2';
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