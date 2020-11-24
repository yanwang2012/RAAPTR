% MBLT Implementation for Multisource
% QYQ 23/11/2019
clear;
tic
%% Set up
simParamsDir = '/work/05884/qyqstc/lonestar/MultiPSO/Task8/searchParams/2bands/superNarrow';
simParamsName = 'searchParams';
inParamsList = dir([simParamsDir,filesep,simParamsName,'*.mat']);
simDataDir = '/work/05884/qyqstc/lonestar/MultiPSO/Task8/BANDEDGE/2bands/superNarrow/simData';
estDataDir = '/work/05884/qyqstc/lonestar/MultiPSO/Task8/BANDEDGE/2bands/superNarrow/Union2_xMBLT/results';
inputFileName = 'GWBsimDataSKASrlz1Nrlz3';

Npara = length(inParamsList);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DON'T FORGET TO CHECK THE NAME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nFile = dir([estDataDir,filesep,'1_',inputFileName,'*.mat']); % count how many iterations are used.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_ite = length(nFile);
% Load the simulated source parameters.
load([simDataDir,filesep,inputFileName,'.mat']);
estTimRes = zeros(simParams.Np,simParams.N);

%% MBLT
outputfiles = dir([estDataDir,filesep,'*',inputFileName,'*.mat']);
NestSrc = length(outputfiles);
Nband1 = NestSrc/2;
outputfilenames = sort_nat({outputfiles.name});
[file,Index]=rassign(estDataDir,outputfilenames,NestSrc,Nband1,simParams,yr);
% disp(["File needs to be skipped: ",file]);
Filename = 'Union2_xMBLT2';
OutputDir = [estDataDir,filesep,Filename];
mkdir(OutputDir);
for i = 1:Npara
    for j = 1:NestSrc
        if j <= (i-1)*num_ite || j > i*num_ite
            if ismember(j,Index) == 1
                continue
            else
%                 disp("j is:"+j);
                path_estData = [estDataDir,filesep,char(outputfilenames(j))];
                disp(['File loaded: ',char(outputfilenames(j))]);
                [srcParams]=ColSrcParams(path_estData);
                [~,estTimRes_tmp] = Amp2Snr(srcParams,simParams,yr);
                estTimRes = estTimRes + estTimRes_tmp;
            end
        end
    end
    newFile = strcat(OutputDir,filesep,num2str(i),'_',inputFileName,'.mat');
    copyfile([simDataDir,filesep,inputFileName,'.mat'],newFile);
    m = matfile(newFile,'Writable',true);
    m.timingResiduals = timingResiduals - estTimRes;
    estTimRes = zeros(simParams.Np,simParams.N);
end


toc
% END
