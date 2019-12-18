function [filename,Index]=rassign(estDataDir,outputfiles,N,simParams,yr)
% [filename,Index]=rassign(estDataDir,outputfiles,N,simParams,yr)
% MBLT data pre-processing
% Reassigning the close point near the edge

% Author Yiqian Qian 2019
%% Set up
%clear;
%tic
% simParamsDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test10/searchParams';
% simParamsName = 'searchParams_Nyquist';
% inParamsList = dir([simParamsDir,filesep,simParamsName,'*.mat']);
% simDataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test10/FullBand/test';
% estDataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test10/FullBand/test/MBLT2/Input';
% outputfiles = dir([estDataDir,filesep,'*.mat']);
% inputFileName = 'GWBsimDataSKASrlz1Nrlz3';
% Npara = length(inParamsList);
% NestSrc = length(outputfiles);
% nFile = dir([estDataDir,filesep,'1_',inputFileName,'*.mat']); % count how many iterations are used.
% num_ite = length(nFile);
% % Load the simulated source parameters.
% load([simDataDir,filesep,inputFileName,'.mat']);
% estTimRes = zeros(simParams.Np,simParams.N);

%% load data
bestRealLoc = zeros(1007,N);
SNR = zeros(N,1);
for i = 1:N
    f = load([estDataDir,filesep,char(outputfiles(i).name)],'bestRealLoc');
    bestRealLoc(:,i) = f.bestRealLoc;
    [SrcParam]=ColSrcParams([estDataDir,filesep,outputfiles(i).name]);
    [SNR(i,1),~]=Amp2Snr(SrcParam,simParams,yr);
end

SrcP = zeros(N,2);

for j = 1:N
    SrcP(j,1) = bestRealLoc(3,j)/(2*pi*365*24*3600); % Frequency
    SrcP(j,2) = SNR(j,1);
end

for k = 1:N
    for m = 1:N
        if abs(SrcP(k,2)-SrcP(m,2)) < 3 && k ~= m % threshold for SNR
            if abs(SrcP(k,1)-SrcP(m,1)) < 5e-10 && k ~= m % threshold for Freq
                if SrcP(k,2) > SrcP(m,2)
                    disp(["File needs to be skipped is : ",outputfiles(m).name,'Index is: ',num2str(m)])
                    filename = outputfiles(m).name;
                    Index = m;
                else
                    disp(["File needs to be skipped is : ",outputfiles(k).name,'Index is: ',num2str(k)])
                end
            end
        end
    end
end




%toc
% EOF