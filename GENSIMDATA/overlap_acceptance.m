%% Overlap-acceptance
% generate a new search parameter files and remove estimated sources in a
% frequency interval.

% QYQ 07/31/2020

%% generate new search parameter file
clear
tic

simParamsDir = '/Users/qyq/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/searchParams/2bands/superNarrow';
outputDir = '/Users/qyq/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/searchParams/2bands/overlap2';
mkdir(outputDir)
searchRange = load([simParamsDir,filesep,'searchParams1.mat']);

times =  2; % parameter to tune the gap
delta = times * (searchRange.searchParams.angular_velocity(1) - searchRange.searchParams.angular_velocity(2))/10; % overlaping frequency interval
range_low = searchRange.searchParams.angular_velocity(1) - delta;
copyfile([simParamsDir,filesep,'searchParams1.mat'],[outputDir,filesep,'searchParams2.mat']);
f = matfile([outputDir,filesep,'searchParams2.mat'],'Writable',true);
searchParams = f.searchParams;
searchParams.angular_velocity(1,1) = searchRange.FreqRange(1); % upper limit of new band
searchParams.angular_velocity(2,1) = range_low - delta; % lower limit of new band

stage = 2; % band #

searchParams.band_num = stage;
f.searchParams = searchParams;

% %% Set up
% simDataDir = '/Users/qyq/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands';
% estDataDir = '/Users/qyq/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands/Results_supNar';
% inputFileName = 'GWBsimDataSKASrlz1Nrlz3';
% ext = '.mat';
% 
% 
% % Load the simulated source parameters.
% load([simDataDir,filesep,inputFileName,ext],'simParams','yr');
% outputfiles = dir([estDataDir,filesep,num2str(stage - 1),'_',inputFileName,'*',ext]);
% NestSrc = length(outputfiles);
% outputfilenames = sort_nat({outputfiles.name});
% 
% estTimRes = zeros(simParams.Np,simParams.N);
% ResCell = {}; % a cell of timing residuals to store all the estimated residuals.
% % SNRarray = [];
% 
% for j = 1:NestSrc
%     
%     
%     path_estData = [estDataDir,filesep,char(outputfilenames(j))];
%     [srcParams]=ColSrcParams(path_estData);
%     if srcParams.omega < searchRange.searchParams.angular_velocity(1,1) && srcParams.omega > range_low
%         [~,estTimRes_tmp] = Amp2Snr(srcParams,simParams,yr);
%         ResCell = [ResCell estTimRes_tmp];
%         %     SNRarray = [SNRarray SNR];
%     else
%         continue
%     end
%     
% end
% 
% % Accumulate the timing residulas for different sources.
% num_src = length(ResCell);
% for nsrc = 1:num_src
%     estTimRes = estTimRes + ResCell{nsrc};
% end
% 
% newFile = [outputDir,filesep,num2str(stage),'_',inputFileName,ext];
% copyfile([simDataDir,filesep,inputFileName,ext],newFile);
% m = matfile(newFile,'Writable',true);
% m.timingResiduals = m.timingResiduals - estTimRes;
% toc


% END