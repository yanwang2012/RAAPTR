% Script to plot frequency gaps between estimated sources in order to place
% the band edge more properly.

% Author: QYQ
% 1/7/2021
%% Data Dir
clear

UnionDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/diff_srlz_xMBLT2_results/Union';
simDataDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/simData/diff_srlz';
simFileName = 'GWBsimDataSKASrlz*Nrlz1';
outDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/searchParams';

simFiles = dir([simDataDir,filesep,simFileName,'*.mat']);
simFilenames = sort_nat({simFiles.name});
Nrlzs = length(simFiles);

for rlz = 1:Nrlzs
    [~,simFile_noExt,~] = fileparts(simFilenames{rlz});
    estSrcFile = dir([UnionDir,filesep,simFile_noExt,filesep,'*.mat']); % for multiple rlzs
%     estSrcFile = dir([UnionDir,filesep,'*.mat']); % for single rlz.
    estSrcFilename = sort_nat({estSrcFile.name});
    Nfiles = length(estSrcFile);
    
    %% load estimated sources
    Np = 100; %1000; % # of pulsars used in simulation
    estSrc = {};
    estFreq = zeros(1,Nfiles);
    
    for i = 1:Nfiles
        path_to_estSrc = [UnionDir,filesep,simFile_noExt,filesep,estSrcFilename{i}]; % multiple
%         path_to_estSrc = [UnionDir,filesep,estSrcFilename{i}]; % single
        estSrc{i} = ColSrcParams(path_to_estSrc,Np);
        estFreq(i) = estSrc{i}.omega/(2*pi*365*24*3600); % convert rad/year -> Hz
    end
    
    estFreq = sort(estFreq); % sort frequencies in ascending order
    
    % calculate joint freq difference
    difreq = zeros(1,Nfiles);
    for j = 1:Nfiles - 1
        difreq(j) = abs(estFreq(j+1) - estFreq(j))/estFreq(j); % relative frequency gap
    end
    
%     threshold = 0.4;
    [~,idx_tmp] = sort(difreq,'descend');
    idx = idx_tmp(2:3); % first 2 big gap, exclude the first 1;
%     idx = find(difreq > 0.3 & difreq < max_dif); % exclude the max difference
    bandedge = (1 + 1/2 * difreq(idx)) .* estFreq(idx);
    bandedge_av = bandedge * 365*24*3600*2*pi; % convert Hz to rad/year
    
    %% Generate search file
    filename = 'Nyquist.json';
    searchParams = jsondecode(fileread(filename));
    NumBands = length(idx) + 1;%5; %10; % n edge splits whole band into n+1 bands
    % bandwidth = bandedge_av;
    FreqRange = searchParams.angular_velocity;
    dest = [outDir,filesep,'Band_opt',filesep,simFile_noExt];
    mkdir(dest);
    
    for i = 1:NumBands % i band edges can split whole band into i+1 segments
        if i == NumBands
            searchParams.angular_velocity(1) = FreqRange(1);
            save([dest,filesep,'searchParams_Nyquist',num2str(i),'.mat'],'searchParams','NumBands','FreqRange');
        else
            searchParams.angular_velocity(1) = bandedge_av(i); % update upper limit
            searchParams.band_num = i;
            save([dest,filesep,'searchParams_Nyquist',num2str(i),'.mat'],'searchParams','NumBands','FreqRange');
            tmp = searchParams.angular_velocity(1); % save upper limit for prev. band
            searchParams.angular_velocity(2) = tmp; % update lower limit
        end
    end
    
    
    %% plot
    saveDir = [dest,filesep,'fig'];
    mkdir(saveDir);
    figure
    plot(estFreq,difreq,'ro')
    xlabel('Frequency [Hz]')
    ylabel('Relative Freq. Diff')
    title('Frequency Gap')
    saveas(gcf,[saveDir,filesep,'FreqGap.png'])
    savefig([saveDir,filesep,'FreqGap'])
end

%END