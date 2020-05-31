% load the entire xmaxmin from a json file then split the frequency into
% multiple bands with random disturbance.

% Author QYQ
%%
clear;
tic
filename = 'Nyquist.json';
searchParams = jsondecode(fileread(filename));
uplim = searchParams.angular_velocity(1);
lowlim = searchParams.angular_velocity(2);
outDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/searchParams/bandTrials';
%intval = 1e-7; % disturbance interval
trials = 100; % number of trials
NumBands = 5; %10; % total number of bands
bandwidth = searchParams.angular_velocity(1)/NumBands;% max frequency divided into 10 bands
FreqRange = searchParams.angular_velocity;
for k = 1:trials
    folder = [outDir,filesep,'trial_',num2str(k)];
    mkdir(folder);
    for i = 1:NumBands
        searchParams.angular_velocity(1) = bandwidth*i + randn;
        if i == 1
            searchParams.angular_velocity(2) = lowlim;
        elseif i == NumBands
            searchParams.angular_velocity(1) = uplim;
            f = matfile([folder,filesep,'searchParams_Nyquist',num2str(i-1),'.mat']);
            temp = f.searchParams;
            searchParams.angular_velocity(2) = temp.angular_velocity(1);
        else
            f = matfile([folder,filesep,'searchParams_Nyquist',num2str(i-1),'.mat']);
            temp = f.searchParams;
            searchParams.angular_velocity(2) = temp.angular_velocity(1);
        end
        searchParams.band_num = i;
        save([folder,filesep,'searchParams_Nyquist',num2str(i),'.mat'],'searchParams','NumBands','FreqRange');
    end
end
toc
%END