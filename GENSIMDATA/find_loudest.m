% script to load simulated NANOGrav data and find out the loudest injected
% sources.

% Yi-qian Qian 2023-04-11

clear;

dataDir = '/Users/yiqianqian/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/Lomb-Scargle';
files = dir([dataDir, filesep,'GWBsimDataSKASrlz*Nrlz1.mat']);
filenames = sort_nat({files.name});
N = length(filenames);
snr = zeros(N);

for file = 1:N
    contents = load([dataDir, filesep, filenames{file}]);
    snr(file) = contents.snr_chr;
end

%% Plot
srlz = 1:1:N;
scatter(srlz,snr)
xlabel('Source Realization')
ylabel('SNR')
title('SNR dist')