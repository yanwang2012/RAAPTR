% Plot simulated noises using realistic sd from NANOGrav 12.5 yr data.

% Author: QYQ
% Date: 11/21/2022

clear;
%% load data
simDataDir = '/Users/yiqianqian/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/SimNANOGrav';
fileName = 'GWBsimDataSKASrlz1Nrlz1';
ext = '.mat';
load([simDataDir, filesep, fileName, ext]);

for n = 1:1:length(psr_data)
    yr = psr_data{n}.yr;
    name = psr_data{n}.psr_name;
    noise = noises{n};
    plot(yr, noise, '*b')
    savefig([simDataDir, filesep, 'noises', filesep, 'noise_', name]);
    saveas(gcf,[simDataDir, filesep, 'noises', filesep, 'noise_', name], 'png');
end
close all % close all open figures

% EOS