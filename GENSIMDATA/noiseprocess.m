%% Noise processing
function [avgnoise] = noiseprocess(estNoiseDir,simNoiseDir,num_ite,num_bands)
% [avgnoise] = noiseprocess(estNoiseDir,simNoiseDir,num_ite,num_bands)
% num_ite: number of iterations used to estimate noise-only data.
% num_bands : number of bands used for noise-only data.
% num_ite and num_bands for noise should match the one used for H1 data.

%% Debug
% estNoiseDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test10/noise/Results20/';
% simNoiseDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test10/noise';
% num_ite = 20;
% num_bands = 5;

%% Main
estNoiseFiles = dir([estNoiseDir,filesep,'*.mat']);
simNoiseFiles = dir([simNoiseDir,filesep,'*.mat']);

numSimNoise = length(simNoiseFiles);

estNoiseName = sort_nat({estNoiseFiles.name});
estNoiseName = reshape(estNoiseName,num_ite*numSimNoise,num_bands);
simNoiseName = sort_nat({simNoiseFiles.name});
noise_tmp = 0;
avgnoise_tmp = zeros(numSimNoise,num_ite);
for m = 1:numSimNoise
    noise = load([simNoiseDir,filesep,char(simNoiseName(m))]);
    for n = 1:num_ite
        for k = 1:num_bands
            pth_estNoise_file = [estNoiseDir,filesep,char(estNoiseName((m-1)*10+n,k))];
            noiseParams = ColSrcParams(pth_estNoise_file);
            [noiseSNR,~] = Amp2Snr(noiseParams,noise.simParams,noise.yr);
            noise_tmp = noise_tmp + noiseSNR;
        end
        avgnoise_tmp(m,n) = noise_tmp/num_bands; % avarage over bands
        noise_tmp = 0;
    end
end

avgnoise = sum(avgnoise_tmp,1)/numSimNoise;

%END