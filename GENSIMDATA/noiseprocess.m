%% Noise processing
function [avgnoise] = noiseprocess(estNoiseDir,simNoiseDir)
%tic
%clear;
%estNoiseDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test10/noise/Results';
estNoiseFiles = dir([estNoiseDir,filesep,'*.mat']);
%simNoiseDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test10/noise';
simNoiseFiles = dir([simNoiseDir,filesep,'*.mat']);
numSimNoise = length(simNoiseFiles);
numEstNoise = length(estNoiseFiles);
num_ite = 10;
num_bands = 5;
estNoiseName = {};
for i = 1:numEstNoise
    estNoiseName = [estNoiseName estNoiseFiles(i).name];
end
simNoiseName = {};
for j = 1:numSimNoise
    simNoiseName = [simNoiseName simNoiseFiles(j).name];
end
estNoiseName = sort_nat(estNoiseName);
estNoiseName = reshape(estNoiseName,num_ite*numSimNoise,num_bands);
simNoiseName = sort_nat(simNoiseName);
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