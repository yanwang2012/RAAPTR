% This is a test.m file to generate the necessary
% parameters and constants for function gensimdata
% Yi-Qian, Sep 18, 2018
%% =========== parameters ============
% Ns: number of GW sources
% Np: number of Pulsars
% Nrlz: number of noise realizations H1
% Nnis: number of realization of noise only cases H0
tic;
clear;
Ns = 200;
Np = 1000; % for the specific pulsar catalog
Nrlz = 1; % number of H1 realization with different noise realizations
Nsrlz = 4; % number of H1 realization with different src realizations
Nnis = 3; % number of H0 realization
%frqRng = [8.26e-8,4.130e-7];% Angular velocity is converted to Hz by divided by 2*pi*365*24*3600
outDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/simData/diff_srlz';
mkdir(outDir);

parameters(Ns,Np,Nrlz,Nsrlz,Nnis,outDir);
% when run 'gensimdata' multiple times need to change 'Srlzs' in it, if rng
% seed is not fixed. When rng seed is fixed, it always gives the same Srlz.
gensimdata([outDir,filesep,'GWBsimDataSKAParams.mat'],...
           'survey_ska.mat',...
            outDir)%,...
%             [5.0468e-09,2e-7]); % frequency range
toc;