% This is a test.m file to generate the necessary
% parameters and constants for function gensimdata
% Yi-Qian, Sep 18, 2018
%% =========== parameters ============
% NumGwsources: number of GW sources
% NumPulsar: number of Pulsars
% NumNoiseReali: number of noise realizations H1
% NumRealiNoise: number of realization of noise only cases H0
NumGWsources = 100;
NumPulsar = 1000;
NumNoiseReali = 5;
NumRealiNoise = 1;
parameters(NumGWsources,NumPulsar,NumNoiseReali,NumRealiNoise)

mkdir('TESTDATA');
gensimdata('parameter.mat','survey_ska.mat','TESTDATA')