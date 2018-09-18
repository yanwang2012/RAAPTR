% This is a test.m file to generate the necessary
% parameters and constants for function gensimdata
% Yi-Qian, Sep 18, 2018
%% ======== constants ===========
pc2ly=3.261563777;  % 1 pc=3.26 ly (Julian)
dy2yr=1.0/365.25;  % 1 day=365.25 yr (Julian)  ??
kilo=1.0*10^3;  % kilo 1000
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
load('parameter.mat')

save('input.mat')