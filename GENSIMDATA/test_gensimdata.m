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
Ns = 100;
Np = 1000; % for the specific pulsar catalog
Nrlz = 5;
Nnis = 3;
parameters(Ns,Np,Nrlz,Nnis)

mkdir('TESTDATA');

gensimdata('GWBsimDataSKAParams.mat',...
           'survey_ska.mat',...
            'TESTDATA')
toc;