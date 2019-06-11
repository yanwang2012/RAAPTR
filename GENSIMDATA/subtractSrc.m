function []=subtractSrc(estTimRes,DataDir,inputFile,outputFile)
% Source subtraction script
% subtract the estimated timing residuals from simulation data and create 
% a new data file.

% QYQ 2019.6.9

% clear;
% timResFile = 'estTimRes.mat';
% load(timResFile);
load([DataDir,filesep,inputFile])
timingResiduals = timingResiduals - estTimRes;
disp("Input file name "+inputFile);
disp("Output file name "+outputFile);
save(outputFile, '-regexp',...
     '^(?!(estTimRes|DataDir|inputFile|outputFile)$).');% save all the variables except estTimRes
