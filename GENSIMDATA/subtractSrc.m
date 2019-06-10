function []=subtractSrc(estTimRes,simDir,filename,outFile)
% Source subtract function
% []=subtractSrc(estTimRes,simDir,filename,outFilename)
% subtract the estimated timing residuals from simulation data and create 
% a new data file with the same name

% QYQ 2019.6.9

load([simDir,filesep,filename]);
timingResiduals = timingResiduals - estTimRes;
save(outFile);