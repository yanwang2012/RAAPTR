function []=subtractSrc(estTimRes,simDir,infilename,outfilename)
% Source subtract function
% []=subtractSrc(estTimRes,simDir,filename,outFilename)
% subtract the estimated timing residuals from simulation data and create 
% a new data file with the same name

% QYQ 2019.6.9

load([simDir,filesep,infilename]);
timingResiduals = timingResiduals - estTimRes;
% disp("File name is "+infilename);
% disp("outFile is "+outfilename);
save(outfilename);