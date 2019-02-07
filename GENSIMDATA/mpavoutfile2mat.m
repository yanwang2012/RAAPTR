function [] = mpavoutfile2mat(inFile,varargin)
%Convert output file in HDF5 format to MAT file with structures restored
%MPAVOUTFILE2MAT(F)
%Translates the .hdf5 output file produced by perfeval_omp.c
%to .mat format. Structures in the .mat file are
%restored. The default name of the output file is the same as the input file
%but with the extension .hdf5 replaced by .mat. If F contains path
%information, the same path is used for the output file.
%
%MPAVINFILE2HDF5(F,O)
%Overrides the default output file name with O. The path information for O
%can be different from that of F.

%Soumya D. Mohanty, Mar 2017

%Read the contents of the input file
inFileInfo = loadhdf5(inFile);

%Restore the structure fields
inFileInfo.id = struct('snr_id',inFileInfo.snr_id,...
                       'loc_id',inFileInfo.loc_id,...
                       'omg_id',inFileInfo.omg_id,...
                       'rlz_id',inFileInfo.rlz_id);
inFileInfo = rmfield(inFileInfo,{'snr_id','loc_id','omg_id','rlz_id'});

%Default output file name
[inFilePath,outFile,~]=fileparts(inFile);
outFile = [inFilePath,filesep,outFile];
if nargin > 1
    outFile = varargin{1};
end

%save output 
save(outFile,'-struct','inFileInfo');

    

