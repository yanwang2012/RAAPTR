function [] = mpavinfile2hdf5(inFile,varargin)
%Convert input data file to HDF5 format file with structures exploded
%MPAVINFILE2HDF5(I)
%Translates an input data file I used in Monte Carlo simulations for
%Max/AvPhase from .mat to .hdf5 format. Structures in the .mat file are
%exploded and the fields are stored out as separate datasets in the .hdf5
%file. The default name of the output file is the same as the input file
%but with the extension .mat replaced by .hdf5. If I contains path
%information, the same path is used for the output file.
%
%MPAVINFILE2HDF5(I,F)
%Overrides the default output file name with F. The path information for F
%can be different from that of I.

%Soumya D. Mohanty, Mar 2017

%Read the contents of the input file
inFileInfo = load(inFile);

%Explode the structure fields
inFileInfo.Np = inFileInfo.simParams.Np;
inFileInfo.N = inFileInfo.simParams.N;
inFileInfo.sd = inFileInfo.simParams.sd;
inFileInfo.alphaP = inFileInfo.simParams.alphaP;
inFileInfo.deltaP = inFileInfo.simParams.deltaP;
inFileInfo.kp = inFileInfo.simParams.kp;

inFileInfo.snr_id = inFileInfo.id.snr_id;
inFileInfo.loc_id = inFileInfo.id.loc_id;
inFileInfo.omg_id = inFileInfo.id.omg_id;
inFileInfo.rlz_id = inFileInfo.id.rlz_id;

inFileInfo = rmfield(inFileInfo,'simParams');
inFileInfo = rmfield(inFileInfo,'id');

%Remove the place holder cell field
inFileInfo = rmfield(inFileInfo,'pname');

%Generate output file name
[inFilePath,inFileName,~] = fileparts(inFile);
if nargin < 2
    %default file name
    outFile = [inFilePath,filesep,inFileName,'.hdf5'];
else
    outFile = varargin{1};
end

%Temporaty .mat file
outFileTmp = [inFilePath,filesep,inFileName,'_tmp.mat'];
save(outFileTmp,'-struct','inFileInfo');

%Convert temporary file to HDF5
matfile2hdf5(outFileTmp,outFile);

%Remove temporary file
delete(outFileTmp);
    

