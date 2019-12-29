function outStruct = loadhdf5(inFile, varargin)
%Load data from an HDF5 file into workspace
%S = LOADHDF5(FILENAME) loads the datasets from an HDF5 file into a
%structure array. Currently only datasets stored under the root ('/') group
%are retrieved.
%
%S = LOADHDF5(FILENAME,VARIABLES) loads only the specified datasets from an
%HDF5 file.

%Soumya D. Mohanty, Jan 2017

%Handle file name without .hdf5 extension
[~,~,fileExt] = fileparts(inFile);
if isempty(fileExt)
    fileExt = '.hdf5';
else
    fileExt = '';
end
inFile = [inFile,fileExt];

fileInfo = h5info(inFile);

%Number of variables to read 
nVar = length(fileInfo.Datasets);

outStruct = struct();
%Count successful reads
countGoodReads = 0;
for lp = 1:nVar
    %Get the names of the variables
    varName = fileInfo.Datasets(lp).Name;
    %Default action is to not read the associated data
    readData = 0;
    if nargin > 1
        for lp2 = 1:(nargin-1)
            if strcmp(varName,varargin{lp2})
                %Read data
                readData = 1;
                countGoodReads = countGoodReads+1;
                break;
            end
        end
    else
        readData = 1;
        countGoodReads = countGoodReads+1;
    end   
    if readData
        %Get the associated data
        varData = h5read(inFile,['/',varName]);
        %Undo matlab's column ordering default
        if isnumeric(varData)
            varData = varData';
        end
        %For some reason, HD5 string appears with a 
        %hidden character when loaded by h5read.
        if ischar(varData)
            varData = varData(1:(end-1));
        end
        outStruct = setfield(outStruct,varName,varData);
    end
end
%Send out empty vector if nothing was read successfully
if ~countGoodReads
    outStruct = [];
end
