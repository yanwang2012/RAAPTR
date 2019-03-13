function [] = matfile2hdf5(inFileName, outFileName)
% Convert input data file to HDF5 format (C convention for matrices)
%MATFILE2HDF5(I,O)
%I is the name of the input .mat file. O is the name of the .h5 file. Only
%scalars, vectors, arrays and strings are transported. Cell and struct
%arrays are ignored. 

%Soumya D. Mohanty, Jan 2016
%Adapted from extinct code MATFILE2HDF5; Now handles strings properly

%How many variables in inFile; Their names
inFileInfo = load(inFileName);
fldNames = fieldnames(inFileInfo);
nVars = length(fldNames);

%Create HDF5 file
fid = H5F.create(outFileName);
H5F.close(fid);

%Loop over variables
for lp = 1:nVars
    %Check: Vectors, matrices, strings only
    isn = isnumeric(inFileInfo.(fldNames{lp}));
    isca = ischar(inFileInfo.(fldNames{lp}));
    if (isn || isca)
        %Store variable in .h5 file; store in row major order -> store transpose
        if(isca)
            %Storing a string requires using the HDF5 package functions. This
            %code is adapted from the examples in the documentation for
            %H5D functions.
            fid = H5F.open(outFileName, 'H5F_ACC_RDWR','H5P_DEFAULT');
%             fid = H5F.create('myfile_strings.hdf5');
%             data = 'soumya.d.mohanty';
            type_id = H5T.copy('H5T_C_S1');
            h5_dims = 1;
            h5_maxdims = h5_dims;
            space_id = H5S.create_simple(1,h5_dims,h5_maxdims);
            dcpl = 'H5P_DEFAULT';
            plist = 'H5P_DEFAULT';
            outStr = inFileInfo.(fldNames{lp});
            if ~isempty(outStr)
                H5T.set_size(type_id,length(outStr));
                %             H5T.set_size(type_id,length(data));
            else
                %Empty string
                H5T.set_size(type_id,1);
                outStr = ' ';
            end           
            dset_id = H5D.create(fid,fldNames{lp},type_id,space_id,dcpl);
            H5D.write(dset_id,'H5ML_DEFAULT','H5S_ALL','H5S_ALL',...
                plist,outStr);
            %             H5D.write(dset_id,'H5ML_DEFAULT','H5S_ALL','H5S_ALL',plist,data);
            H5S.close(space_id);
            H5D.close(dset_id);
            H5T.close(type_id);
            H5F.close(fid);
        else
            nDim = length(size((inFileInfo.(fldNames{lp}))'));
            baseSz = ones(1, nDim);
            h5create(outFileName,['/',fldNames{lp}],max([size((inFileInfo.(fldNames{lp}))'); baseSz]));
            if ~isempty(inFileInfo.(fldNames{lp}))
                h5write(outFileName,['/',fldNames{lp}],(inFileInfo.(fldNames{lp}))');
            else
                h5write(outFileName,['/',fldNames{lp}],NaN);
            end
        end
    else
        disp(['ignoring non-numeric/string variable ',fldNames{lp}]);
    end

end