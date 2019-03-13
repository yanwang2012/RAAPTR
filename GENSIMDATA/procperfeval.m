%% Post-processing of the output from PERFEVAL
% Process the the output files from PERFEVAL to produce summary
% results and store them in the files described below.
%
% * 'AllInfo': Relevant information extracted from each file and stored in
%              a struct array. There is one struct array element per
%              PERFEVAL output file. The structure fields contain the
%              following info.
%  *            Relevant meta data from associated output file 
%  *            Fitness value for best of M runs for each noise realization. 
%  *            Min and max fitness values over M runs for each noise realization.
%  *            Estimated signal parameters for best of M runs for each noise realization.
%  *            Mean root deviation of each parameter from best value over M runs for each noise realization.
%  *            Number of function evaluations for best of M runs for each noise realization.
%  *            Min and max of number of function evaluations over M runs for each noise realization.
%  *            Number of iterations for best of M runs for each noise realization.
%  *            Min and max of number of iterations over M runs for each noise realization. 
%
%  * 'H0stats': consolidated info across noise realizations from output files containing tag 'H0'.
%  * 'H1stats_<*>': consolidated info across noise realizations from output files containing tag 'H1' and signal 
%   parameter id tag *.
%  Each file contains the following info.
%  *            List of PERFEVAL output files used
%  *            Best fitness value for each noise realization.
%  *            Best estimated signal parameters for
%               each noise realization.
%  *            Number of functions evaluations used in the best run for
%               each noise realization.
% 
%
% * 'DetStats_<*>': This file contains detection performance related
%    information (Effective SNR and ROC)
% * 'EstStats_<*>': Estimation performance related information (Mean and
%    covariance of estimated parameters)
%  

%% Read files and gather info
%Relevant data from input files is loaded into elements of a struct array.
%The struct array elements are then processed to produce the desired output
%files.
%simDataDir = 'simData1';
%simDataDir = '/Users/ywang/Research/PULSARTIMING/MultiCW/simData17_snr123_loc3/';
%simDataDir = '/Users/ywang/Research/PULSARTIMING/AvPhase/AvMax_MatlabCodes/AvPhase/simData17_snr123_loc3/';
%simDataDir = '/Users/ywang/Research/PULSARTIMING/MultiCW/simDataSKA_snr0123_loc1_test/';
%simDataDir = '/Users/ywang/Research/PULSARTIMING/MultiCW/simDataSKA_snr0123_loc3/';
simDataDir = '/Users/ywang/Research/PULSARTIMING/MultiCW/simDataSKA_loc6/';
%simDataDir = '/Users/ywang/Research/PULSARTIMING/MultiCW/simDataSKA_loc_PG1302/';
%simDataDir = '/Users/ywang/Research/PULSARTIMING/MultiCW/simData9_snr0123_loc3/';

outDir = [simDataDir,filesep,'results_LS5_simDataSKA_loc6_avPhase',filesep,'summary'];
mkdir(outDir);
inDataDir = [simDataDir,filesep,'results_LS5_simDataSKA_loc6_avPhase'];
% List of input data files (we will only process .mat files)
inputFiles = dir([inDataDir,filesep,'*.mat']);
% Number of input Data files 
nInputFiles = length(inputFiles);
%% Output struct array
outStruct = struct('dummy',cell(1,nInputFiles),...
                   'metaData',struct('fileName','',...
                                     'genHypothesis','',...
                                     'id',struct()...
                                     ),...
                   'fitnessBOR',[],...
                   'minFitness',[],...
                   'maxFitness',[],...
                   'estSigParams',[],...
                   'rssDevSigParams',[],...
                   'funcEvalsBOR',[],...
                   'minFuncEvals',[],...
                   'maxFuncEvals',[],...
                   'nIterBOR',[],...
                   'minIter',[],...
                   'maxIter',[]...
                   );
%% Fill structure array
for lpfiles = 1:nInputFiles
    inFileName = inputFiles(lpfiles).name;
    inFileContents = load([inDataDir,filesep,inFileName]);
    nRuns = inFileContents.nRuns;
    outStruct(lpfiles).metaData.fileName = [inDataDir,filesep,inFileName];
    outStruct(lpfiles).metaData.genHypothesis = inFileContents.genHypothesis;
    outStruct(lpfiles).metaData.id = inFileContents.id;
    bestRun = inFileContents.bestRun;
    fitnessVals = inFileContents.fitnessVals;
    outStruct(lpfiles).fitnessBOR = fitnessVals(bestRun);
    outStruct(lpfiles).minFitness = min(fitnessVals);
    outStruct(lpfiles).maxFitness = max(fitnessVals);
    bestLocVals = inFileContents.bestLocVals;
    bestLoc = bestLocVals(bestRun,:);
    outStruct(lpfiles).estSigParams = inFileContents.bestRealLoc;
    %TODO complete the calculation of std dev once the real coordinates for
    %all runs are being stored.
    %outStruct(lpfiles).rssDevSigParams = sqrt((1/(max(nRuns-1,1)))*...
    %                                          sum(inFileContents.bestLocRealCVals-
    funcEvals = inFileContents.numFitEvals;
    outStruct(lpfiles).funcEvalsBOR = funcEvals(bestRun);
    outStruct(lpfiles).minFuncEvals = min(funcEvals);
    outStruct(lpfiles).maxFuncEvals = max(funcEvals);
    %TODO store number of iterations so that this field can be filled.
    %numIterVec = inFileContents.numIter;
    %outStruct(lpfiles).nIterBOR = numIterVec(bestRun);
    %outStruct(lpfiles).minIter = min(numIterVec);
    %outStruct(lpfiles).maxIter = max(numIterVec);
end
save([outDir,filesep,'allInfo'],'outStruct');

%%
% Identify distinct types of files (noise/signal and signal id's). 
% 'elements': list of file numbers that match the current id (snr, loc and
% omg should match).
fileTypes = struct(...
                   'elements',1,...
                   'id',outStruct(1).metaData.id,...
                   'hypothesis',outStruct(1).metaData.genHypothesis);
typeCount = 1;
%Match each file against all current file types. If a match happens, add
%the file to corresponding file type's list. If not, increment the number
%of file types.
for lpfiles = 2:nInputFiles
    for lptype = 1:typeCount
        noMatch = 1;
        if outStruct(lpfiles).metaData.id.snr_id==fileTypes(lptype).id.snr_id && ...
           outStruct(lpfiles).metaData.id.loc_id==fileTypes(lptype).id.loc_id && ...
           outStruct(lpfiles).metaData.id.omg_id==fileTypes(lptype).id.omg_id
               fileTypes(lptype).elements = [fileTypes(lptype).elements,lpfiles];
               if ~strcmp(fileTypes(lptype).hypothesis, outStruct(lpfiles).metaData.genHypothesis)
                   warning(['Something wrong with ', outStruct(lpfiles).metaData.fileName]);
               end
               noMatch = 0;
               break;
        end
    end
    if noMatch
        typeCount = typeCount+1;
        fileTypes(typeCount).id = outStruct(lpfiles).metaData.id;
        fileTypes(typeCount).hypothesis = outStruct(lpfiles).metaData.genHypothesis;
        fileTypes(typeCount).elements = [fileTypes(typeCount).elements,lpfiles];
    end
end
%% Detection and estimation results
% List of files containing null and alternative hypotheses results.
h0File = '';
h1Files = {};
for lptype = 1:length(fileTypes)
    hypType = fileTypes(lptype).hypothesis;
    structElements = fileTypes(lptype).elements;
    dummyStruct = struct('FilesUsed',{},...
                'id',struct(),...
                'fitnessVals',[],...
                'estSigParams',[]...
                );
    for lpel = 1:length(structElements)
        dummyStruct(1).FilesUsed = [dummyStruct.FilesUsed,...
                                   {outStruct(structElements(lpel)).metaData.fileName}];
        dummyStruct(1).fitnessVals = [dummyStruct.fitnessVals;...
                                   outStruct(structElements(lpel)).fitnessBOR];
        dummyStruct(1).estSigParams = [dummyStruct.estSigParams;...
                                   outStruct(structElements(lpel)).estSigParams];        
    end
    dummyStruct.id = fileTypes(lptype).id;
    idTag = [...
             'snr',num2str(dummyStruct.id.snr_id),...
             '_loc',num2str(dummyStruct.id.loc_id),...
             '_omg',num2str(dummyStruct.id.omg_id)...
             ];
    switch hypType
        case 'H1 data'
            outSumryFile = [outDir,filesep,'H1stats_',idTag];
            h1Files = [h1Files,{outSumryFile}];
        case 'H0 data'
            outSumryFile = [outDir,filesep,'H0stats'];
            h0File = outSumryFile;
    end
    save(outSumryFile,'-struct','dummyStruct');
    
end
nh1Files = length(h1Files);

%% Effective SNR and other metrics
%  $X_{H1/H0}$ detection statistic (negative of fitness
% values) under the alternative and null hypotheses.
% Effective SNR:
%
% $$\frac{\overline{X}_{H1}-\overline{X}_{H0}}{\sqrt{{\rm var}(X_{H0})}}$$
%
if ~isempty(h0File)
    xh0 = load(h0File,'fitnessVals');
    meanXh0 = mean(-xh0.fitnessVals);
    stdevXh0 = std(-xh0.fitnessVals);
end
for lpfiles = 1:nh1Files
    h1Dat = load(h1Files{lpfiles});
    % Detection metrics
    xh1 = h1Dat.fitnessVals;
    meanXh1 = mean(-xh1);
    stdevXh1 = std(-xh1);
    if ~isempty(h0File)
     effSNR = (meanXh1-meanXh0)/stdevXh0;
    else
        effSNR = NaN;
    end
    % Estimation metrics
    estParams = h1Dat.estSigParams;
    meanParams = mean(estParams);
    stdevParams = std(estParams);
    % Store metrics
    id = h1Dat.id;
    idTag = [...
             'snr',num2str(id.snr_id),...
             '_loc',num2str(id.loc_id),...
             '_omg',num2str(id.omg_id)...
             ];
    outFile = [outDir,filesep,'DetStats_',idTag];
    save(outFile,'meanXh1','stdevXh1','effSNR',...
        'meanParams','stdevParams',...
        'id');
end


