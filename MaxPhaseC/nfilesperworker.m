function nFilesList = nfilesperworker(nWorkers, nFiles)
%List of number of files per worker
%N = NFILESPERWORKER(W, F)
%Generates the start and stop indices of files to assign to each of W workers given F
%files to process. N(i,1) and N(i,2) are the start and end files for the
%i'th worker.

%Soumya D. Mohanty, Jun 2016

remFiles  = mod(nFiles,nWorkers);
baseNumFiles = (nFiles - remFiles)/nWorkers;
nFilesList = zeros(nWorkers,2);
nFilesList(:,2) = baseNumFiles:baseNumFiles:(nFiles-remFiles);
nFilesList(:,1) = nFilesList(:,2) - (baseNumFiles-1);
for lp = 1:remFiles
    nFilesList(lp:nWorkers,2) = nFilesList(lp:nWorkers,2) +1;
    nFilesList((lp+1):nWorkers,1) = nFilesList((lp+1):nWorkers,1) +1;
end
fprintf(1,'%d %d\n',nFilesList');

