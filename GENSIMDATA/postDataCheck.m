inFileList = dir(['~/Research/PulsarTiming/SimDATA/MultiSource/Multi_Results_MAT/Srlz1Nrlz1',filesep,'*.mat']);
nFiles = length(inFileList);
N = 10; % number of bands
bestRealLoc = zeros(7,N); % collect the first 7 entries from the output file
xmaxmin = zeros(7,2,N); % search parameters used in pso
freq = zeros(10,1);
for lpc = 1:nFiles
    inFile = ['~/Research/PulsarTiming/SimDATA/MultiSource/Multi_Results_MAT/Srlz1Nrlz1',filesep,inFileList(lpc).name];
    inParam = ['~/Research/PulsarTiming/SimDATA/MultiSource/searchParams_GWBsimDataSKA',filesep,inFileList(lpc).name];
    postfile = load(inFile);
    param = load(inParam);
    bestRealLoc(:,lpc) = postfile.bestRealLoc(1:7);
    freq(lpc,1) = bestRealLoc(3,lpc)/(2*pi*31536000); % convert the angular frequency rad/yr to herz
    xmaxmin(:,:,lpc) = param.xmaxmin;
end
OutFile = 'summary.mat';
save(OutFile,'bestRealLoc','freq','xmaxmin')