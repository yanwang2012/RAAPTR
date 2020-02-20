% Source removing script.
% 2019.11.05 QYQ
clear
tic
%% Dir settings
dataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands';
outDataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands/FreqRM';
srchParamDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/searchParams/2bands/superNarrow';
FileName = 'GWBsimDataSKASrlz1Nrlz3';
searchParamName = 'searchParams';
ext = '.mat';

%% Config
bandNum = 1; % band where the source needs to be removed
numSrc = 5; % number of sources need to be removed
n = 2; % start from nth loudest source.

inputFile = strcat(dataDir,filesep,FileName,ext);
simFile = inputFile;
newFileName = {};

%% Main
for j = 1:numSrc
    load(inputFile);
    searchParams = strcat(srchParamDir,filesep,searchParamName,num2str(bandNum),ext);
    load(searchParams);
    Index = find(omega >= searchParams.angular_velocity(2) & ...
        omega <= searchParams.angular_velocity(1));
    binsrcsnr = snr_chr(Index);
    Np = simParams.Np;
    N = simParams.N;
    %     [~,I] = max(binsrcsnr); % subtract from Max SNR
    [~,sortedIndex] = sort(binsrcsnr,'descend');
    
    %I = 33; % debug
    %% calcu. timing residual
    phiI = zeros(Np,1);
    snr_chr2_tmp=zeros(Np,1);  % squared characteristic snr for each pulsar and source
    tmp=zeros(1,N); % noiseless timing residuals from a source
    timingResiduals_tmp = zeros(Np,N);
    for i=1:1:Np  % number of pulsar
        % GW sky location in Cartesian coordinate
        k=zeros(1,3);  % unit vector pointing from SSB to source
        k(1)=cos(delta(Index(sortedIndex(n))))*cos(alpha(Index(sortedIndex(n))));
        k(2)=cos(delta(Index(sortedIndex(n))))*sin(alpha(Index(sortedIndex(n))));
        k(3)=sin(delta(Index(sortedIndex(n))));
        theta=acos(k*simParams.kp(i,:)');
        phiI(i)=mod(phi0(Index(sortedIndex(n)))-0.5*omega(Index(sortedIndex(n)))*simParams.distP(i)*(1-cos(theta)), pi); % modulus after division, YW 04/30/14 check original def. of phiI
        
        
        tmp = FullResiduals(alpha(Index(sortedIndex(n))),delta(Index(sortedIndex(n))),omega(Index(sortedIndex(n))),phi0(Index(sortedIndex(n))),phiI(i),simParams.alphaP(i),...
            simParams.deltaP(i),Amp(Index(sortedIndex(n))),iota(Index(sortedIndex(n))),thetaN(Index(sortedIndex(n))),theta,yr);
        
        timingResiduals_tmp(i,:) = tmp';
        snr_chr2_tmp(i,1) = dot(tmp,tmp)/simParams.sd(i)^2;
        
    end
    
    snr_chr_tmp=sqrt(sum(snr_chr2_tmp,1));  % sum of elements in each column
    
    %% substitute timing residual
    timingResiduals_tmp1 = timingResiduals - timingResiduals_tmp;
    newFile = strcat(outDataDir,filesep,FileName,'_rm',num2str(j),ext);
    newFileName = [newFileName newFile];
    copyfile(inputFile,newFile);
    m = matfile(newFile,'Writable',true);
    m.timingResiduals = timingResiduals_tmp1;
    %     m.snr_chr(1,Index(I)) = snr_chr(Index(I)) - snr_chr_tmp;
    snr_chr_new = m.snr_chr;
    snr_chr_new(Index(sortedIndex(n))) = snr_chr(Index(sortedIndex(n))) - snr_chr_tmp;
    m.snr_chr = snr_chr_new;
    inputFile = newFile;
    n = n+1;
end
n = n - 1;
disp(['Remove up to ',num2str(n),'th source.']);
%% Plot check
ori = load(simFile);
new = load(newFileName{2});

ox = ori.snr_chr;
oy = ori.omega/(2*pi*365*24*3600);

nx = new.snr_chr;
ny = new.omega/(2*pi*365*24*3600);

plot(ox,oy,'or',nx,ny,'sb');
legend('Origin','Removed')






toc
% END