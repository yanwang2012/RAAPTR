% Removing sources in simulation data to investigate the performance of
% PSO.
% 2019.11.05 QYQ
clear
tic
%% Dir settings
dataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands';
outDataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands/FreqRM';
mkdir(outDataDir);
srchParamDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/searchParams/2bands/superNarrow';
FileName = 'GWBsimDataSKASrlz1Nrlz3';
searchParamName = 'searchParams';
ext = '.mat';

%% Config
bandNum = 1; % band where the source needs to be removed
numSrc = 5; % number of sources need to be removed
srcN = 2; % The starting source.
startpoint = srcN; % set up the

inputFile = strcat(dataDir,filesep,FileName,ext);
simFile = inputFile;
newFileNameC = {};
newFileNameI = {};

%% Main
load(inputFile);
searchParams = strcat(srchParamDir,filesep,searchParamName,num2str(bandNum),ext);
load(searchParams);
Index = find(omega >= searchParams.angular_velocity(2) & ...
    omega <= searchParams.angular_velocity(1));
binsrcsnr = snr_chr(Index);
Np = simParams.Np;
N = simParams.N;
[~,sortedIndex] = sort(binsrcsnr,'descend'); % sort sources according to SNR.

%% calcu. timing residual
snr_chr_acu = snr_chr; % Removing sources continuously.
snr_chr_indv = snr_chr; % Removing source individually.
phiI = zeros(Np,1);
snr_chr2_tmp=zeros(Np,1);  % squared characteristic snr for each pulsar and source
tmp=zeros(1,N); % noiseless timing residuals from a source
timingResiduals_tmp = zeros(Np,N);
timingResiduals_tmpc = zeros(Np,N); % accumulated 
timingResiduals_tmpi = zeros(Np,N); % individually

for lp = 1:numSrc
    for j = 1:lp
        for i=1:1:Np  % number of pulsar
            % GW sky location in Cartesian coordinate
            k=zeros(1,3);  % unit vector pointing from SSB to source
            k(1)=cos(delta(Index(sortedIndex(startpoint))))*cos(alpha(Index(sortedIndex(startpoint))));
            k(2)=cos(delta(Index(sortedIndex(startpoint))))*sin(alpha(Index(sortedIndex(startpoint))));
            k(3)=sin(delta(Index(sortedIndex(startpoint))));
            theta=acos(k*simParams.kp(i,:)');
            phiI(i)=mod(phi0(Index(sortedIndex(startpoint)))-0.5*omega(Index(sortedIndex(startpoint)))*simParams.distP(i)*(1-cos(theta)), pi); % modulus after division, YW 04/30/14 check original def. of phiI
            
            
            tmp = FullResiduals(alpha(Index(sortedIndex(startpoint))),delta(Index(sortedIndex(startpoint))),omega(Index(sortedIndex(startpoint))),phi0(Index(sortedIndex(startpoint))),phiI(i),simParams.alphaP(i),...
                simParams.deltaP(i),Amp(Index(sortedIndex(startpoint))),iota(Index(sortedIndex(startpoint))),thetaN(Index(sortedIndex(startpoint))),theta,yr);
            
            timingResiduals_tmp(i,:) = tmp';
            snr_chr2_tmp(i,1) = dot(tmp,tmp)/simParams.sd(i)^2;
            
        end
        
        snr_chr_tmp=sqrt(sum(snr_chr2_tmp,1));  % sum of elements in each column
        
        %% substitute timing residual
        timingResiduals_tmpc = timingResiduals_tmpc + timingResiduals_tmp; % accumulate timing residuals, remove all the sources from start point.
        timingResiduals_tmpi = timingResiduals - timingResiduals_tmp; % only remove one source each time.
        snr_chr_acu(Index(sortedIndex(startpoint))) = snr_chr(Index(sortedIndex(startpoint))) - snr_chr_tmp;
        if j == lp
            snr_chr_indv(Index(sortedIndex(startpoint))) = snr_chr(Index(sortedIndex(startpoint))) - snr_chr_tmp;
        end
        startpoint = startpoint + 1;
    end
    newFileC = strcat(outDataDir,filesep,FileName,'_rm',num2str(startpoint-1),ext);
    newFileI = strcat(outDataDir,filesep,FileName,'_rm',num2str(startpoint-1),'only',ext);
    newFileNameC = [newFileNameC newFileC];
    newFileNameI = [newFileNameI newFileI];
    copyfile(inputFile,newFileC);
    copyfile(inputFile,newFileI);
    m = matfile(newFileC,'Writable',true);
    mm = matfile(newFileI,'Writable',true);
    m.timingResiduals = timingResiduals - timingResiduals_tmpc;
    mm.timingResiduals = timingResiduals_tmpi;
    m.snr_chr = snr_chr_acu;
    mm.snr_chr = snr_chr_indv;
    % reinitialization
    startpoint = srcN; % initialize startpoint again.
    snr_chr_indv = snr_chr;
    timingResiduals_tmpc = zeros(Np,N);
end


%% Plot check
figure(1)
ori = load(simFile);
new = load(newFileNameC{4});

ox = ori.snr_chr;
oy = ori.omega/(2*pi*365*24*3600);

nx = new.snr_chr;
ny = new.omega/(2*pi*365*24*3600);

plot(ox,oy,'or',nx,ny,'sb');
legend('Origin','Removed')

figure(2)

ori = load(simFile);
new = load(newFileNameI{4});

ox = ori.snr_chr;
oy = ori.omega/(2*pi*365*24*3600);

nx = new.snr_chr;
ny = new.omega/(2*pi*365*24*3600);

plot(ox,oy,'or',nx,ny,'sb');
legend('Origin','Removed')



toc
% END