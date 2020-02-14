%function []=rmsrc(dataDir,outDir,srchParamDir,FileName,searchParamName,ext,bandNum,numSrc)
%[]=rmsrc(dataDir,FileName,searchParamName,ext,bandNum,numSrc)
% Subtract sources from Max/Min SNR  recursively and output a new data file.
% dataDir: Input data dirctory.
% outDir: Output directory.
% srchParamDir: search parameter file directory.
% fileName: file wants to be subtracted without extension.
% searchParamName: search parameter file name.
% ext: File type.
% bandNum: Where the sources are going to be removed.
% numSrc: number of sources want to be subtracted.

%2019.11.05 QYQ
clear
tic
%% Test
dataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands/';
outDataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands/RMSubsets/Int100';
mkdir(outDataDir);
srchParamDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/searchParams/Whole';
FileName = 'GWBsimDataSKASrlz1Nrlz3';
searchParamName = 'searchParams';
ext = '.mat';
bandNum = 1;

SNRcut = [0 100];% Start removing sources between this SNR interval.
SNRmax = 700; % maximum SNR
SNRint = SNRcut(2) - SNRcut(1);
steps = (SNRmax - SNRcut(1))/SNRint;

inputFile = strcat(dataDir,filesep,FileName,ext);
% c = 0;
for lp = 1:steps  % Maximum SNR
%     c = c+1;
    load(inputFile);
    searchParams = strcat(srchParamDir,filesep,searchParamName,num2str(bandNum),ext);
    load(searchParams);
    Index = find(omega >= searchParams.angular_velocity(2) & ...
        omega <= searchParams.angular_velocity(1));
    binsrcsnr = snr_chr(Index);
    Np = simParams.Np;
    N = simParams.N;
    I = find(SNRcut(1) <= binsrcsnr & binsrcsnr <= SNRcut(2));
    ite = length(I);
    %I = 33; % debug
    %% calcu. timing residual
    phiI = zeros(Np,1);
    snr_chr2_tmp=zeros(Np,1);  % squared characteristic snr for each pulsar and source
    tmp=zeros(1,N); % noiseless timing residuals from a source
    timingResiduals_tmp = zeros(Np,N);
    timingResiduals_tmp1 = zeros(Np,N);
    for j = 1:ite
        for i=1:1:Np  % number of pulsar
            % GW sky location in Cartesian coordinate
            k=zeros(1,3);  % unit vector pointing from SSB to source
            k(1)=cos(delta(Index(I(j))))*cos(alpha(Index(I(j))));
            k(2)=cos(delta(Index(I(j))))*sin(alpha(Index(I(j))));
            k(3)=sin(delta(Index(I(j))));
            theta=acos(k*simParams.kp(i,:)');
            phiI(i)=mod(phi0(Index(I(j)))-0.5*omega(Index(I(j)))*simParams.distP(i)*(1-cos(theta)), pi); % modulus after division, YW 04/30/14 check original def. of phiI
            
            
            tmp = FullResiduals(alpha(Index(I(j))),delta(Index(I(j))),omega(Index(I(j))),phi0(Index(I(j))),phiI(i),simParams.alphaP(i),...
                simParams.deltaP(i),Amp(Index(I(j))),iota(Index(I(j))),thetaN(Index(I(j))),theta,yr);
            
            timingResiduals_tmp(i,:) = tmp';
            snr_chr2_tmp(i,1) = dot(tmp,tmp)/simParams.sd(i)^2;
            
        end
        
        snr_chr_tmp=sqrt(sum(snr_chr2_tmp,1));  % sum of elements in each column
        
        %% substitute timing residual and SNR
        timingResiduals_tmp1 = timingResiduals_tmp1 + timingResiduals_tmp;
        snr_chr(Index(I(j))) = snr_chr(Index(I(j))) - snr_chr_tmp;
    end
    
    newFile = strcat(outDataDir,filesep,FileName,'_rm_SNR_',num2str(SNRcut(1)),'_',num2str(SNRcut(2)),ext);
    copyfile(inputFile,newFile);
    m = matfile(newFile,'Writable',true);
    m.timingResiduals = m.timingResiduals - timingResiduals_tmp1;
    m.snr_chr = snr_chr;
    SNRcut = SNRcut + SNRint;
end
toc
% EOF