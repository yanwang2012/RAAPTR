%% Subtract specific sources
% Given one specific parameter set, subtract it from the simulation data

%2019.11.05 QYQ
clear;
tic
%% load sim data
dataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test10/FullBand';
dataFile = 'GWBsimDataSKASrlz1Nrlz3';
load([dataDir,filesep,dataFile,'.mat']);
Np = simParams.Np;
N = simParams.N;
[M,I] = max(snr_chr);
%% calcu. timing residual
phiI = zeros(Np,1);
snr_chr2_tmp=zeros(Np,1);  % squared characteristic snr for each pulsar and source
tmp=zeros(1,N); % noiseless timing residuals from a source
timingResiduals_tmp = zeros(Np,N);
for i=1:1:Np  % number of pulsar
        % GW sky location in Cartesian coordinate
        k=zeros(1,3);  % unit vector pointing from SSB to source
        k(1)=cos(delta(I))*cos(alpha(I));
        k(2)=cos(delta(I))*sin(alpha(I));
        k(3)=sin(delta(I));
        theta=acos(k*simParams.kp(i,:)');
        phiI(i)=mod(phi0(I)-0.5*omega(I)*simParams.distP(i)*(1-cos(theta)), pi); % modulus after division, YW 04/30/14 check original def. of phiI
        
        
        tmp = FullResiduals(alpha(I),delta(I),omega(I),phi0(I),phiI(i),simParams.alphaP(i),...
            simParams.deltaP(i),Amp(I),iota(I),thetaN(I),theta,yr);

        timingResiduals_tmp(i,:) = tmp';
        snr_chr2_tmp(i,1) = dot(tmp,tmp)/simParams.sd(i)^2;

end

snr_chr_tmp=sqrt(sum(snr_chr2_tmp,1));  % sum of elements in each column

%% substitute timing residual
timingResiduals = timingResiduals - timingResiduals_tmp;
newFile = strcat(dataFile,'_rm','.mat');
copyfile([dataDir,filesep,dataFile,'.mat'],[dataDir,filesep,newFile]);
m = matfile([dataDir,filesep,newFile],'Writable',true);
m.timingResiduals = timingResiduals;
m.snr_chr(1,I) = snr_chr(I) - snr_chr_tmp;

toc
% EOF