function []=rmsrc(dataDir,srchParamDir,FileName,searchParamName,ext,bandNum,numSrc)
%[]=rmsrc(dataDir,FileName,searchParamName,ext,bandNum,numSrc)
% Subtract sources according to SNR recusively and output a new file.
% dataDir: data dirctory.
% srchParamDir: search parameter file directory.
% fileName: file wants to be subtracted without extension.
% searchParamName: search parameter file name.
% ext: File type.
% bandNum: Band number.
% numSrc: number of sources want to be subtracted.

%2019.11.05 QYQ
tic
%% Test
% dataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test10/FullBand';
% srchParamDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test10/searchParams';
% FileName = 'GWBsimDataSKASrlz1Nrlz3';
% searchParamName = 'searchParams_Nyquist';
% ext = '.mat';
% bandNum = 1;
% numSrc = 5;

inputFile = strcat(dataDir,filesep,FileName,ext);
for j = 1:numSrc
    load(inputFile);
    searchParams = strcat(srchParamDir,filesep,searchParamName,num2str(bandNum),ext);
    load(searchParams);
    Index = find(omega >= searchParams.angular_velocity(2) & ...
        omega <= searchParams.angular_velocity(1));
    binsrcsnr = snr_chr(Index);
    Np = simParams.Np;
    N = simParams.N;
    [~,I] = max(binsrcsnr);
    %I = 33; % debug
    %% calcu. timing residual
    phiI = zeros(Np,1);
    snr_chr2_tmp=zeros(Np,1);  % squared characteristic snr for each pulsar and source
    tmp=zeros(1,N); % noiseless timing residuals from a source
    timingResiduals_tmp = zeros(Np,N);
    for i=1:1:Np  % number of pulsar
        % GW sky location in Cartesian coordinate
        k=zeros(1,3);  % unit vector pointing from SSB to source
        k(1)=cos(delta(Index(I)))*cos(alpha(Index(I)));
        k(2)=cos(delta(Index(I)))*sin(alpha(Index(I)));
        k(3)=sin(delta(Index(I)));
        theta=acos(k*simParams.kp(i,:)');
        phiI(i)=mod(phi0(Index(I))-0.5*omega(Index(I))*simParams.distP(i)*(1-cos(theta)), pi); % modulus after division, YW 04/30/14 check original def. of phiI
        
        
        tmp = FullResiduals(alpha(Index(I)),delta(Index(I)),omega(Index(I)),phi0(Index(I)),phiI(i),simParams.alphaP(i),...
            simParams.deltaP(i),Amp(Index(I)),iota(Index(I)),thetaN(Index(I)),theta,yr);
        
        timingResiduals_tmp(i,:) = tmp';
        snr_chr2_tmp(i,1) = dot(tmp,tmp)/simParams.sd(i)^2;
        
    end
    
    snr_chr_tmp=sqrt(sum(snr_chr2_tmp,1));  % sum of elements in each column
    
    %% substitute timing residual
    timingResiduals_tmp1 = timingResiduals - timingResiduals_tmp;
    newFile = strcat(dataDir,filesep,FileName,'_rm',num2str(j),ext);
    copyfile(inputFile,newFile);
    m = matfile(newFile,'Writable',true);
    m.timingResiduals = timingResiduals_tmp1;
    m.snr_chr(1,Index(I)) = snr_chr(Index(I)) - snr_chr_tmp;
    inputFile = newFile;
end
toc
% EOF