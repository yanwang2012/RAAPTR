function [snr_chr]=Amp2Snr(srcParams,psrParams,phiI,yr)
% A function to convert amplitude to SNR
% [SNR]=Amp2Snr(sourceParams,pulsarParams)

% QYQ 9th April, 2019
%% calculate SNR
%phiI = psrParams.phiI;
Np = psrParams.Np;% number of pulsars
N = psrParams.N;
snr_chr2_tmp=zeros(Np,1);  % squared characteristic srn for each pulsar and source
tmp=zeros(1,N); % noiseless timing residuals from a source
for j=1:1:Np  % number of pulsar
    %for j=1:1:Ns  % number of GW source
        
        % GW sky location in Cartesian coordinate
        k=zeros(1,3);  % unit vector pointing from SSB to source
        k(1)=cos(srcParams.delta)*cos(srcParams.alpha);
        k(2)=cos(srcParams.delta)*sin(srcParams.alpha);
        k(3)=sin(srcParams.delta);
        theta=acos(k*psrParams.kp(j,:)');
        %phiI(j)=mod(srcParams.phi0-0.5*srcParams.omega*psrParams.distP(j)*(1-cos(theta)), pi);  % modulus after division, YW 04/30/14 check original def. of phiI
        
        
        tmp = FullResiduals(srcParams.alpha,srcParams.delta,srcParams.omega,srcParams.phi0,phiI(j),psrParams.alphaP(j),...
            psrParams.deltaP(j),srcParams.Amp,srcParams.iota,srcParams.thetaN,theta,yr);

        
        snr_chr2_tmp(j,1) = dot(tmp,tmp)/psrParams.sd(j)^2;

end

snr_chr=sqrt(sum(snr_chr2_tmp,1));  % sum of elements in each column