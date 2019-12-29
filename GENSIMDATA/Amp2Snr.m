function [snr_chr,timingResiduals]=Amp2Snr(srcParams,psrParams,yr)
% A function to convert amplitude to SNR
% [SNR,timingResiduals]=Amp2Snr(sourceParams,pulsarParams)

% QYQ 9th April, 2019
%% calculate SNR
%phiI = psrParams.phiI;
Np = psrParams.Np;% number of pulsars
N = psrParams.N;% observation period 5 yr biweekly
snr_chr2_tmp=zeros(Np,1);  % squared characteristic snr for each pulsar and source
tmp=zeros(1,N); % noiseless timing residuals from a source
timingResiduals = zeros(Np,N);
for i=1:1:Np  % number of pulsar
        % GW sky location in Cartesian coordinate
        k=zeros(1,3);  % unit vector pointing from SSB to source
        k(1)=cos(srcParams.delta)*cos(srcParams.alpha);
        k(2)=cos(srcParams.delta)*sin(srcParams.alpha);
        k(3)=sin(srcParams.delta);
        theta=acos(k*psrParams.kp(i,:)');
        %phiI(j)=mod(srcParams.phi0-0.5*srcParams.omega*psrParams.distP(j)*(1-cos(theta)), pi);  % modulus after division, YW 04/30/14 check original def. of phiI
        
        
        tmp = FullResiduals(srcParams.alpha,srcParams.delta,srcParams.omega,srcParams.phi0,srcParams.phiI(i),psrParams.alphaP(i),...
            psrParams.deltaP(i),srcParams.Amp,srcParams.iota,srcParams.thetaN,theta,yr);

        timingResiduals(i,:) = tmp';
        snr_chr2_tmp(i,1) = dot(tmp,tmp)/psrParams.sd(i)^2;

end

snr_chr=sqrt(sum(snr_chr2_tmp,1));  % sum of elements in each column
% END of function