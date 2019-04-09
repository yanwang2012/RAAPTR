% core function of converting Amp to SNR
function [snr_chr]=convertAmp2snr(sourceParams,pulsarParams)
src = sourceParams;
psr = pulsarParams;
%% calculate SNR
phiI = psr.phiI;
Np = psr.Np;% number of pulsars
N = psr.N;
snr_chr2_tmp=zeros(Np,1);  % squared characteristic srn for each pulsar and source
tmp=zeros(1,N); % noiseless timing residuals from a source
for j=1:1:Np  % number of pulsar
    %for j=1:1:Ns  % number of GW source
        
        % GW sky location in Cartesian coordinate
        k=zeros(1,3);  % unit vector pointing from SSB to source
        k(1)=cos(src.delta)*cos(src.alpha);
        k(2)=cos(src.delta)*sin(src.alpha);
        k(3)=sin(src.delta);
        theta=acos(k*psr.kp(j,:)');
        phiI(j)=mod(src.phi0-0.5*src.omega*psr.distP(j)*(1-cos(theta)), pi);  % modulus after division, YW 04/30/14 check original def. of phiI
        
        
        tmp = FullResiduals(src.alpha,src.delta,src.omega,src.phi0,phiI(j),psr.alphaP(j),...
            psr.deltaP(j),src.Amp,src.iota,src.thetaN,theta,psr.yr);

        
        snr_chr2_tmp(j,1) = dot(tmp,tmp)/psr.sd(j)^2;

end

snr_chr=sqrt(sum(snr_chr2_tmp,1));  % sum of elements in each column