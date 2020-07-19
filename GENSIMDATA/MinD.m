function [d,id_tsrc] = MinD(band,NtsrcBand,SrcAlpha,SrcDelta,SrcOmega,SrcPhi0,...
    SrcIota,SrcThetaN,SrcAmp,EstSrc,simParams,yr)
% A function returns the index of true source which has the minimum
% distance to the estimated source.
% [id_tsrc] = MinD(band,SrcAlpha,SrcDelta,SrcOmega,SrcPhi0,...
%    SrcIota,SrcThetaN,SrcAmp,EstSrc,simParams,yr)
% band: number of bands
% NtsrcBand: number of true sources in the band
% SrcX: parameters of true sources
% EstSrc: parameters of estimated sources
% simParams: pulsar parameters
% yr: observation span
% returns: 
% d: distance array for each est. source
% id_tsrc: the index of the closest true source w.r.t the est. source

% Author: QYQ 07/16/2020



%% Pulsar configuration
Np = simParams.Np; % number of pulsars
deltaP = simParams.deltaP;
alphaP = simParams.alphaP;
kp = simParams.kp;
distP = simParams.distP;

% per-pulsar distance for est. sources and true sources.
dp = zeros(NtsrcBand,Np);

d = zeros(NtsrcBand,1); % averaged per-pulsar distance for each true source.



for tsrc = 1:NtsrcBand
    for psr = 1:Np
        % GW sky location in Cartesian coordinate
        k=zeros(1,3);  % unit vector pointing from SSB to source
        k(1)=cos(SrcDelta{band}(tsrc))*cos(SrcAlpha{band}(tsrc));
        k(2)=cos(SrcDelta{band}(tsrc))*sin(SrcAlpha{band}(tsrc));
        k(3)=sin(SrcDelta{band}(tsrc));
        theta=acos(k*kp(psr,:)');
        %sprintf('%d pulsar theta=%g',i,theta)
        %phiI(i)=mod(phi0-omega*distP(i)*(1-cos(theta)), 2*pi);  % modulus after division
        %phiI(i)=mod(2*phi0-omega_tmp(l)*distP(i)*(1-cos(theta)), pi);  % modulus after division, YW 09/10/13
        phiI(psr)=mod(SrcPhi0{band}(tsrc)-0.5*SrcOmega{band}(tsrc)*distP(psr)*(1-cos(theta)), pi);  % modulus after division, YW 04/30/14 check original def. of phiI
        
        tmp = FullResiduals(SrcAlpha{band}(tsrc),SrcDelta{band}(tsrc),SrcOmega{band}(tsrc),SrcPhi0{band}(tsrc),phiI(psr),alphaP(psr),deltaP(psr),...
            SrcAmp{band}(tsrc),SrcIota{band}(tsrc),SrcThetaN{band}(tsrc),theta,yr); % timing residuals for true src
        
        tmp_est = FullResiduals(EstSrc.alpha,EstSrc.delta,EstSrc.omega,EstSrc.phi0,EstSrc.phiI(psr),alphaP(psr),deltaP(psr),...
            EstSrc.Amp,EstSrc.iota,EstSrc.thetaN,theta,yr); % timing residuals for estimated source
        
        dp(tsrc,psr) = norm(tmp-tmp_est); % per-pulsar distance.
    end
end

d = sum(dp,2)/Np;

[~,id_tsrc] = min(d);

% EOF
