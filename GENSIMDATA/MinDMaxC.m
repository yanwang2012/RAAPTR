function [rho,rho_max,dif_freq_max,dif_ra_max,dif_dec_max,id_max,estSNR] = MinDMaxC(Nband,NestsrcBand,SrcAlpha,SrcDelta,SrcOmega,SrcPhi0,...
    SrcIota,SrcThetaN,SrcAmp,EstSrc,simParams,yr)
% A function calculates cross-correlation coefficients using max weighted
% algorithm.
% [rho,rho_max,dif_freq_max,dif_ra_max,dif_dec_max,id_max] =
% MinDMaxC(Nband,NestsrcBand,SrcAlpha,SrcDelta,SrcOmega,SrcPhi0,SrcIota,SrcThetaN,SrcAmp,EstSrc,simParams,yr)
% rho: cross-correlaltion coefficient matrix
% rho_max: max cross-correlation coefficients matrix
% weighted parameter.
% dif_freq_max: error in frequency.
% dif_ra_max: error in RA.
% dif_dec_max: error in DEC.
% id_max: index of sources which reach the maximum of coefficients.
% estSNR: SNR value for estimated sources.
% Nband: number of band.
% NestsrcBand: number of estimated source in each band.
% SrcXx: parameters of simulated source.
% EstSrc: cell collects all the parameters of estimated source.
% simParams: pulsar config.
% yr: observation span.

% Author: QYQ 07/16/2020


%% Pulsar configuration
Np = simParams.Np; % number of pulsars
deltaP = simParams.deltaP;
alphaP = simParams.alphaP;
kp = simParams.kp;
distP = simParams.distP;

rho_tmp = zeros(Np,1);

rho = {}; % cross-corelation coefficients matrix
rho_max = {};

% investigations
estSNR = zeros(Nband,NestsrcBand);
dif_freq = {}; % frequency difference
dif_ra = {};
dif_dec = {};
% d = {};

id_max = zeros(NestsrcBand,Nband); % index of max. cc
dif_freq_max = zeros(NestsrcBand,Nband);
dif_ra_max = zeros(NestsrcBand,Nband);
dif_dec_max = zeros(NestsrcBand,Nband);

for band = 1:Nband
    NtsrcBand = length(SrcAlpha{band}); % number of true sources in each band
    rho_tmp2 = zeros(NestsrcBand,NtsrcBand);
    for src = 1:NestsrcBand
        [snr,~] = Amp2Snr(EstSrc{band,src},simParams,yr); % get SNR for estimated source
        estSNR(band,src) = snr;
        [~,id_tsrc] = MinD(band,NtsrcBand,SrcAlpha,SrcDelta,SrcOmega,SrcPhi0,...
            SrcIota,SrcThetaN,SrcAmp,EstSrc{band,src},simParams,yr); % get the index of true source whic is closest to the est. source.
%         disp(["Band: ",band, "Est. source: ",src, "True source match by: ",id_tsrc]);
        for psr = 1:Np
            % GW sky location in Cartesian coordinate
            k=zeros(1,3);  % unit vector pointing from SSB to source
            k(1)=cos(SrcDelta{band}(id_tsrc))*cos(SrcAlpha{band}(id_tsrc));
            k(2)=cos(SrcDelta{band}(id_tsrc))*sin(SrcAlpha{band}(id_tsrc));
            k(3)=sin(SrcDelta{band}(id_tsrc));
            theta=acos(k*kp(psr,:)');
            %sprintf('%d pulsar theta=%g',i,theta)
            %phiI(i)=mod(phi0-omega*distP(i)*(1-cos(theta)), 2*pi);  % modulus after division
            %phiI(i)=mod(2*phi0-omega_tmp(l)*distP(i)*(1-cos(theta)), pi);  % modulus after division, YW 09/10/13
            phiI(psr)=mod(SrcPhi0{band}(id_tsrc)-0.5*SrcOmega{band}(id_tsrc)*distP(psr)*(1-cos(theta)), pi);  % modulus after division, YW 04/30/14 check original def. of phiI
            
            tmp = FullResiduals(SrcAlpha{band}(id_tsrc),SrcDelta{band}(id_tsrc),SrcOmega{band}(id_tsrc),SrcPhi0{band}(id_tsrc),phiI(psr),alphaP(psr),deltaP(psr),...
                SrcAmp{band}(id_tsrc),SrcIota{band}(id_tsrc),SrcThetaN{band}(id_tsrc),theta,yr); % timing residuals for true src
            
            tmp_est = FullResiduals(EstSrc{band,src}.alpha,EstSrc{band,src}.delta,EstSrc{band,src}.omega,EstSrc{band,src}.phi0,EstSrc{band,src}.phiI(psr),alphaP(psr),deltaP(psr),...
                EstSrc{band,src}.Amp,EstSrc{band,src}.iota,EstSrc{band,src}.thetaN,theta,yr); % timing residuals for estimated source
            
            rho_tmp(psr,1) = abs(tmp' * tmp_est/(norm(tmp)*norm(tmp_est))); % cross-correlation
        end
        
        rho_tmp2(src,id_tsrc) = max(rho_tmp);
        
        dif_freq{band}(src,1) = abs(SrcOmega{band}(id_tsrc) - EstSrc{band,src}.omega);
        dif_ra{band}(src,1) = abs(SrcAlpha{band}(id_tsrc) - EstSrc{band,src}.alpha);
        dif_dec{band}(src,1) = abs(SrcDelta{band}(id_tsrc) - EstSrc{band,src}.delta);
    end
    
    rho{band} = rho_tmp2;
    [rho_max{band},id_max(:,band)] = max(rho{band},[],2); % get index of true sources when rho reaches maximum
    %     dif_freq_max(:,band) = abs(arrayfun(@(x) EstSrc{band,x}.omega, 1:NestsrcBand) - SrcOmega{band}(id_max(:,band)')) / (365*24*3600*2*pi); % convert to Hz
    dif_freq_max(:,band) = abs(arrayfun(@(x) EstSrc{band,x}.omega, 1:NestsrcBand) - SrcOmega{band}(id_max(:,band)')) * 100 ./ SrcOmega{band}(id_max(:,band)'); % error as percentage
    %     dif_ra_max(:,band) = abs(arrayfun(@(x) EstSrc{band,x}.alpha, 1:NestsrcBand) - SrcAlpha{band}(id_max(:,band)'));
    dif_ra_max(:,band) = abs(arrayfun(@(x) EstSrc{band,x}.alpha, 1:NestsrcBand) - SrcAlpha{band}(id_max(:,band)')) * 100 ./ SrcAlpha{band}(id_max(:,band)'); % error as percentage
    %     dif_dec_max(:,band) = abs(arrayfun(@(x) EstSrc{band,x}.delta, 1:NestsrcBand) - SrcDelta{band}(id_max(:,band)'));
    dif_dec_max(:,band) = abs(arrayfun(@(x) EstSrc{band,x}.delta, 1:NestsrcBand) - SrcDelta{band}(id_max(:,band)')) * 100 ./ abs(SrcDelta{band}(id_max(:,band)')); % error as percentage
end



% EOF