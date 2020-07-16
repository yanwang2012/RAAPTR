function [rho,rho_max,dif_freq_max,dif_ra_max,dif_dec_max,id_max,estSNR] = MWC(Nband,NestsrcBand,SrcAlpha,SrcDelta,SrcOmega,SrcPhi0,...
    SrcIota,SrcThetaN,SrcAmp,SrcSNR,EstSrc,simParams,yr,weight)
% A function calculates cross-correlation coefficients using max weighted
% algorithm.
% [rho,rho_max,dif_freq_max,dif_ra_max,dif_dec_max,id_max] =
% MWC(Nband,NestsrcBand,SrcAlpha,SrcDelta,SrcOmega,SrcPhi0,SrcIota,SrcThetaN,SrcAmp,SrcSNR,EstSrc,simParams,yr,weight)
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
% weight: choose weight factor, can be 'omega', 'snr' and 0.

% Author: QYQ 5/24/2020
%%

% display weight factor used
if strcmp(weight,'omega') == 1
    disp("Using Freq. as weight factor")
elseif strcmp(weight,'snr') == 1
    disp("Using SNR as weight factor")
elseif weight == 0
    disp("Using no weight factor")
end


% Pulsar configuration
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

id_max = zeros(NestsrcBand,Nband); % index of max. cc
dif_freq_max = zeros(NestsrcBand,Nband);
dif_ra_max = zeros(NestsrcBand,Nband);
dif_dec_max = zeros(NestsrcBand,Nband);

for band = 1:Nband
    NtsrcBand = length(SrcAlpha{band}); % number of true sources in each band
    rho_tmp2 = zeros(NestsrcBand,NtsrcBand);
    rho_tmp3 = zeros(NestsrcBand,NtsrcBand);
    for src = 1:NestsrcBand
        [snr,~] = Amp2Snr(EstSrc{band,src},simParams,yr); % get SNR for estimated source
        estSNR(band,src) = snr;
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
                
                tmp_est = FullResiduals(EstSrc{band,src}.alpha,EstSrc{band,src}.delta,EstSrc{band,src}.omega,EstSrc{band,src}.phi0,EstSrc{band,src}.phiI(psr),alphaP(psr),deltaP(psr),...
                    EstSrc{band,src}.Amp,EstSrc{band,src}.iota,EstSrc{band,src}.thetaN,theta,yr); % timing residuals for estimated source
                
                rho_tmp(psr,tsrc) = abs(tmp' * tmp_est/(norm(tmp)*norm(tmp_est))); % cross-correlation
                %                 rho_tmp(psr,1) = tmp' * tmp_est/(norm(tmp)*norm(tmp_est));
                
            end
            % choose weight factor
            if strcmp(weight,'omega') == 1
%                 disp("Using Freq. as weight factor")
                w = 1/(1+abs(SrcOmega{band}(tsrc)-EstSrc{band,src}.omega) / SrcOmega{band}(tsrc)); % weight factor: Freq
            elseif strcmp(weight,'snr') == 1
%                 disp("Using SNR as weight factor")
                w = 1/(1+abs(SrcSNR{band}(tsrc)-snr)/SrcSNR{band}(tsrc)); % weight factor: SNR
            elseif weight == 0
%                 disp("Using no weight factor")
                w = 1; % no weight used
            end
            
            rho_tmp2(src,tsrc) = max(rho_tmp(:,tsrc)) * w;
            
            dif_freq{band}(src,tsrc) = abs(SrcOmega{band}(tsrc) - EstSrc{band,src}.omega);
            dif_ra{band}(src,tsrc) = abs(SrcAlpha{band}(tsrc) - EstSrc{band,src}.alpha);
            dif_dec{band}(src,tsrc) = abs(SrcDelta{band}(tsrc) - EstSrc{band,src}.delta);
        end
        
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