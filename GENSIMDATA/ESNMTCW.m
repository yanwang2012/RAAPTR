function [gamma,rho,id_max,estSNR1,estSNR2] = ESNMTCW(EstSrc1,EstSrc2,simParams,yr,threshold)
% A function calculates cross-correlation coefficients using NMTC above a chosen
% threshold for different set of est. sources without spliting into
% different bands.
% [gamma,rho,dif_freq_max,dif_ra_max,dif_dec_max,id_max,estSNR1,estSNR2] =
% ESNMTCW(Nband,NestsrcBand,EstSrc1,EstSrc2,simParams,yr,threshold)
% gamma: cross-correlation coefficient matrix.
% rho: max(gamma).
% dif_freq_max: error in frequency.
% dif_ra_max: error in RA.
% dif_dec_max: error in DEC.
% id_max: index of sources which reach the maximum of coefficients.
% estSNRx: SNR value for estimated sources.
% Nestsrc: number of est. sources in each set of data realization.
% EstSrcx: cell collects all the parameters of estimated source.
% simParams: pulsar config.
% yr: observation span.
% threshold: a threshold chosen by user.
% output rho in the shape of [estSrc1,estSrc2].

% Author: QYQ 09/14/2020

%% Cross-Corelation
Np = simParams.Np; % number of pulsars
deltaP = simParams.deltaP;
alphaP = simParams.alphaP;
kp = simParams.kp;
% distP = simParams.distP;

rho_tmp = zeros(Np,1);

Nestsrc1 = length(EstSrc1);
Nestsrc2 = length(EstSrc2);

rho = zeros(Nestsrc1,Nestsrc2); % cross-correlation coefficients matrix
gamma = zeros(Nestsrc1,Nestsrc2); % averaged cross-correlation coefficients

%% investigations
estSNR1 = zeros(Nestsrc1,1);
estSNR2 = zeros(Nestsrc2,1);
% dif_freq = {}; % frequency difference
% dif_ra = {};
% dif_dec = {};

id_max = zeros(Nestsrc1,1); % index of max. cc
% dif_freq_max = zeros(Nestsrc);
% dif_ra_max = zeros(Nestsrc);
% dif_dec_max = zeros(Nestsrc);


%%
for src1 = 1:Nestsrc1
    [snr1,~] = Amp2Snr(EstSrc1{src1},simParams,yr); % get SNR for estimated source
    estSNR1(src1) = snr1;
    
    for src2 = 1:Nestsrc2
        [snr2,~] = Amp2Snr(EstSrc2{src2},simParams,yr); % get SNR for estimated source
        estSNR2(src2) = snr2;
        %             tmp_true = 0; % for gamma star
        %             tmp_est1 = 0;
        %             tmp2 = 0;
        for psr = 1:Np
            % GW sky location in Cartesian coordinate
            k=zeros(1,3);  % unit vector pointing from SSB to source
            k(1)=cos(EstSrc1{src1}.delta)*cos(EstSrc1{src1}.alpha);
            k(2)=cos(EstSrc1{src1}.delta)*sin(EstSrc1{src1}.alpha);
            k(3)=sin(EstSrc1{src1}.delta);
            theta=acos(k*kp(psr,:)');
            %sprintf('%d pulsar theta=%g',i,theta)
            %phiI(i)=mod(phi0-omega*distP(i)*(1-cos(theta)), 2*pi);  % modulus after division
            %phiI(i)=mod(2*phi0-omega_tmp(l)*distP(i)*(1-cos(theta)), pi);  % modulus after division, YW 09/10/13
            %                 phiI(psr)=mod(EstSrc1{band,tsrc}.phi0-0.5*EstSrc1{band,tsrc}.omega*distP(psr)*(1-cos(theta)), pi);  % modulus after division, YW 04/30/14 check original def. of phiI
            
            tmp = FullResiduals(EstSrc1{src1}.alpha,EstSrc1{src1}.delta,EstSrc1{src1}.omega,EstSrc1{src1}.phi0,EstSrc1{src1}.phiI(psr),alphaP(psr),deltaP(psr),...
                EstSrc1{src1}.Amp,EstSrc1{src1}.iota,EstSrc1{src1}.thetaN,theta,yr); % timing residuals for true src
            
            tmp_est = FullResiduals(EstSrc2{src2}.alpha,EstSrc2{src2}.delta,EstSrc2{src2}.omega,EstSrc2{src2}.phi0,EstSrc2{src2}.phiI(psr),alphaP(psr),deltaP(psr),...
                EstSrc2{src2}.Amp,EstSrc2{src2}.iota,EstSrc2{src2}.thetaN,theta,yr); % timing residuals for estimated source
            
            %                 tmp_est1 = tmp_est1 + norm(tmp_est);
            %                 tmp_true = tmp_true + norm(tmp); % for gamma star
            
            rho_tmp(psr,src2) = abs(tmp' * tmp_est/(norm(tmp)*norm(tmp_est))); % cross-correlation
            %                 rho_tmp(psr,tsrc) = tmp' * tmp_est/(norm(tmp)*norm(tmp_est)); % don't use abs.
            
            %                 tmp2 = tmp2 + rho_tmp(psr,tsrc) * norm(tmp_est);
            %                 tmp2 = tmp2 + rho_tmp(psr,tsrc) * norm(tmp); % for gamma star
            
        end
        
    end
    above_threshold = sum(rho_tmp > threshold);
    [val,id_max(src1)] = max(above_threshold);
    threshold_tmp = threshold;
        % fix the bug when max(above_threshold) is 0, it automatically maps
        % estimated source to the first true source.
        while val == 0
            threshold_tmp = threshold_tmp - 0.1;
            above_threshold = sum(rho_tmp > threshold_tmp);
            [val,id_max(src1)] = max(above_threshold);
        end
    
    %         gamma{band}(src,id_max(src,band)) = max(rho_tmp(:,id_max(src,band))); % maximize method
    gamma(src1,id_max(src1)) = sum(rho_tmp(:,id_max(src1))) / Np; % nomalized over Np (1000) pulsars.
end
rho = max(gamma,[],2); % get index of true sources when rho reaches maximum


% %     dif_freq_max(:,band) = abs(arrayfun(@(x) EstSrc{band,x}.omega, 1:NestsrcBand) - SrcOmega{band}(id_max(:,band)')) / (365*24*3600*2*pi); % convert to Hz
% dif_freq_max(:,band) = arrayfun(@(x) abs(EstSrc2{band,x}.omega - EstSrc1{band,id_max(x,band)}.omega) * 100 / EstSrc1{band,id_max(x,band)}.omega, 1:NestsrcBand); % error as percentage
% %     dif_ra_max(:,band) = abs(arrayfun(@(x) EstSrc{band,x}.alpha, 1:NestsrcBand) - SrcAlpha{band}(id_max(:,band)'));
% dif_ra_max(:,band) = arrayfun(@(x) abs(EstSrc2{band,x}.alpha - EstSrc1{band,id_max(x,band)}.alpha) * 100 / EstSrc1{band,id_max(x,band)}.alpha, 1:NestsrcBand); % error as percentage
% %     dif_dec_max(:,band) = abs(arrayfun(@(x) EstSrc{band,x}.delta, 1:NestsrcBand) - SrcDelta{band}(id_max(:,band)'));
% dif_dec_max(:,band) = arrayfun(@(x) abs(EstSrc2{band,x}.delta - EstSrc1{band,id_max(x,band)}.delta) * 100 / EstSrc1{band,id_max(x,band)}.delta, 1:NestsrcBand); % error as percentage