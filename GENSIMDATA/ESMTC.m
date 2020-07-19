function [rho,rho_max,dif_freq_max,dif_ra_max,dif_dec_max,id_max,estSNR1,estSNR2] = ESMTC(Nband,NestsrcBand,EstSrc1,EstSrc2,simParams,yr,threshold)
% A function calculates cross-correlation coefficients above a chosen
% threshold for different estimated sources.
% [rho,rho_max,dif_freq_max,dif_ra_max,dif_dec_max,estSNR] =
% MWC(Nband,NestsrcBand,EstSrc1,EstSrc2,simParams,yr,threshold)
% rho: cross-correlation coefficient matrix
% rho_max: maximum value of rho.
% dif_freq_max: error in frequency.
% dif_ra_max: error in RA.
% dif_dec_max: error in DEC.
% id_max: index of sources which reach the maximum of coefficients.
% estSNRx: SNR value for estimated sources.
% Nband: number of band.
% NestsrcBand: number of estimated source in each band.
% EstSrcx: cell collects all the parameters of estimated source.
% simParams: pulsar config.
% yr: observation span.
% threshold: a threshold chosen by user.

% Author: QYQ 5/24/2020
%%

% display weight factor used
% if strcmp(weight,'omega') == 1
%     disp("Using Freq. as weight factor")
% elseif strcmp(weight,'snr') == 1
%     disp("Using SNR as weight factor")
% elseif weight == 0
%     disp("Using no weight factor")
% end

%% Cross-Corelation
Np = simParams.Np; % number of pulsars
deltaP = simParams.deltaP;
alphaP = simParams.alphaP;
kp = simParams.kp;
% distP = simParams.distP;

rho_tmp = zeros(Np,1);

rho = {}; % cross-correlation coefficients matrix
% rho_max = {};
gamma = cell(1,2); % averaged cross-correlation coefficients

% investigations
estSNR1 = zeros(Nband,NestsrcBand);
estSNR2 = zeros(Nband,NestsrcBand);
% dif_freq = {}; % frequency difference
% dif_ra = {};
% dif_dec = {};

id_max = zeros(NestsrcBand,Nband); % index of max. cc
dif_freq_max = zeros(NestsrcBand,Nband);
dif_ra_max = zeros(NestsrcBand,Nband);
dif_dec_max = zeros(NestsrcBand,Nband);


for band = 1:Nband
    NestSrc1 = length(EstSrc1); % number of true sources in each band
    gamma{band} = zeros(NestsrcBand,NestSrc1); % initialize the gamma cell.
    for src = 1:NestsrcBand
        [snr1,~] = Amp2Snr(EstSrc1{band,src},simParams,yr); % get SNR for estimated source
        estSNR1(band,src) = snr1;
        
        for tsrc = 1:NestSrc1
            [snr2,~] = Amp2Snr(EstSrc2{band,tsrc},simParams,yr); % get SNR for estimated source
            estSNR2(band,tsrc) = snr2;
            %             tmp_true = 0; % for gamma star
            %             tmp_est1 = 0;
            %             tmp2 = 0;
            for psr = 1:Np
                % GW sky location in Cartesian coordinate
                k=zeros(1,3);  % unit vector pointing from SSB to source
                k(1)=cos(EstSrc1{band,tsrc}.delta)*cos(EstSrc1{band,tsrc}.alpha);
                k(2)=cos(EstSrc1{band,tsrc}.delta)*sin(EstSrc1{band,tsrc}.alpha);
                k(3)=sin(EstSrc1{band,tsrc}.delta);
                theta=acos(k*kp(psr,:)');
                %sprintf('%d pulsar theta=%g',i,theta)
                %phiI(i)=mod(phi0-omega*distP(i)*(1-cos(theta)), 2*pi);  % modulus after division
                %phiI(i)=mod(2*phi0-omega_tmp(l)*distP(i)*(1-cos(theta)), pi);  % modulus after division, YW 09/10/13
                %                 phiI(psr)=mod(EstSrc1{band,tsrc}.phi0-0.5*EstSrc1{band,tsrc}.omega*distP(psr)*(1-cos(theta)), pi);  % modulus after division, YW 04/30/14 check original def. of phiI
                
                tmp = FullResiduals(EstSrc1{band,tsrc}.alpha,EstSrc1{band,tsrc}.delta,EstSrc1{band,tsrc}.omega,EstSrc1{band,tsrc}.phi0,EstSrc1{band,tsrc}.phiI(psr),alphaP(psr),deltaP(psr),...
                    EstSrc1{band,tsrc}.Amp,EstSrc1{band,tsrc}.iota,EstSrc1{band,tsrc}.thetaN,theta,yr); % timing residuals for true src
                
                tmp_est = FullResiduals(EstSrc2{band,src}.alpha,EstSrc2{band,src}.delta,EstSrc2{band,src}.omega,EstSrc2{band,src}.phi0,EstSrc2{band,src}.phiI(psr),alphaP(psr),deltaP(psr),...
                    EstSrc2{band,src}.Amp,EstSrc2{band,src}.iota,EstSrc2{band,src}.thetaN,theta,yr); % timing residuals for estimated source
                
                %                 tmp_est1 = tmp_est1 + norm(tmp_est);
                %                 tmp_true = tmp_true + norm(tmp); % for gamma star
                
                rho_tmp(psr,tsrc) = abs(tmp' * tmp_est/(norm(tmp)*norm(tmp_est))); % cross-correlation
                %                 rho_tmp(psr,tsrc) = tmp' * tmp_est/(norm(tmp)*norm(tmp_est)); % don't use abs.
                
                %                 tmp2 = tmp2 + rho_tmp(psr,tsrc) * norm(tmp_est);
                %                 tmp2 = tmp2 + rho_tmp(psr,tsrc) * norm(tmp); % for gamma star
                
            end
         
        end
        above_threshold = sum(rho_tmp > threshold);
        [~,id_max(src,band)] = max(above_threshold);
        
        %         gamma{band}(src,id_max(src,band)) = max(rho_tmp(:,id_max(src,band))); % maximize method
        gamma{band}(src,id_max(src,band)) = sum(rho_tmp(rho_tmp(:,id_max(src,band)) > threshold,id_max(src,band))) / 1000; % nomalized over Np (1000) pulsars.
    end
    rho{band} = gamma{band};
    rho_max{band} = max(rho{band},[],2); % get index of true sources when rho reaches maximum
    
    
    %     dif_freq_max(:,band) = abs(arrayfun(@(x) EstSrc{band,x}.omega, 1:NestsrcBand) - SrcOmega{band}(id_max(:,band)')) / (365*24*3600*2*pi); % convert to Hz
    dif_freq_max(:,band) = arrayfun(@(x) abs(EstSrc2{band,x}.omega - EstSrc1{band,id_max(x,band)}.omega) * 100 / EstSrc1{band,id_max(x,band)}.omega, 1:NestsrcBand); % error as percentage
    %     dif_ra_max(:,band) = abs(arrayfun(@(x) EstSrc{band,x}.alpha, 1:NestsrcBand) - SrcAlpha{band}(id_max(:,band)'));
    dif_ra_max(:,band) = arrayfun(@(x) abs(EstSrc2{band,x}.alpha - EstSrc1{band,id_max(x,band)}.alpha) * 100 / EstSrc1{band,id_max(x,band)}.alpha, 1:NestsrcBand); % error as percentage
    %     dif_dec_max(:,band) = abs(arrayfun(@(x) EstSrc{band,x}.delta, 1:NestsrcBand) - SrcDelta{band}(id_max(:,band)'));
    dif_dec_max(:,band) = arrayfun(@(x) abs(EstSrc2{band,x}.delta - EstSrc1{band,id_max(x,band)}.delta) * 100 / EstSrc1{band,id_max(x,band)}.delta, 1:NestsrcBand); % error as percentage
end