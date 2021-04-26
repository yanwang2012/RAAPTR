function [gamma,rho,rho_max,id_max] = TNMTC(Nband,simSrc,simParams,yr,threshold)
% A function calculates cross-correlation coefficients above a chosen
% threshold for true sources vs. true sources
% [rho,rho_max,dif_freq_max,dif_ra_max,dif_dec_max,estSNR] =
% TNMTC(Nband,simSrc,simParams,yr,threshold)
% rho: cross-correlation coefficient matrix
% rho_max: maximum value of rho.
% dif_freq_max: error in frequency.
% dif_ra_max: error in RA.
% dif_dec_max: error in DEC.
% id_max: index of sources which reach the maximum of coefficients.
% Nband: number of band.
% simSrc: simulated sources.
% simParams: pulsar config.
% yr: observation span.
% threshold: a threshold chosen by user.

% Author: QYQ 04/26/2021
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
distP = simParams.distP;

SrcAlpha = {simSrc.SrcAlpha};
SrcDelta = {simSrc.SrcDelta};
SrcOmega = {simSrc.SrcOmega};
SrcPhi0 = {simSrc.SrcPhi0};
SrcIota = {simSrc.SrcIota};
SrcThetaN = {simSrc.SrcThetaN};
SrcAmp = {simSrc.SrcAmp};

rho_tmp = zeros(Np,1);

rho = {}; % cross-correlation coefficients matrix
% rho_max = {};
gamma = cell(1,Nband); % averaged cross-correlation coefficients

% investigations
% NestsrcBand = max(EstSrcBand);
% estSNR = zeros(Nband,NestsrcBand);
% % dif_freq = {}; % frequency difference
% % dif_ra = {};
% % dif_dec = {};
% 
% id_max = zeros(NestsrcBand,Nband); % index of max. cc
% dif_freq_max = zeros(NestsrcBand,Nband);
% dif_ra_max = zeros(NestsrcBand,Nband);
% dif_dec_max = zeros(NestsrcBand,Nband);


for band = 1:Nband
    NtsrcBand = length(SrcAlpha{band}); % number of true sources in each band
%     NestsrcBand = EstSrcBand(band);
    gamma{band} = zeros(NtsrcBand,NtsrcBand); % initialize the gamma cell.
    % search along x-axis
    for src1 = 1:NtsrcBand
%         [snr,~] = Amp2Snr(EstSrc{band,src1},simParams,yr); % get SNR for estimated source
%         estSNR(band,src1) = snr;
        for src2 = 1:NtsrcBand
            %             tmp_true = 0; % for gamma star
            %             tmp_est1 = 0;
            %             tmp2 = 0;
            for psr = 1:Np
                % GW sky location in Cartesian coordinate
                k=zeros(1,3);  % unit vector pointing from SSB to source
                k(1)=cos(SrcDelta{band}(src2))*cos(SrcAlpha{band}(src2));
                k(2)=cos(SrcDelta{band}(src2))*sin(SrcAlpha{band}(src2));
                k(3)=sin(SrcDelta{band}(src2));
                theta=acos(k*kp(psr,:)');
                %sprintf('%d pulsar theta=%g',i,theta)
                %phiI(i)=mod(phi0-omega*distP(i)*(1-cos(theta)), 2*pi);  % modulus after division
                %phiI(i)=mod(2*phi0-omega_tmp(l)*distP(i)*(1-cos(theta)), pi);  % modulus after division, YW 09/10/13
                phiI(psr)=mod(SrcPhi0{band}(src2)-0.5*SrcOmega{band}(src2)*distP(psr)*(1-cos(theta)), pi);  % modulus after division, YW 04/30/14 check original def. of phiI
                
                tmp = FullResiduals(SrcAlpha{band}(src2),SrcDelta{band}(src2),SrcOmega{band}(src2),SrcPhi0{band}(src2),phiI(psr),alphaP(psr),deltaP(psr),...
                    SrcAmp{band}(src2),SrcIota{band}(src2),SrcThetaN{band}(src2),theta,yr); % timing residuals for true src
                
                tmp_src1 = FullResiduals(SrcAlpha{band}(src1),SrcDelta{band}(src1),SrcOmega{band}(src1),SrcPhi0{band}(src1),phiI(psr),alphaP(psr),deltaP(psr),...
                    SrcAmp{band}(src1),SrcIota{band}(src1),SrcThetaN{band}(src1),theta,yr); % timing residuals for estimated source
                
                %                 tmp_est1 = tmp_est1 + norm(tmp_est);
                %                 tmp_true = tmp_true + norm(tmp); % for gamma star
                
                rho_tmp(psr,src2) = abs(tmp' * tmp_src1/(norm(tmp)*norm(tmp_src1))); % cross-correlation
                
                %                 tmp2 = tmp2 + rho_tmp(psr,tsrc) * norm(tmp_est);
                %                 tmp2 = tmp2 + rho_tmp(psr,tsrc) * norm(tmp); % for gamma star
                
            end
            % choose weight factor
            %             if strcmp(weight,'omega') == 1
            %                 %                 disp("Using Freq. as weight factor")
            %                 w = 1/(1+abs(SrcOmega{band}(tsrc)-EstSrc{band,src}.omega) / SrcOmega{band}(tsrc)); % weight factor: Freq
            %             elseif strcmp(weight,'snr') == 1
            %                 %                 disp("Using SNR as weight factor")
            %                 w = 1/(1+abs(SrcSNR{band}(tsrc)-snr)/SrcSNR{band}(tsrc)); % weight factor: SNR
            %             elseif weight == 0
            %                 %                 disp("Using no weight factor")
            %                 w = 1; % no weight used
            %             end
            
            %             gamma{band}(src,tsrc) = tmp2 * w / tmp_true; % averaged c-c coefficients using gamma star.
            %             gamma{band}(src,tsrc) = tmp2 * w / tmp_est1; % averaged c-c coefficients
            
            %             dif_freq{band}(src,tsrc) = abs(SrcOmega{band}(tsrc) - EstSrc{band,src}.omega);
            %             dif_ra{band}(src,tsrc) = abs(SrcAlpha{band}(tsrc) - EstSrc{band,src}.alpha);
            %             dif_dec{band}(src,tsrc) = abs(SrcDelta{band}(tsrc) - EstSrc{band,src}.delta);
        end
        above_threshold = sum(rho_tmp > threshold); % calculate how many CC. above the threshold.
        [~,id_max(src1,band)] = max(above_threshold);
        gamma{band}(src1,:) = sum(rho_tmp)/Np; % cc values of src1 for each src2
        %         gamma{band}(src,id_max(src,band)) = max(rho_tmp(:,id_max(src,band))); % Maximized CC
        rho{band}(src1,id_max(src1,band)) = sum(rho_tmp(:,id_max(src1,band))) / Np; % nomalized over Np (1000) pulsars.
        %         gamma{band}(src,id_max(src,band)) = sum(rho_tmp(:,id_max(src,band)) > threshold) / 1000;
    end
    
    rho_tmp = zeros(Np,1); % needs to re-init. rho_tmp when change band, since size of rho_tmp is not constant.
    rho_max{band} = max(rho{band},[],2); % get the maximum CC. value over true sources.
    % make id_max unique
    %     for src1 = 1:NestsrcBand - 1
    %         for src2 = src1+1:NestsrcBand
    %             if id_max(src1,band) == id_max(src2,band)
    %                 [~,I] = sort(rho{band}(src2,:),'descend');
    %                 id_max(src2,band) = I(2); % take the index of second maximum.
    %                 rho_max{band}(src2) = rho{band}(src2,id_max(src2,band));
    %             end
    %         end
    %     end
    
    %     dif_freq_max(:,band) = abs(arrayfun(@(x) EstSrc{band,x}.omega, 1:NestsrcBand) - SrcOmega{band}(id_max(:,band)')) / (365*24*3600*2*pi); % convert to Hz
%     dif_freq_max(1:NestsrcBand,band) = arrayfun(@(x) abs(EstSrc{band,x}.omega - SrcOmega{band}(id_max(x,band))) * 100 / SrcOmega{band}(id_max(x,band)), 1:NestsrcBand); % error as percetage
%     dif_ra_max(1:NestsrcBand,band) = arrayfun(@(x) abs(EstSrc{band,x}.alpha - SrcAlpha{band}(id_max(x,band))) * 100 / SrcAlpha{band}(id_max(x,band)), 1:NestsrcBand); % error as percentage
%     %     dif_dec_max(:,band) = abs(arrayfun(@(x) EstSrc{band,x}.delta, 1:NestsrcBand) - SrcDelta{band}(id_max(:,band)'));
%     dif_dec_max(1:NestsrcBand,band) = arrayfun(@(x) abs(EstSrc{band,x}.delta - SrcDelta{band}(id_max(x,band))) * 100 / SrcDelta{band}(id_max(x,band)), 1:NestsrcBand); % error as percentage
end