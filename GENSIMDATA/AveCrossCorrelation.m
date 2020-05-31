% Averaged Cross-Correlation Coefficients Matrix
% for debug use.
% Author: QYQ
% 05/22/2020

clear;
tic

%% Dir settings
simParamsDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/searchParams/2bands/superNarrow';
simdataDir = '/Users/qianyiqian/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands/';
estdataDir = '/Users/qianyiqian/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands/supperNarrow_iMBLT_final';
Filename = 'GWBsimDataSKASrlz1Nrlz3';
ext = '.mat';

%% Files
paraFile = dir([simParamsDir,filesep,'searchParams','*.mat']);
Nband = length(paraFile);
simFile = [simdataDir,filesep,Filename,ext];
estFile = dir([estdataDir,filesep,'*',Filename,'*',ext]);
Nestsrc = length(estFile);
paraFilename = sort_nat({paraFile.name});
estFilename = sort_nat({estFile.name});
load(simFile);

%% Seperate sources into different bands
% Ntsrc = length(alpha); % Number of true sources.
SrcSNR = {};
SrcAlpha = {};
SrcAmp = {};
SrcDelta = {};
SrcIota = {};
SrcOmega = {};
SrcPhi0 = {};
SrcThetaN = {};

for i = 1:Nband
    load([simParamsDir,filesep,char(paraFilename(i))]);
    Indx = find(omega >= searchParams.angular_velocity(2) & ...
        omega <= searchParams.angular_velocity(1));
    
    SrcSNR{i} = snr_chr(Indx);
    SrcAlpha{i} = alpha(Indx);
    SrcDelta{i} = delta(Indx);
    SrcAmp{i} = Amp(Indx);
    SrcIota{i} = iota(Indx);
    SrcOmega{i} = omega(Indx);
    SrcPhi0{i} = phi0(Indx);
    SrcThetaN{i} = thetaN(Indx);
    
end

%% Sort sources in different bands

for j = 1:Nband
    [~,id] = sort(SrcSNR{j},'descend'); % sort true sources according to SNR value
    SrcSNR{j} = SrcSNR{j}(id);
    SrcAlpha{j} = SrcAlpha{j}(id);
    SrcDelta{j} = SrcDelta{j}(id);
    SrcAmp{j} = SrcAmp{j}(id);
    SrcIota{j} = SrcIota{j}(id);
    SrcOmega{j} = SrcOmega{j}(id);
    SrcPhi0{j} = SrcPhi0{j}(id);
    SrcThetaN{j} = SrcThetaN{j}(id);
end


%% Get estimated sources info
NestsrcBand = Nestsrc/Nband; % number of sources in a band.
EstSrc = {};
for band = 1:Nband
    for k = 1:NestsrcBand
        path_to_estimatedData = [estdataDir,filesep,char(estFilename((band - 1) * NestsrcBand + k))];
        EstSrc{band,k} = ColSrcParams(path_to_estimatedData);
    end
end


%% Cross-Corelation
Np = simParams.Np; % number of pulsars
deltaP = simParams.deltaP;
alphaP = simParams.alphaP;
kp = simParams.kp;
distP = simParams.distP;

rho_tmp = zeros(Np,1);

rho = {}; % cross-correlation coefficients matrix
rho_max = {};
gamma = {}; % averaged cross-correlation coefficients

% investigations
estSNR = zeros(Nband,NestsrcBand);
dif_freq = {}; % frequency difference
dif_ra = {};
dif_dec = {};

id_max = zeros(NestsrcBand,Nband); % index of max. cc
dif_freq_max = zeros(NestsrcBand,Nband);
dif_ra_max = zeros(NestsrcBand,Nband);
dif_dec_max = zeros(NestsrcBand,Nband);

tmp3 = zeros(Np,1);
tmp4 = zeros(Np,1);
for band = 1:Nband
    NtsrcBand = length(SrcAlpha{band}); % number of true sources in each band

    for src = 1:1%NestsrcBand
        [snr,~] = Amp2Snr(EstSrc{band,src},simParams,yr); % get SNR for estimated source
        estSNR(band,src) = snr;
        for tsrc = 1:NtsrcBand
            tmp_est1 = 0;
            tmp2 = 0;
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
                tmp_est1 = tmp_est1 + norm(tmp_est);
                
                rho_tmp(psr,tsrc) = abs(tmp' * tmp_est/(norm(tmp)*norm(tmp_est))); % cross-correlation
                tmp2 = tmp2 + rho_tmp(psr,tsrc) * norm(tmp_est);
                tmp3(psr,tsrc) = rho_tmp(psr,tsrc) * norm(tmp_est);
                tmp4(psr,tsrc) = norm(tmp);
                
            end
            %             w = 1/(1+abs(SrcSNR{band}(tsrc)-snr)/SrcSNR{band}(tsrc)); % weight factor: SNR
            w = 1/(1+abs(SrcOmega{band}(tsrc)-EstSrc{band,src}.omega) / SrcOmega{band}(tsrc)); % weight factor: Freq
            %             w = 1; % no weight used
            %             rho_tmp2(src,tsrc) = max(rho_tmp) * 1/(1+abs(SrcSNR{band}(tsrc)-snr)/SrcSNR{band}(tsrc)); % add a weight function using SNR
            gamma{band}(src,tsrc) = tmp2 * w / tmp_est1; % averaged c-c coefficients
            
            dif_freq{band}(src,tsrc) = abs(SrcOmega{band}(tsrc) - EstSrc{band,src}.omega);
            dif_ra{band}(src,tsrc) = abs(SrcAlpha{band}(tsrc) - EstSrc{band,src}.alpha);
            dif_dec{band}(src,tsrc) = abs(SrcDelta{band}(tsrc) - EstSrc{band,src}.delta);
        end
        
    end
    rho{band} = gamma{band};
    [rho_max{band},id_max(:,band)] = max(rho{band},[],2); % get index of true sources when rho reaches maximum
    %     dif_freq_max(:,band) = abs(arrayfun(@(x) EstSrc{band,x}.omega, 1:NestsrcBand) - SrcOmega{band}(id_max(:,band)')) / (365*24*3600*2*pi); % convert to Hz
    dif_freq_max(:,band) = abs(arrayfun(@(x) EstSrc{band,x}.omega, 1:NestsrcBand) - SrcOmega{band}(id_max(:,band)')) * 100 ./ SrcOmega{band}(id_max(:,band)'); % error as percentage
    %     dif_ra_max(:,band) = abs(arrayfun(@(x) EstSrc{band,x}.alpha, 1:NestsrcBand) - SrcAlpha{band}(id_max(:,band)'));
    dif_ra_max(:,band) = abs(arrayfun(@(x) EstSrc{band,x}.alpha, 1:NestsrcBand) - SrcAlpha{band}(id_max(:,band)')) * 100 ./ SrcAlpha{band}(id_max(:,band)'); % error as percentage
    %     dif_dec_max(:,band) = abs(arrayfun(@(x) EstSrc{band,x}.delta, 1:NestsrcBand) - SrcDelta{band}(id_max(:,band)'));
    dif_dec_max(:,band) = abs(arrayfun(@(x) EstSrc{band,x}.delta, 1:NestsrcBand) - SrcDelta{band}(id_max(:,band)')) * 100 ./ abs(SrcDelta{band}(id_max(:,band)')); % error as percentage
end

% set coefficients of sources which are not maximum to zero
% for band = 1:Nband
%     NtsrcBand = length(SrcAlpha{band});
%     for src = 1:NestsrcBand
%         for tsrc = 1:NtsrcBand
%             if tsrc ~= id_max(src,band)
%                 rho{band}(src,tsrc) = 0;
%             end
%         end
%     end
% end



%% Plotting
prefix = [estdataDir,filesep,'fig',filesep,Filename];
mkdir(prefix);
figname = 'Ave-Cross-Corelation';

for fig = 1:Nband
    figure
    imagesc(rho{fig});
    colorbar
    xlabel('True sources')
    ylabel('Estimated sources')
    title(['Band ',num2str(fig)])
    saveas(gcf,[prefix,filesep,figname,'Band ',num2str(fig)],'png');
    savefig([prefix,filesep,figname,'Band ',num2str(fig)]);
end

figname2 = 'AveCC_estSNR';
for fig2 = 1:Nband
    figure
    plot(estSNR(fig2,:),rho_max{fig2},'ob')
    xlabel('Estimated SNR')
    ylabel('Ave. Weighted CC-freq')
    title(['Band ',num2str(fig2)])
    saveas(gcf,[prefix,filesep,figname2,'Band ',num2str(fig2)],'png');
    savefig([prefix,filesep,figname2,'Band ',num2str(fig2)]);
end

% figname3 = 'Freq_AveCC';
% for fig3 = 1:Nband
%     for n = 1:NestsrcBand
%         figure
%         plot(dif_freq{fig3}(n,:),rho{fig3}(n,:),'ob')
%         xlabel('Freq. difference')
%         ylabel('Weighted cross-correlation')
%         title(['Estimated source ',num2str(n)])
%         saveas(gcf,[prefix,filesep,figname3,'EstSrc ',num2str(n)],'png');
%         savefig([prefix,filesep,figname3,'EstSrc ',num2str(n)]);
%     end
% end
% 
% figname4 = 'RA_AveCC';
% for fig3 = 1:Nband
%     for n = 1:NestsrcBand
%         figure
%         plot(dif_ra{fig3}(n,:),rho{fig3}(n,:),'ob')
%         xlabel('RA difference')
%         ylabel('Weighted cross-correlation')
%         title(['Estimated source ',num2str(n)])
%         saveas(gcf,[prefix,filesep,figname4,'EstSrc ',num2str(n)],'png');
%         savefig([prefix,filesep,figname4,'EstSrc ',num2str(n)]);
%     end
% end
% 
% 
% figname5 = 'DEC_AveCC';
% for fig3 = 1:Nband
%     for n = 1:NestsrcBand
%         figure
%         plot(dif_dec{fig3}(n,:),rho{fig3}(n,:),'ob')
%         xlabel('DEC difference')
%         ylabel('Weighted cross-correlation')
%         title(['Estimated source ',num2str(n)])
%         saveas(gcf,[prefix,filesep,figname5,'EstSrc ',num2str(n)],'png');
%         savefig([prefix,filesep,figname5,'EstSrc ',num2str(n)]);
%     end
% end


figname6 = 'Ave-CC-Freq';

for fig = 1:Nband
    figure
    plot(dif_freq_max(:,fig),rho_max{fig},'ob')
    xlabel('Difference of Freq. Percentage (%)')
    ylabel('Ave. Weighted CC-freq')
    title(['Band ',num2str(fig)])
    saveas(gcf,[prefix,filesep,figname6,'Band ',num2str(fig)],'png');
    savefig([prefix,filesep,figname6,'Band ',num2str(fig)]);
end

figname7 = 'Ave-CC-RA';

for fig = 1:Nband
    figure
    plot(dif_ra_max(:,fig),rho_max{fig},'ob')
    xlabel('Difference of RA Percentage (%)')
    ylabel('Ave. Weighted CC-freq')
    title(['Band ',num2str(fig)])
    saveas(gcf,[prefix,filesep,figname7,'Band ',num2str(fig)],'png');
    savefig([prefix,filesep,figname7,'Band ',num2str(fig)]);
end

figname8 = 'Ave-CC-DEC';

for fig = 1:Nband
    figure
    plot(dif_dec_max(:,fig),rho_max{fig},'ob')
    xlabel('Difference of DEC Percentage (%)')
    ylabel('Ave. Weighted CC-freq')
    title(['Band ',num2str(fig)])
    saveas(gcf,[prefix,filesep,figname8,'Band ',num2str(fig)],'png');
    savefig([prefix,filesep,figname8,'Band ',num2str(fig)]);
end




toc