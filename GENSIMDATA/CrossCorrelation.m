% Cross-Correlation Coefficients Matrix
% Author: QYQ
% 05/13/2020

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

% Max Weighted CC
% [rho,rho_max,dif_freq_max,dif_ra_max,dif_dec_max,id_max,estSNR] = MAC(Nband,NestsrcBand,SrcAlpha,SrcDelta,SrcOmega,SrcPhi0,SrcIota,SrcThetaN,SrcAmp,SrcSNR,EstSrc,simParams,yr,'snr');

% Max Weighted Ave. CC
[rho,rho_max,dif_freq_max,dif_ra_max,dif_dec_max,id_max,estSNR] = MWAC(Nband,NestsrcBand,SrcAlpha,SrcDelta,SrcOmega,SrcPhi0,SrcIota,SrcThetaN,SrcAmp,SrcSNR,EstSrc,simParams,yr,0);



% Max coefficients of sources
for band = 1:Nband
    NtsrcBand = length(SrcAlpha{band});
    for src = 1:NestsrcBand
        for tsrc = 1:NtsrcBand
            if tsrc ~= id_max(src,band)
                rho{band}(src,tsrc) = 0;
            end
        end
    end
end



%% Plotting
prefix = [estdataDir,filesep,'fig',filesep,Filename];
mkdir(prefix);
figname = 'MWAC_0_unique';

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


% figname2 = 'MWAC_SNR_snr';
% for fig2 = 1:Nband
%     figure
%     plot(estSNR(fig2,:),rho_max{fig2},'ob')
%     xlabel('Estimated SNR')
%     ylabel('MWAC_{snr}')
%     title(['Band ',num2str(fig2)])
%     saveas(gcf,[prefix,filesep,figname2,'Band ',num2str(fig2)],'png');
%     savefig([prefix,filesep,figname2,'Band ',num2str(fig2)]);
% end

% figname3 = 'Freq_CC';
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
% figname4 = 'RA_CC';
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
% figname5 = 'DEC_CC';
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


% figname6 = 'MWAC_freq_snr';
% 
% for fig = 1:Nband
%     figure
%     plot(dif_freq_max(:,fig),rho_max{fig},'ob')
%     xlabel('Difference of Freq. Percentage (%)')
%     ylabel('MWAC_{snr}')
%     title(['Band ',num2str(fig)])
%     saveas(gcf,[prefix,filesep,figname6,'Band ',num2str(fig)],'png');
%     savefig([prefix,filesep,figname6,'Band ',num2str(fig)]);
% end
% 
% figname7 = 'MWAC_RA_snr';
% 
% for fig = 1:Nband
%     figure
%     plot(dif_ra_max(:,fig),rho_max{fig},'ob')
%     xlabel('Difference of RA Percentage (%)')
%     ylabel('MWAC_{snr}')
%     title(['Band ',num2str(fig)])
%     saveas(gcf,[prefix,filesep,figname7,'Band ',num2str(fig)],'png');
%     savefig([prefix,filesep,figname7,'Band ',num2str(fig)]);
% end
% 
% figname8 = 'MWAC_DEC_snr';
% 
% for fig = 1:Nband
%     figure
%     plot(dif_dec_max(:,fig),rho_max{fig},'ob')
%     xlabel('Difference of DEC Percentage (%)')
%     ylabel('MWAC_{snr}')
%     title(['Band ',num2str(fig)])
%     saveas(gcf,[prefix,filesep,figname8,'Band ',num2str(fig)],'png');
%     savefig([prefix,filesep,figname8,'Band ',num2str(fig)]);
% end

%% plot est. src vs. true src.
figname9 = 'Freq-LogSNR-W0';
cst = 1/(365*24*3600*2*pi); % constant to convert angular frequency to Hz.
EstFreq = zeros(Nband,NestsrcBand);
EstRA = zeros(Nband,NestsrcBand);
EstDEC = zeros(Nband,NestsrcBand);

for fig = 1:Nband
    EstFreq(fig,:) = arrayfun(@(r) EstSrc{fig,r}.omega * cst, 1:NestsrcBand); % get the freq. of est. src.
    figure
    plot(log(SrcSNR{fig}(id_max(:,fig))),SrcOmega{fig}(id_max(:,fig))*cst,'ob',log(estSNR(fig,:)),EstFreq(fig,:),'sr')
    xlabel('Log(SNR)')
    ylabel('Frequency')
    title(['Band ',num2str(fig)])
    legend('True Sources','Estimated Sources')
    saveas(gcf,[prefix,filesep,figname9,'Band ',num2str(fig)],'png');
    savefig([prefix,filesep,figname9,'Band ',num2str(fig)]);
end

figname10 = 'skyLoc-W0';

for fig = 1:Nband
    EstRA(fig,:) = arrayfun(@(r) EstSrc{fig,r}.alpha, 1:NestsrcBand); % get the RA of est. src.
    EstDEC(fig,:) = arrayfun(@(r) EstSrc{fig,r}.delta, 1:NestsrcBand); % get the DEC of est. src.
    figure
    plot(SrcAlpha{fig}(id_max(:,fig)),SrcDelta{fig}(id_max(:,fig)),'ob',EstRA(fig,:),EstDEC(fig,:),'sr')
    xlabel('RA')
    ylabel('DEC')
    title(['Band ',num2str(fig)])
    legend('True Sources','Estimated Sources')
    saveas(gcf,[prefix,filesep,figname10,'Band ',num2str(fig)],'png');
    savefig([prefix,filesep,figname10,'Band ',num2str(fig)]);
end


toc