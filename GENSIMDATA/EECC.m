% Cross-Correlation Coefficients Matrix for Est. & Est. sources.

% Author: QYQ
% 05/13/2020

clear;
tic

%% Dir settings
simParamsDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/searchParams/2bands/superNarrow';
simdataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands';
estSrc1Dir = '/Users/qyq/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands/SuperNarrow/Results_supNar';
estsrc1 = 'initial';
estSrc2Dir = '/Users/qyq/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands/SuperNarrow/SupNar_xMBLT_iMBLT20/iMBLT20'; 
estsrc2 = 'xMBLT-iMBLT';
Filename = 'GWBsimDataSKASrlz1Nrlz3';
ext = '.mat';

%% Files
paraFile = dir([simParamsDir,filesep,'searchParams*',ext]);
simFile = [simdataDir,filesep,Filename,ext];
estSrc1File = dir([estSrc1Dir,filesep,'*',Filename,'*',ext]);
estSrc2File = dir([estSrc2Dir,filesep,'*',Filename,'*',ext]);
Nestsrc = length(estSrc2File);

paraFilename = sort_nat({paraFile.name});
exp = 'searchParams\d.mat'; % regular expressions for desire file names
paraFilename = regexp(paraFilename,exp,'match');
paraFilename = paraFilename(~cellfun(@isempty,paraFilename)); % get rid of empty cells
Nband = length(paraFilename);

estSrc2Filename = sort_nat({estSrc2File.name});
estSrc1Filename = sort_nat({estSrc1File.name});
load(simFile);

%% pre-process true sources
% Seperate sources into different bands
Ntsrc = length(alpha); % Number of true sources.
SrcSNR = {};
SrcAlpha = {};
SrcAmp = {};
SrcDelta = {};
SrcIota = {};
SrcOmega = {};
SrcPhi0 = {};
SrcThetaN = {};

for i = 1:Nband
    load([simParamsDir,filesep,char(paraFilename{i})]);
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

% Sort sources in different bands
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
EstSrc2 = {};
EstSrc1 = {};
for band = 1:Nband
    for k = 1:NestsrcBand
        path_to_estimatedDataestSrc2 = [estSrc2Dir,filesep,char(estSrc2Filename((band - 1) * NestsrcBand + k))];
        path_to_estimatedDataestSrc1 = [estSrc1Dir,filesep,char(estSrc1Filename((band - 1) * NestsrcBand + k))];
        EstSrc2{band,k} = ColSrcParams(path_to_estimatedDataestSrc2);
        EstSrc1{band,k} = ColSrcParams(path_to_estimatedDataestSrc1);
    end
end


%% Cross-Corelation

% Max Weighted CC
% [rho,rho_max,dif_freq_max,dif_ra_max,dif_dec_max,id_max,estSNR] = MAC(Nband,NestsrcBand,SrcAlpha,SrcDelta,SrcOmega,SrcPhi0,SrcIota,SrcThetaN,SrcAmp,SrcSNR,EstSrc,simParams,yr,'snr');

% Max Weighted Ave. CC
% [rho,rho_max,dif_freq_max,dif_ra_max,dif_dec_max,id_max,estSNR] = MWAC(Nband,NestsrcBand,SrcAlpha,SrcDelta,SrcOmega,SrcPhi0,SrcIota,SrcThetaN,SrcAmp,SrcSNR,EstSrc,simParams,yr,0);

% Max over Threshold CC
[gamma,rho,dif_freq_max,dif_ra_max,dif_dec_max,id_max,estSNR1,estSNR2] = ESNMTC(Nband,NestsrcBand,EstSrc1,EstSrc2,simParams,yr,0.90);


% Max coefficients of sources
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

%% Eliminating spurious sources
t = 0.70; % NMTC threshold used to identify sources.
isrc = {}; % identified sources.
r = {}; % rows
c = {}; % columns
for b = 1:Nband
    [r{b},c{b},~] = find(gamma{b} > t); % r is the row of rho, c is the column of rho.
    % in gamma, rows correspond to EstSrc2, columns correspond to EstSrc1.
    % select the identified sources from est. sources.
    for rr = 1:length(r{b})
        isrc{b,rr} = EstSrc2{b,r{b}(rr)};
    end
end


%% Plotting
metric = 'NMTC';
methods = 'xMBLT-iMBLT-vs-initial';
prefix = [estSrc2Dir,filesep,'fig',filesep,metric,'-',methods];
mkdir(prefix);


figname = 'NMTC';

for fig = 1:Nband
    figure
    imagesc(gamma{fig});
    colorbar
    xlabel(estsrc1)
    ylabel(estsrc2)
    title(['Band ',num2str(fig)])
    saveas(gcf,[prefix,filesep,figname,'Band ',num2str(fig)],'png');
    savefig([prefix,filesep,figname,'Band ',num2str(fig)]);
end


% figname2 = 'NMTC_SNR';
% for fig2 = 1:Nband
%     figure
%     plot(estSNR(fig2,:),rho_max{fig2},'ob')
%     xlabel('Estimated SNR')
%     ylabel('NMTC')
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


% figname6 = 'NMTC_freq';
% 
% for fig = 1:Nband
%     figure
%     plot(dif_freq_max(:,fig),rho_max{fig},'ob')
%     xlabel('Difference of Freq. Percentage (%)')
%     ylabel('NMTC')
%     title(['Band ',num2str(fig)])
%     saveas(gcf,[prefix,filesep,figname6,'Band ',num2str(fig)],'png');
%     savefig([prefix,filesep,figname6,'Band ',num2str(fig)]);
% end
% 
% figname7 = 'NMTC_RA';
% 
% for fig = 1:Nband
%     figure
%     plot(dif_ra_max(:,fig),rho_max{fig},'ob')
%     xlabel('Difference of RA Percentage (%)')
%     ylabel('NMTC')
%     title(['Band ',num2str(fig)])
%     saveas(gcf,[prefix,filesep,figname7,'Band ',num2str(fig)],'png');
%     savefig([prefix,filesep,figname7,'Band ',num2str(fig)]);
% end
% 
% figname8 = 'NMTC_DEC';
% 
% for fig = 1:Nband
%     figure
%     plot(dif_dec_max(:,fig),rho_max{fig},'ob')
%     xlabel('Difference of DEC Percentage (%)')
%     ylabel('NMTC')
%     title(['Band ',num2str(fig)])
%     saveas(gcf,[prefix,filesep,figname8,'Band ',num2str(fig)],'png');
%     savefig([prefix,filesep,figname8,'Band ',num2str(fig)]);
% end

figname9 = 'identified sources';

for fig = 1:Nband
    ifreq = arrayfun(@(x) isrc{fig,x}.omega/(2*pi*365*24*3600), 1:length(r{fig}));
    figure
    plot(SrcSNR{fig},SrcOmega{fig}/(2*pi*365*24*3600),'ob',estSNR2(fig,r{fig}),ifreq,'sr')
    text(SrcSNR{fig}+0.5,SrcOmega{fig}/(2*pi*365*24*3600), num2str((1:numel(SrcSNR{fig}))'), 'Color', '#0072BD')
    text(estSNR2(fig,r{fig})-2.5,ifreq, num2str(r{fig}), 'HorizontalAlignment','right', 'Color', '#D95319')
    title(['Identified Sources Band ',num2str(fig)])
    xlabel('SNR')
    ylabel('Frequency(Hz)')
    legend('True Source','Identified Source')
    saveas(gcf,[prefix,filesep,figname9,' Band ',num2str(fig)],'png');
    savefig([prefix,filesep,figname9,' Band ',num2str(fig)]);
end
    


toc