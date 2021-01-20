% Cross-Correlation Coefficients Matrix
% True sources vs. Identified sources

% Author: QYQ
% 1/6/2021

clear;
tic

%% Dir settings
searchParamsDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/searchParams/2bands/superNarrow';
simdataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11';
identifydataDir = '/Users/qyq/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands/SuperNarrow/Union2_xMBLT/Union2_iMBLT_after20';
Filename = 'GWBsimDataSKASrlz1Nrlz3';
identifyFilename = 'identifiedSrc';
ext = '.mat';

%% Files
paraFile = dir([searchParamsDir,filesep,'searchParams*',ext]);
simFile = [simdataDir,filesep,Filename,ext];
idFile = [identifydataDir,filesep,identifyFilename,ext];

paraFilename = sort_nat({paraFile.name});
exp = 'searchParams\d.mat'; % regular expressions for desire file names
paraFilename = regexp(paraFilename,exp,'match');
paraFilename = paraFilename(~cellfun(@isempty,paraFilename)); % get rid of empty cells
Nband = length(paraFilename);


load(simFile);
load(idFile);

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
    load([searchParamsDir,filesep,char(paraFilename{i})]);
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
simSrc = struct('SrcSNR',SrcSNR,'SrcAlpha',SrcAlpha,'SrcDelta',SrcDelta,'SrcAmp',SrcAmp,...
    'SrcIota',SrcIota,'SrcOmega',SrcOmega,'SrcPhi0',SrcPhi0,'SrcThetaN',SrcThetaN); % Simulated sources parameters

%% Get identified sources info
% idsrcBand1 = sum(~cellfun('isempty',idsrc(1,:))); % number of sources in a band.
% idsrcBand2 = sum(~cellfun('isempty',idsrc(2,:)));
% idsrcBand = struct('Band1',idsrcBand1,'Band2',idsrcBand2);
idsrcBand = NidsrcBand;

%% Cross-Corelation

% Max Weighted CC
% [rho,rho_max,dif_freq_max,dif_ra_max,dif_dec_max,id_max,estSNR] = MWC(Nband,NestsrcBand,SrcAlpha,SrcDelta,SrcOmega,SrcPhi0,SrcIota,SrcThetaN,SrcAmp,SrcSNR,EstSrc,simParams,yr,'snr');

% Max Weighted Ave. CC
% [rho,rho_max,dif_freq_max,dif_ra_max,dif_dec_max,id_max,estSNR] = MWAC(Nband,NestsrcBand,SrcAlpha,SrcDelta,SrcOmega,SrcPhi0,SrcIota,SrcThetaN,SrcAmp,SrcSNR,EstSrc,simParams,yr,'snr');

% Max over Threshold CC
% [rho,rho_max,dif_freq_max,dif_ra_max,dif_dec_max,id_max,estSNR] = MTC(Nband,NestsrcBand,SrcAlpha,SrcDelta,SrcOmega,SrcPhi0,SrcIota,SrcThetaN,SrcAmp,EstSrc,simParams,yr,0.85);

% Normalized MTC
[rho,rho_max,dif_freq_max,dif_ra_max,dif_dec_max,id_max,estSNR] = NMTC(Nband,idsrcBand,simSrc,idsrc,simParams,yr,0.90);


% Minimum distance Maximum CC.
% [rho,rho_max,dif_freq_max,dif_ra_max,dif_dec_max,id_max,estSNR] = MinDMaxC(Nband,NestsrcBand,SrcAlpha,SrcDelta,SrcOmega,SrcPhi0,SrcIota,SrcThetaN,SrcAmp,EstSrc,simParams,yr);

%% Save matched true sources
% save sky locations
matched_alpha = []; % right ascension
matched_dec = []; % declination
matched_snr = [];

for band = 1:Nband
    matched_alpha = [matched_alpha SrcAlpha{band}(id_max(id_max(:,band) ~= 0, band))]; % exclude 0 elements
    matched_dec = [matched_dec SrcDelta{band}(id_max(id_max(:,band) ~= 0, band))];
    matched_snr = [matched_snr SrcSNR{band}(id_max(id_max(:,band) ~= 0, band))];
end

save([identifydataDir,filesep,'Matched_Sources.mat'],'id_max','matched_alpha','matched_dec','matched_snr',...
    'SrcAlpha','SrcDelta');
%% Plotting
metric = 'NMTC';
methods = 'True vs Union2-xMBLT-iMBLT-identified';
prefix = [identifydataDir,filesep,'fig',filesep,metric,'-',methods];
mkdir(prefix);

figname1 = metric;
for fig = 1:Nband
    figure
    imagesc(rho{fig});
    a = colorbar;
    xlabel('True sources')
    ylabel('Identified sources')
    ylabel(a,'Corss-Correlation Coefficients')
    title(['Band ',num2str(fig)])
    saveas(gcf,[prefix,filesep,figname1,'Band ',num2str(fig)],'png');
    savefig([prefix,filesep,figname1,'Band ',num2str(fig)]);
end


figname2 = [metric,'-SNR'];
for fig2 = 1:Nband
    switch fig2
        case 1
            N = idsrcBand1;
        case 2
            N = idsrcBand2;
    end
    figure
    plot(estSNR(fig2,1:N),rho_max{fig2},'ob')
    xlabel('Estimated SNR')
    ylabel(metric)
    title(['Band ',num2str(fig2)])
    saveas(gcf,[prefix,filesep,figname2,'Band ',num2str(fig2)],'png');
    savefig([prefix,filesep,figname2,'Band ',num2str(fig2)]);
end

% figname3 = [metric,'identified sources'];
%
% for fig = 1:Nband
%     ifreq = arrayfun(@(x) isrc{fig,x}.omega/(2*pi*365*24*3600), 1:length(r{fig}));
%     figure
%     plot(SrcSNR{fig},SrcOmega{fig}/(2*pi*365*24*3600),'ob',estSNR(fig,r{fig}),ifreq,'sr')
%     text(SrcSNR{fig}+0.5,SrcOmega{fig}/(2*pi*365*24*3600), num2str((1:numel(SrcSNR{fig}))'), 'Color', '#0072BD')
%     text(estSNR(fig,r{fig})-2,ifreq, num2str(r{fig}), 'HorizontalAlignment','right', 'Color', '#D95319')
%     title(['Identified Sources Band ',num2str(fig)])
%     xlabel('SNR')
%     ylabel('Frequency(Hz)')
%     legend('True Source','Identified Source')
%     saveas(gcf,[prefix,filesep,figname3,'Band ',num2str(fig)],'png');
%     savefig([prefix,filesep,figname3,'Band ',num2str(fig)]);
% end

% figname6 = [metric,'-freq'];
% 
% for fig = 1:Nband
%     switch fig
%         case 1
%             N = idsrcBand1;
%         case 2
%             N = idsrcBand2;
%     end
%     figure
%     plot(dif_freq_max(1:N,fig),rho_max{fig},'ob')
%     xlabel('Difference of Freq. Percentage (%)')
%     ylabel(metric)
%     title(['Band ',num2str(fig)])
%     saveas(gcf,[prefix,filesep,figname6,'Band ',num2str(fig)],'png');
%     savefig([prefix,filesep,figname6,'Band ',num2str(fig)]);
% end
% 
% figname7 = [metric,'-RA'];
% 
% for fig = 1:Nband
%     switch fig
%         case 1
%             N = idsrcBand1;
%         case 2
%             N = idsrcBand2;
%     end
%     figure
%     plot(dif_ra_max(1:N,fig),rho_max{fig},'ob')
%     xlabel('Difference of RA Percentage (%)')
%     ylabel(metric)
%     title(['Band ',num2str(fig)])
%     saveas(gcf,[prefix,filesep,figname7,'Band ',num2str(fig)],'png');
%     savefig([prefix,filesep,figname7,'Band ',num2str(fig)]);
% end
% 
% figname8 = [metric,'-DEC'];
% 
% for fig = 1:Nband
%     switch fig
%         case 1
%             N = idsrcBand1;
%         case 2
%             N = idsrcBand2;
%     end
%     figure
%     plot(dif_dec_max(1:N,fig),rho_max{fig},'ob')
%     xlabel('Difference of DEC Percentage (%)')
%     ylabel(metric)
%     title(['Band ',num2str(fig)])
%     saveas(gcf,[prefix,filesep,figname8,'Band ',num2str(fig)],'png');
%     savefig([prefix,filesep,figname8,'Band ',num2str(fig)]);
% end


toc