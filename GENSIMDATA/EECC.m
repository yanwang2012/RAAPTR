% Cross-Correlation Coefficients Matrix for Est. & Est. sources.

% Author: QYQ
% 05/13/2020

clear;
tic

%% Dir settings
serchParamsDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/searchParams/Band_opt/New';
simdataDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/Band_opt/simData';
estSrc1Dir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/Band_opt/results_New_xMBLT';
estsrc1 = 'Band-opt-xMBLT';
estSrc2Dir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/Band_opt/iMBLT/results'; 
estsrc2 = 'Band-opt-iMBLT';
Filename = 'GWBsimDataSKASrlz1Nrlz3';
ext = '.mat';

%% Files
paraFile = dir([serchParamsDir,filesep,'searchParams*',ext]);
simFile = [simdataDir,filesep,Filename,ext];
estSrc1File = dir([estSrc1Dir,filesep,'*',Filename,'*',ext]);
estSrc2File = dir([estSrc2Dir,filesep,'*',Filename,'*',ext]);
Nestsrc2 = length(estSrc2File);

paraFilename = sort_nat({paraFile.name});
exp = 'searchParams_Nyquist\d.mat'; % regular expressions for desire file names
paraFilename = regexp(paraFilename,exp,'match');
paraFilename = paraFilename(~cellfun(@isempty,paraFilename)); % get rid of empty cells
Nband = length(paraFilename);

estSrc2Filename = sort_nat({estSrc2File.name});
estSrc1Filename = sort_nat({estSrc1File.name});
NestSrc1band1 = sum(startsWith(estSrc1Filename,'1_'));
NestSrc1band2 = sum(startsWith(estSrc1Filename,'2_'));
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
    load([serchParamsDir,filesep,char(paraFilename{i})]);
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
simSrc = struct('SrcSNR',SrcSNR,'SrcAlpha',SrcAlpha,'SrcDelta',SrcDelta,'SrcAmp',SrcAmp,...
    'SrcIota',SrcIota,'SrcOmega',SrcOmega,'SrcPhi0',SrcPhi0,'SrcThetaN',SrcThetaN); % Simulated sources parameters

%% Get estimated sources info
NestSrc2Band = Nestsrc2/Nband; % number of sources in a band.
BandSrc = struct('NestSrc1band1',NestSrc1band1,'NestSrc1band2',NestSrc1band2,'NestSrc2Band',NestSrc2Band); % for Union which contains uneven No. of sources in each band.
EstSrc2 = {};
EstSrc1 = {};
for band = 1:Nband
    for k = 1:NestSrc2Band
        path_to_estimatedDataestSrc2 = [estSrc2Dir,filesep,char(estSrc2Filename((band - 1) * NestSrc2Band + k))];
        path_to_estimatedDataestSrc1 = [estSrc1Dir,filesep,char(estSrc1Filename((band - 1) * NestSrc2Band + k))];
        
        EstSrc2{band,k} = ColSrcParams(path_to_estimatedDataestSrc2,simParams.Np);
        EstSrc1{band,k} = ColSrcParams(path_to_estimatedDataestSrc1,simParams.Np);
    end
end

% Union only!! for uneven No. of band src.
% for band = 1:Nband
%     switch band
%         case 1
%             for k = 1:NestSrc1band1
%             path_to_estimatedDataestSrc1 = [estSrc1Dir,filesep,char(estSrc1Filename(k))];
%             EstSrc1{band,k} = ColSrcParams(path_to_estimatedDataestSrc1,simParams.Np);
%             end
%         case 2
%             for k = 1:NestSrc1band2
%                 path_to_estimatedDataestSrc1 = [estSrc1Dir,filesep,char(estSrc1Filename(k + NestSrc1band1))];
%                 EstSrc1{band,k} = ColSrcParams(path_to_estimatedDataestSrc1,simParams.Np);
%             end
%     end
% end
            


%% Cross-Corelation

% Max Weighted CC
% [rho,rho_max,dif_freq_max,dif_ra_max,dif_dec_max,id_max,estSNR] = MAC(Nband,NestsrcBand,SrcAlpha,SrcDelta,SrcOmega,SrcPhi0,SrcIota,SrcThetaN,SrcAmp,SrcSNR,EstSrc,simParams,yr,'snr');

% Max Weighted Ave. CC
% [rho,rho_max,dif_freq_max,dif_ra_max,dif_dec_max,id_max,estSNR] = MWAC(Nband,NestsrcBand,SrcAlpha,SrcDelta,SrcOmega,SrcPhi0,SrcIota,SrcThetaN,SrcAmp,SrcSNR,EstSrc,simParams,yr,0);

% Max over Threshold CC
[gamma,rho,dif_freq_max,dif_ra_max,dif_dec_max,id_max,estSNR1,estSNR2] = ESNMTC(Nband,BandSrc,EstSrc1,EstSrc2,simParams,yr,0.90);
% shape of gamma is (EstSrc2,EstSrc1) and we check all the EstSrc2 for each
% EstSrc1, i.e. we are doing ESNMTC along the column direction.

%% Eliminating spurious sources
t = 0.70; % NMTC threshold used to identify sources.
idsrc = {}; % identified sources.
r = {}; % rows
c = {}; % columns
for b = 1:Nband
    [r{b},c{b},~] = find(gamma{b} > t); % r is the row of gamma, c is the column of gamma.
    % in gamma, rows correspond to EstSrc2, columns correspond to EstSrc1.
    % select the identified sources from est. sources.
    for rr = 1:length(r{b})
        idsrc{b,rr} = EstSrc2{b,r{b}(rr)};
    end
end

NidsrcBand = zeros(Nband,1);
for idb = 1:Nband
    NidsrcBand(idb) = sum(~cellfun('isempty',idsrc(idb,:))); % # of identified sources in each band
end

idMethods = ['id-',estsrc1,'-vs-',estsrc2];
idDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/Band_opt';
idFolder = [idDir,filesep,idMethods];
mkdir(idFolder);
save([idFolder,filesep,'IdentifiedSrc.mat'],'idsrc','NidsrcBand')

%% Plotting
metric = 'NMTC';
methods = 'Band-opt-xMBLT vs Band-opt-iMBLT';
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
    ifreq = arrayfun(@(x) idsrc{fig,x}.omega/(2*pi*365*24*3600), 1:length(r{fig}));
    figure
    plot(SrcSNR{fig},SrcOmega{fig}/(2*pi*365*24*3600),'ob',estSNR2(fig,r{fig}),ifreq,'sr')
    text(SrcSNR{fig}+0.5,SrcOmega{fig}/(2*pi*365*24*3600), num2str((1:numel(SrcSNR{fig}))'), 'Color', '#0072BD')
    text(estSNR2(fig,r{fig})-2.5,ifreq, num2str((1:numel((r{fig})))'), 'HorizontalAlignment','right', 'Color', '#D95319')
    title(['Identified Sources Band ',num2str(fig)])
    xlabel('SNR')
    ylabel('Frequency(Hz)')
    legend('True Source','Identified Source')
    saveas(gcf,[prefix,filesep,figname9,' Band ',num2str(fig)],'png');
    savefig([prefix,filesep,figname9,' Band ',num2str(fig)]);
end
    


toc