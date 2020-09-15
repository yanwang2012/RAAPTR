% Cross-Correlation Coefficients analysis for two sets of est. sources
% without spliting into different bands and take the union of both set of est.
% sources by combining those highly correlated sources as well keep those not
% correlated sources.

% Author: QYQ
% 09/14/2020

clear;
tic

%% Dir settings
simdataDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands';
estSrc1Dir = '/Users/qyq/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands/Results_supNar/GWBsimDataSKASrlz1Nrlz3_xMBLT/results';
estsrc1 = 'supNarxMBLT';
estSrc2Dir = '/Users/qyq/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/2bands/Results_supNar_rand1/GWBsimDataSKASrlz1Nrlz3_xMBLT/results';
estsrc2 = 'supNarxMBLTRand1';
Filename = 'GWBsimDataSKASrlz1Nrlz3';
ext = '.mat';

%% Files
estSrc1File = dir([estSrc1Dir,filesep,'*',Filename,'*',ext]);
estSrc2File = dir([estSrc2Dir,filesep,'*',Filename,'*',ext]);
Nestsrc = length(estSrc2File);
estSrc2Filename = sort_nat({estSrc2File.name});
estSrc1Filename = sort_nat({estSrc1File.name});
simFile = [simdataDir,filesep,Filename,ext];
load(simFile);
%% Get estimated sources info
EstSrc2 = {};
EstSrc1 = {};

for k = 1:Nestsrc
    path_to_estimatedDataestSrc2 = [estSrc2Dir,filesep,char(estSrc2Filename(k))];
    path_to_estimatedDataestSrc1 = [estSrc1Dir,filesep,char(estSrc1Filename(k))];
    EstSrc2{k} = ColSrcParams(path_to_estimatedDataestSrc2);
    EstSrc1{k} = ColSrcParams(path_to_estimatedDataestSrc1);
end

% Cross-Corelation
[gamma,rho,id_max,estSNR1,estSNR2] = ESNMTCW(Nestsrc,EstSrc1,EstSrc2,simParams,yr,0.90);


%% find highly correlated sources
t = 0.80; % NMTC threshold used to select sources.
[r,c,~] = find(gamma > t); % r is the row of rho, c is the column of rho.
% in gamma, rows correspond to EstSrc2, columns correspond to EstSrc1.
map = [r c]; % create a map between EstSrc1 & EstSrc2.

%% Take union of two Est. sets
EstSrc1d = EstSrc1(setdiff(1:Nestsrc,c)); % sources in EstSrc1 which indices are not in c
EstSrc2d = EstSrc2(setdiff(1:Nestsrc,r)); % ...........EstSrc2..........................r
estSNR1d = estSNR1(setdiff(1:Nestsrc,c));
estSNR2d = estSNR2(setdiff(1:Nestsrc,r));

UnSrc = [EstSrc1d EstSrc2d EstSrc2(r)]; % Union of 2 sets of est. sources with the choice of EstSrc2(r) as the combined highly correlated sources.
UnSNR = cat(1,estSNR1d,estSNR2d,estSNR2(r));

%% Plotting
metric = 'NMTC';
methods = 'supNarxMBLT-supNarxMBLTRand1-ALL';
prefix = [estSrc2Dir,filesep,'fig',filesep,metric,'-',methods];
mkdir(prefix);

figname = 'NMTC';

figure
imagesc(gamma);
colorbar
xlabel(estsrc1)
ylabel(estsrc2)
title(figname)
saveas(gcf,[prefix,filesep,figname],'png');
savefig([prefix,filesep,figname]);

figname = 'Uinon vs True';

figure
plot(snr_chr,omega/(365*24*3600*2*pi),'ob')
hold on
for i = 1:length(UnSrc)
   plot(UnSNR(i),UnSrc{i}.omega/(365*24*3600*2*pi),'sr') 
end
hold off

title('Union vs. True')
xlabel('SNR')
ylabel('Frequency')
legend('True','Union')
saveas(gcf,[prefix,filesep,figname],'png')
savefig([prefix,filesep,figname])


toc
%END