% script to plot identified sources vs. True sources
% (1) 2D skymap
% (2) 3D skymap

% Author: QYQ
% 1/7/2021
clear;
%% Load data
simDataDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/Band_opt/simData';
idDataDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/Band_opt/iMBLT/results';
Filename = 'GWBsimDataSKASrlz1Nrlz3';
identifyFilename = 'SelectedSrc';
ext = '.mat';

simFile = [simDataDir,filesep,Filename,ext];
idFile = [idDataDir,filesep,identifyFilename,ext];
id2true = [idDataDir,filesep,'Matched_Sources.mat'];
load(simFile);
load(idFile);
load(id2true);

%% Sky location
simRA_nm = []; % not matched true sources
simDec_nm = [];

idRA = [];
idDec = [];
idSNR = [];
Nband = length(NidsrcBand);
idBandSrc = zeros(Nband,1);

% get No. of identified sources in each band
for band = 1:Nband
idBandSrc(band) = sum(~cellfun('isempty',idsrc(band,:)));
end

for b = 1:Nband
    N = idBandSrc(b);
    idx = setdiff(1:length(SrcAlpha{b}),id_max(:,b)); % get the complementary index
    for i = 1:N
        idRA = [idRA idsrc{b,i}.alpha];
        idDec = [idDec idsrc{b,i}.delta];
        [idSNR_tmp,~] = Amp2Snr(idsrc{b,i},simParams,yr);
        idSNR = [idSNR idSNR_tmp];
        simRA_nm = [simRA_nm SrcAlpha{b}(idx)];
        simDec_nm = [simDec_nm SrcDelta{b}(idx)];
    end
end

%% plot
figure
% plot(simRA,simDec,'o','MarkerEdgeColor','##D4D8D9') % plot true sources
scatter(simRA_nm,simDec_nm,'o','MarkerEdgeColor','#D4D8D9') % use scatter to make marker size change

hold on

% plot(idRA,idDec,'rs') % plot identified sources
% scatter(idRA,idDec,idSNR,'rs') % with different size
scatter(idRA,idDec,[],idSNR, 'filled') % with different color

% pause

% plot(matched_alpha,matched_dec,'ob') % plot truly matched true sources
% scatter(matched_alpha,matched_dec,matched_snr,'ob')
scatter(matched_alpha,matched_dec,[],matched_snr,'s')

% pause

for j = 1:length(idRA)
    plot([idRA(j),matched_alpha(j)],[idDec(j),matched_dec(j)],'Color','m') % connect identified and matched sources
%     pause
end


xlabel('RA')
ylabel('DEC')
title('Sky Location')
colormap turbo
c = colorbar;
ylabel(c,'SNR','FontSize',14)
ylim([-2 inf])
legend({'True Sources', 'Identified Sources','Matched True Source','Matched & Identi.'},'Location','southeast')
saveas(gcf,[idDataDir,filesep,'fig',filesep,'SkyLocationC.png'])
savefig([idDataDir,filesep,'fig',filesep,'SkyLocationC'])

%% 3D sphere plot
figure
[x,y,z] = sphere; % generate a unit sphere

% convert sources coord. to spherical
[id_x,id_y,id_z] = sph2cart(idRA,idDec,1);
[sim_x,sim_y,sim_z] = sph2cart(simRA_nm,simDec_nm,1);
[matched_x,matched_y,matched_z] = sph2cart(matched_alpha,matched_dec,1);

surf(x,y,z,'LineStyle','--','FaceColor','white')
axis equal
% pause
hold on
scatter3(id_x,id_y,id_z,40,idSNR,'filled')
scatter3(sim_x,sim_y,sim_z,'o','MarkerEdgeColor','#888E8F')
scatter3(matched_x,matched_y,matched_z,[],matched_snr,'s')
c = colorbar;
ylabel(c,'SNR','FontSize',14)

% connect identified sources with matched true sources with great circle on
% sphere
center = zeros(1,3); % center of sphere
for src = 1:length(idRA)
[v] = GreatCircle(center,[id_x(src),id_y(src),id_z(src)],[matched_x(src),matched_y(src),matched_z(src)],1);
plot3(v(1,:),v(2,:),v(3,:),'m')
pause
end

legend({'Grid','Identified Sources','True Sources','Matched True Source','Matched & Identi.'},'Location','southeast')
saveas(gcf,[idDataDir,filesep,'fig',filesep,'SkyLocation3D.png'])
savefig([idDataDir,filesep,'fig',filesep,'SkyLocation3D'])

%END