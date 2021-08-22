% Plot xMBLT, iMBLT and estmate sources after cross validation
% Author QYQ
% 08/09/2021
clear;
%% load datasets
simdataDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/simData/Band_opt_diff';
xMBLTDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/results_diff_opt_xMBLT';
iMBLTDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/results_diff_opt_iMBLT';
CVDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/results_diff_opt_iMBLT';
% get the source realizations
files = dir(xMBLTDir);
folderFlag = [files.isdir];
folderNames = files(folderFlag);
folderNames = sort_nat({folderNames.name});
exp = 'GWBsimData.*';
folderNames = regexp(folderNames,exp,'match');
folderNames = folderNames(~cellfun(@isempty,folderNames));
NSrlz = length(folderNames);

% plot
for Srlz = 1:NSrlz
    load([simdataDir,filesep,folderNames{Srlz}{1},'.mat'])
    xMBLT = load([xMBLTDir,filesep,folderNames{Srlz}{1},filesep,'EstSrc.mat']);
    iMBLT = load([iMBLTDir,filesep,folderNames{Srlz}{1},filesep,'EstSrc.mat']);
    CVSrc = load([CVDir,filesep,folderNames{Srlz}{1},filesep,'Confirmed_Src_xMBLT_SNR20.mat']);
    xMBLT_SNR = xMBLT.sx;
    xMBLT_Freq = xMBLT.sy;
    iMBLT_SNR = iMBLT.sx;
    iMBLT_Freq = iMBLT.sy;
    Nsrc = length(CVSrc.confirm_src);
    CV_SNR = [];
    CV_Freq = [];
    for src = 1:Nsrc
    [CV_SNR_tmp,~] = Amp2Snr(CVSrc.confirm_src{src},simParams,yr);
    CV_SNR = [CV_SNR CV_SNR_tmp];
    CV_Freq_tmp = CVSrc.confirm_src{src}.omega/(2*pi*365*24*3600);
    CV_Freq = [CV_Freq CV_Freq_tmp];
    end
    %% Plot
    figure
    [ha,~] = tight_subplot(1,3,[.03 .01],[.1 .1], [.1 .05]);
    axes(ha(1))
    plot(xMBLT_SNR(xMBLT_SNR>20),xMBLT_Freq(xMBLT_SNR>20),'.','MarkerSize',15)
    hold on
    plot(xMBLT_SNR(xMBLT_SNR<=20), xMBLT_Freq(xMBLT_SNR<=20),'o')
    hold off
    %title('xBSE')
    xlabel('SNR');
    ylabel('Frequency [Hz]')
    
    axes(ha(2))
    plot(iMBLT_SNR(iMBLT_SNR>20),iMBLT_Freq(iMBLT_SNR>20),'.','MarkerSize',15)
    hold on
    plot(iMBLT_SNR(iMBLT_SNR<=20), iMBLT_Freq(iMBLT_SNR<=20),'o')
    hold off
    %title('iBSE')
    xlabel('SNR')
    %ylabel('Frequency [Hz]')
    
    axes(ha(3))
    plot(CV_SNR(CV_SNR>20),CV_Freq(CV_SNR>20),'.','MarkerSize',15)
    hold on
    plot(CV_SNR(CV_SNR<=20), CV_Freq(CV_SNR<=20),'o')
    hold off
    %title('CV')
    xlabel('SNR')
    %ylabel('Frequency [Hz]')
    set(ha(2:3),'YTickLabel','')
    %linkaxes([ax1,ax2,ax3],'xy')
    
    saveas(gcf,[iMBLTDir,filesep,'fig',filesep,folderNames{Srlz}{1},filesep,'Together.png'])
    savefig([iMBLTDir,filesep,'fig',filesep,folderNames{Srlz}{1},filesep,'Together'])
    
    close all;
end