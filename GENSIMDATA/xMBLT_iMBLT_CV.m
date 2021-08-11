% Plot xMBLT, iMBLT and estmate sources after cross validation
% Author QYQ
% 08/09/2021
clear;
%% load datasets
simdataDir = '/Users/yiqianqian/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/simData/Band_opt_diff';
xMBLTDir = '/Users/yiqianqian/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/results_diff_opt_xMBLT';
iMBLTDir = '/Users/yiqianqian/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/results_diff_opt_iMBLT';
CVDir = '/Users/yiqianqian/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/results_diff_opt_iMBLT';
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
    Nsrc = length(CVSrc.CnfrmSrc_SNR);
    CV_SNR = [];
    CV_Freq = [];
    for src = 1:Nsrc
    [CV_SNR_tmp,~] = Amp2Snr(CVSrc.CnfrmSrc_SNR{src},simParams,yr);
    CV_SNR = [CV_SNR CV_SNR_tmp];
    CV_Freq_tmp = CVSrc.CnfrmSrc_SNR{src}.omega/(2*pi*365*24*3600);
    CV_Freq = [CV_Freq CV_Freq_tmp];
    end
    % figure
    figure
    ax1 = subplot(1,3,1);
    plot(xMBLT_SNR,xMBLT_Freq,'o')
    title('xBSE')
    xlabel('SNR');
    ylabel('Frequency [Hz]')
    
    ax2 = subplot(1,3,2);
    plot(iMBLT_SNR,iMBLT_Freq,'o')
    title('iBSE')
    xlabel('SNR')
    ylabel('Frequency [Hz]')
    
    ax3 = subplot(1,3,3);
    plot(CV_SNR,CV_Freq,'o')
    title('CV')
    xlabel('SNR')
    ylabel('Frequency [Hz]')
    linkaxes([ax1,ax2,ax3],'xy')
    
    saveas(gcf,[iMBLTDir,filesep,'fig',filesep,folderNames{Srlz}{1},filesep,'Together.png'])
    savefig([iMBLTDir,filesep,'fig',filesep,folderNames{Srlz}{1},filesep,'Together'])
    
    close all;
end