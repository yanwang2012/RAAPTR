% Plot xBSE + iBSE - 20 SNRCut for 6 source realizations

% Author QYQ 08/20/2021
clear;

%% Config
DataDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/results_diff_opt_iMBLT';
simDataDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/simData/Band_opt_diff';
CVfileName = 'Confirmed_Src_xMBLT_SNR20.mat';
% get the source realizations
files = dir(simDataDir);
folderNames = sort_nat({files.name});
exp = '^GWBsimData.*';
folderNames = regexp(folderNames,exp,'match');
folderNames = folderNames(~cellfun(@isempty,folderNames));
NSrlz = length(folderNames);

%% Plot
figure
[ha,~] = tight_subplot(2,3,[.1 .05],[.1 .1], [.1 .05]);
for rlz = 1:NSrlz
    [~,baseName,~] = fileparts(folderNames{rlz}{1});
    load([simDataDir,filesep,baseName,'.mat']);
    CV = load([DataDir,filesep,baseName,filesep,CVfileName]);
    confirm_src = CV.confirm_src;
    CV_Freq = [];
    CV_SNR = [];
    Nsrc = length(confirm_src);
    for src = 1:Nsrc
       [CV_SNR_tmp,~] = Amp2Snr(confirm_src{src},simParams,yr); 
       CV_SNR = [CV_SNR CV_SNR_tmp];
       CV_Freq = [CV_Freq confirm_src{src}.omega/(2*pi*24*365*3600)];
    end
    axes(ha(rlz))
    plot(snr_chr,omega/(2*pi*24*365*3600),'bo')
    hold on
    plot(CV_SNR,CV_Freq,'r.','MarkerSize',10)
    hold off
    xlabel('SNR')
end
yticklabels(ha([2,3,5,6]),'')
ylabel(ha([1,4]),'Frequency [Hz]')
saveas(gcf,[DataDir,filesep,'fig',filesep,'xBSE_iBSE_No_SNR.png'])
savefig([DataDir,filesep,'fig',filesep,'xBSE_iBSE_No_SNR'])

% EOS