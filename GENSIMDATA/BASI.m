% Illustration for Band Automate Selection
% Author QYQ
% 08/09/2021
%% load data
clear;
xMBLT1Data = '/Users/yiqianqian/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/diff_srlz_cos_xMBLT1_results';
xMBLT2Data = '/Users/yiqianqian/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/diff_srlz_cos_xMBLT2_results';
UnionData = '/Users/yiqianqian/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/diff_srlz_cos_xMBLT2_results/Union';
searchParamsDir = '/Users/yiqianqian/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/searchParams/2bands/superNarrow';
BandOptDir = '/Users/yiqianqian/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/searchParams/Band_opt_cos';
simDataDir = '/Users/yiqianqian/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/simData/diff_srlz_cos';
FileName = 'GWBsimDataSKASrlz*Nrlz1';
simFiles = dir([simDataDir,filesep,FileName,'.mat']);
simFileNames = sort_nat({simFiles.name});

Nrlz = length(simFileNames);

searchFiles = dir([searchParamsDir,filesep,'searchParams*.mat']);
searchFiles = {searchFiles.name};
exp1 = 'searchParams\d\.mat';
searchFiles1 = regexp(searchFiles,exp1,'match');
searchFiles1 = searchFiles1(~cellfun(@isempty,searchFiles1));

exp2 = 'searchParamsRand\d\.mat';
searchFiles2 = regexp(searchFiles,exp2,'match');
searchFiles2 = searchFiles2(~cellfun(@isempty,searchFiles2));

nband = length(searchFiles1);
Band1 = []; % bands for set 1
Band2 = []; % bands for set 2

for band = 1:nband
    
    search1 = load([searchParamsDir,filesep,searchFiles1{band}{1}]);
    search2 = load([searchParamsDir,filesep,searchFiles2{band}{1}]);
    
    Band1 = [Band1 search1.searchParams.angular_velocity(1)/(2*pi*365*24*3600)];
    Band2 = [Band2 search2.searchParams.angular_velocity(1)/(2*pi*365*24*3600)];
    
end

for rlz=1:Nrlz
    [~,baseName,~] = fileparts(simFileNames{rlz});
    load([simDataDir,filesep,baseName,'.mat'])
    
    files1 = dir([xMBLT1Data,filesep,'*_',baseName,'*.mat']);
    files1Name = sort_nat({files1.name});
    files2 = dir([xMBLT2Data,filesep,'*_',baseName,'*.mat']);
    files2Name = sort_nat({files2.name});
    filesU = dir([UnionData,filesep,baseName,filesep,'*',baseName,'*.mat']);
    filesUName = sort_nat({filesU.name});
    
    nFiles = length(files1);
    nFilesU = length(filesU);
    
    searchFilesU = dir([BandOptDir,filesep,baseName,filesep,'searchParams_Nyquist*.mat']);
    searchFilesU = {searchFilesU.name};
    expU = 'searchParams_Nyquist\d\.mat';
    searchFilesU = regexp(searchFilesU,expU,'match');
    searchFilesU = searchFilesU(~cellfun(@isempty,searchFilesU));
    BandU = [];
    nBandU = length(searchFilesU);
    
    for band = 1:nBandU
        searchU = load([BandOptDir,filesep,baseName,filesep,searchFilesU{band}{1}]);
        BandU = [BandU searchU.searchParams.angular_velocity(1)/(2*pi*365*24*3600)];
    end
    
    snr1 = [];
    freq1 = [];
    snr2 = [];
    freq2 = [];
    snrU = [];
    freqU = [];
    
    for n = 1:nFiles
        estSrc1 = ColSrcParams([xMBLT1Data,filesep,files1Name{n}],simParams.Np);
        estSrc2 = ColSrcParams([xMBLT2Data,filesep,files2Name{n}],simParams.Np);
        [snr1_tmp,~] = Amp2Snr(estSrc1, simParams, yr);
        snr1 = [snr1 snr1_tmp];
        freq1 = [freq1 estSrc1.omega/(2*pi*365*24*3600)];
        
        [snr2_tmp,~] = Amp2Snr(estSrc2, simParams, yr);
        snr2 = [snr2 snr2_tmp];
        freq2 = [freq2 estSrc2.omega/(2*pi*365*24*3600)];
    end
    
    for u = 1:nFilesU
        estSrcU = ColSrcParams([UnionData,filesep,baseName,filesep,filesUName{u}],simParams.Np);
        [snrU_tmp,~] = Amp2Snr(estSrcU, simParams, yr);
        snrU = [snrU snrU_tmp];
        freqU = [freqU estSrcU.omega/(2*pi*365*24*3600)];
    end
    
    %% Plot
    figure
    ax1 = subplot(1,3,1);
    plot(snr1,freq1,'o')
    hold on
    for band = 1:nband
        plot([0, max(cat(2,snr1,snr2,snrU))],[Band1(band), Band1(band)],'r')
    end
    hold off
    title('Band set 1')
    xlabel('SNR')
    ylabel('Frequency [Hz]')
    xlim([0 max(cat(2,snr1,snr2,snrU))])
    
    ax2 = subplot(1,3,2);
    plot(snr2,freq2,'o')
    hold on
    for band = 1:nband
        plot([0, max(cat(2,snr1,snr2,snrU))],[Band2(band), Band2(band)],'r')
    end
    hold off
    title('Band set 2')
    xlabel('SNR')
    ylabel('Frequency [Hz]')
    xlim([0 max(cat(2,snr1,snr2,snrU))])
    
    ax3 = subplot(1,3,3);
    plot(snrU,freqU,'o')
    hold on
    for band = 1:nBandU
        plot([0, max(cat(2,snr1,snr2,snrU))],[BandU(band), BandU(band)],'r')
    end
    hold off
    title('Union')
    xlabel('SNR')
    ylabel('Frequency [Hz]')
    xlim([0 max(cat(2,snr1,snr2,snrU))])
    
    linkaxes([ax1,ax2,ax3],'xy')
    
    mkdir([UnionData,filesep,'fig'])
    saveas(gcf,[UnionData,filesep,'fig',filesep,'BASI-',num2str(rlz),'.png'])
    save([UnionData,filesep,'fig',filesep,'BASI-',num2str(rlz)])
    
    close all
    
end