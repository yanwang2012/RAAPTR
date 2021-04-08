% Cross-Correlation Coefficients analysis for two sets of est. sources
% without spliting into different bands and take the union of both set of est.
% sources by combining those highly correlated sources as well keep those not
% correlated sources.

% Author: QYQ
% 09/14/2020

clear;
tic

%% Dir settings
simdataDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/simData';
estSrc1Dir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/xMBLT_set1/results';
estsrc1 = 'supNarxMBLT';
estSrc2Dir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/xMBLT_set2/results';
estsrc2 = 'supNarxMBLTRand1';
Filename = 'GWBsimDataSKASrlz1Nrlz';
ext = '.mat';

%% Files
estSrc1File = dir([estSrc1Dir,filesep,'*',Filename,'*',ext]);
estSrc2File = dir([estSrc2Dir,filesep,'*',Filename,'*',ext]);
simFile = dir([simdataDir,filesep,Filename,'*',ext]);
estSrc2Filename = sort_nat({estSrc2File.name});
estSrc1Filename = sort_nat({estSrc1File.name});
simFilename = sort_nat({simFile.name});

Nestsrc = length(estSrc1File);
Nrlzs = length(simFile);
Nite = Nestsrc/Nrlzs;

for rlz = 1:4%Nrlzs
    load([simdataDir,filesep,simFilename{rlz}]);
    %% Get estimated sources info
    EstSrc2 = {};
    EstSrc1 = {};
    exp = ['\d+_',Filename,num2str(rlz),'(?=_|\.mat)_?\d{0,2}.mat']; % regular expression for multiple realizations
    estSrc1File_tmp = regexp(estSrc1Filename,exp,'match');
    estSrc1File_tmp = estSrc1File_tmp(~cellfun(@isempty, estSrc1File_tmp));
    estSrc2File_tmp = regexp(estSrc2Filename,exp,'match');
    estSrc2File_tmp = estSrc2File_tmp(~cellfun(@isempty, estSrc2File_tmp));
    for k = 1:Nite
        path_to_estimatedDataestSrc2 = [estSrc2Dir,filesep,estSrc2File_tmp{k}{1}];
        path_to_estimatedDataestSrc1 = [estSrc1Dir,filesep,estSrc1File_tmp{k}{1}];
        EstSrc2{k} = ColSrcParams(path_to_estimatedDataestSrc2, simParams.Np);
        EstSrc1{k} = ColSrcParams(path_to_estimatedDataestSrc1, simParams.Np);
    end
    
    % Cross-Corelation
    [gamma,rho,id_max,estSNR1,estSNR2] = ESNMTCW(Nite,EstSrc1,EstSrc2,simParams,yr,0.90);
    
    
    %% find highly correlated sources
    t = 0.80; % NMTC threshold used to select sources.
    [r,c,~] = find(gamma > t); % r is the row of rho, c is the column of rho.
    % in gamma, rows correspond to EstSrc2, columns correspond to EstSrc1.
    % map = [r c]; % create a map between EstSrc1 & EstSrc2.
    
    %% Take union of two Est. sets
    EstSrc1d = EstSrc1(setdiff(1:Nite,c)); % sources in EstSrc1 which indices are not in c
    EstSrc2d = EstSrc2(setdiff(1:Nite,r)); % ...........EstSrc2..........................r
    estSNR1d = estSNR1(setdiff(1:Nite,c));
    estSNR2d = estSNR2(setdiff(1:Nite,r));
    
    UnSrc = [EstSrc1d EstSrc2]; % Union of 2 sets of est. sources with the choice of EstSrc2(r) as the combined highly correlated sources.
    UnSNR = cat(1,estSNR1d,estSNR2);
    
    %% copy files together
    folderName = 'Union';
    [~,simFile_noExt,~] = fileparts(simFilename{rlz});
    NewFolder = [estSrc2Dir,filesep,folderName,filesep,simFile_noExt];
    mkdir(NewFolder);
    for files = 1:Nite
        copyfile([estSrc2Dir,filesep,estSrc2File_tmp{files}{1}],NewFolder); % copy all files from est. sources 2
    end
    
    % Copy EstSrc1 which are not highly correlated with EstSrc2.
    id_diff = setdiff(1:Nite,c); % get indices of est. srouces 1 not in map.
    for i = 1:length(id_diff)
        newFilename = [NewFolder,filesep,'diff_',estSrc1File_tmp{id_diff(i)}{1}];
        copyfile([estSrc1Dir,filesep,estSrc1File_tmp{id_diff(i)}{1}],newFilename);
    end
    
    
    %% Plotting
    % Uncomment to check for a specific realization.
    
    % metric = 'NMTC';
    % methods = 'supNarxMBLT-supNarxMBLTRand1-ALL';
    % prefix = [estSrc2Dir,filesep,folderName,filesep,'fig',filesep,metric,'-',methods];
    % mkdir(prefix);
    %
    % figname = 'NMTC';
    %
    % figure
    % imagesc(gamma);
    % colorbar
    % xlabel(estsrc1)
    % ylabel(estsrc2)
    % title(figname)
    % saveas(gcf,[prefix,filesep,figname],'png');
    % savefig([prefix,filesep,figname]);
    %
    % figname = 'Uinon vs True';
    %
    % figure
    % plot(snr_chr,omega/(365*24*3600*2*pi),'ob')
    % hold on
    % for i = 1:length(UnSrc)
    %     plot(UnSNR(i),UnSrc{i}.omega/(365*24*3600*2*pi),'sr')
    % end
    % hold off
    %
    % title('Union vs. True')
    % xlabel('SNR')
    % ylabel('Frequency')
    % legend('True','Union')
    % saveas(gcf,[prefix,filesep,figname],'png')
    % savefig([prefix,filesep,figname])
end

toc

%END