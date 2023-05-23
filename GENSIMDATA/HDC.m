% calculate the average CorrelationCoefficients between every two pulsar

% author: QYQ
% 2022/04/14
%% Main
clear;
tic
simDataDir = '/Users/yiqianqian/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/HellingsDowns_test';
files = dir([simDataDir,filesep,'GWBsimDataSKASrlz*.mat']);
fileNames = sort_nat({files.name});
Np = 1000; % number of pulsars
Nrlz = length(fileNames);

a=zeros(1,Np*(Np-1));% there are 17 pulsars in use so there are 16*17/2 coefficients
b=zeros(1,Np*(Np-1));
binSize = 1; % bin size
bins = 0:binSize:180;
bin_avg_CE = zeros(1,3); % averaged correlation coefficients over each bin

parfor i = 1:1:Nrlz % calculate 100 times
    [CE,thetaC,ct]=CorrelationCoefficient(simDataDir, fileNames{i});
    a(i,:) = CE;
    b(i,:)=thetaC;
end

AC=mean(a); % averaged cross-correlation coefficients over realizations
AT=mean(b);

% average over bins
parfor bin = 0:binSize:180
    index = find(AT>bin & AT<bin+binSize);
    bin_avg_CE(bin+1) = mean(AC(index));
end

% average for each realization
bin_avg_CE_prlz = zeros(1,3);
for n = 1:Nrlz
    for bin = 0:binSize:180
        idx = find(b(n,:) > bin & b(n,:) < bin+binSize );
        bin_avg_CE_prlz(n,bin+1) = mean(a(n,idx));
    end
end

% calculate the variance
for bin = 0:binSize:180
    var(bin+1) = sum((bin_avg_CE_prlz(:,bin+1)-bin_avg_CE(bin+1)).^2)/(length(bin_avg_CE)-1);
end

std = sqrt(var);

% Analytical
theta = linspace(0,180,1000);
x = (1-cos(theta*pi/180))/2;
y = 3/2 .* x .* log(x) - 1/4 .* x + 1/2;

% calculate the derivative
drvtv = diff(bin_avg_CE); % numerical
% convert bins from degree to rad
bins_rad = bins * pi/180;
drvtv_aly = 5/8 .* sin(bins_rad) + 3/4 * log((1-cos(bins_rad))/2) .* sin(bins_rad); % analytical

%% Plot
if isfolder([simDataDir,filesep,'fig']) == 0
    mkdir([simDataDir,filesep,'fig'])
end

% load([simDataDir,filesep,'fig',filesep,'HDC.mat'])

% figure 1 HD curve
figure(1)
% fa = 0.5/bin_avg_CE(1); % normalization factor
hold on
for n = 1:Nrlz
%     fa_rlz = 0.5/bin_avg_CE_prlz(n,1);
%     lines = plot(bins, fa_rlz*bin_avg_CE_prlz(n,:), '-', 'color',...
%     '#C1BEBE '); % with calibration factor fa_rlz
    lines = plot(bins, bin_avg_CE_prlz(n,:), '-', 'color', '#C1BEBE');
end
% avg_fig = plot(bins,fa*bin_avg_CE,'r-');
avg_fig = plot(bins,bin_avg_CE,'r-');
theo_fig = plot(theta,y,'b--');
% var_fig = shadedErrorBar(bins, bin_avg_CE, std,'LineProps','r');
boundedline(bins,bin_avg_CE, std, 'r');
hold off
legend([lines, avg_fig, theo_fig],{'Realizations','bin averaged','theoretical'})


exportgraphics(gcf,[simDataDir,filesep,'fig',filesep,'Hellings-Downs-Curve-Var.png'], 'Resolution',300)
savefig([simDataDir,filesep,'fig',filesep,'Hellings-Downs-Curve-Var'])

% figure 2 Standard deviation vs. derivative of HD curve
figure(2)
plot(drvtv_aly(1:180), var(1:180),'r+') % exclude the last point since it is NAN
xlabel('Derivative')
ylabel('Variance')
title('var vs. derivative')
exportgraphics(gcf,[simDataDir,filesep,'fig',filesep,'var_vs_derivative_aly.png'], 'Resolution',300)
savefig([simDataDir,filesep,'fig',filesep,'var_vs_derivative_aly'])

figure(3)
plot(drvtv(1:180), var(1:180),'r+') % exclude the last point since it is NAN
xlabel('Derivative')
ylabel('Variance')
title('var vs. derivative')
exportgraphics(gcf,[simDataDir,filesep,'fig',filesep,'var_vs_derivative.png'], 'Resolution',300)
savefig([simDataDir,filesep,'fig',filesep,'var_vs_derivative'])


%% save workspace
save([simDataDir,filesep,'fig',filesep,'HDC.mat'])




toc


%% Save


% END OF SCRIPT