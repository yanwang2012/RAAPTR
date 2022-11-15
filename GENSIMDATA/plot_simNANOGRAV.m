% Plot Simulated and Estimated sources using NANOGrav 12.5 year pulsars.

% Author: Yiqian Qian
% Date: 11/15/2022

clear;
%% calculate SNR
%phiI = psrParams.phiI;
DataDir = '/Users/yiqianqian/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/SimNANOGrav';
ResultsDir = '/Users/yiqianqian/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/SimNANOGrav/Results';
resultFiles = dir([ResultsDir,filesep,'*.mat']);
resultFileNames = {resultFiles.name};
Ns = length(resultFileNames); % Number of estimated sources.
psrFile = dir([DataDir,filesep,'NANOGrav12.5yv4.mat']);

psrParams = load([DataDir,filesep,psrFile.name]);

Np = psrParams.Np;% number of pulsars
N = psrParams.N;% observation period 5 yr biweekly
snr_chr2_tmp=zeros(Np,1);  % squared characteristic snr for each pulsar and source

% pulsar parameters
alphaP = psrParams.alphaP;
deltaP = psrParams.deltaP;
sd = psrParams.sd;
yr = psrParams.yr;

% sky location of pulsars in Cartesian coordinate
for i=1:1:Np
    kp(i,1)=cos(deltaP(i))*cos(alphaP(i));
    kp(i,2)=cos(deltaP(i))*sin(alphaP(i));
    kp(i,3)=sin(deltaP(i));
end

omega_est = zeros(Ns,1);
snr_est = zeros(Ns,1);
for src = 1:1:Ns
    srcParams = load([ResultsDir,filesep,resultFileNames{src}]);
    % asign source parameters
    alpha = srcParams.bestRealLoc(1);
    delta = srcParams.bestRealLoc(2);
    omega = srcParams.bestRealLoc(3); % need to check the unit
    omega_est(src) = omega; % save frequency for later use
    phi0 = srcParams.bestRealLoc(4);
    Amp = srcParams.bestRealLoc(5);
    iota = srcParams.bestRealLoc(6);
    thetaN = srcParams.bestRealLoc(7);
    phiI = srcParams.bestRealLoc(8:end);

    timingResiduals = {};
    for i=1:1:Np  % number of pulsar
        % GW sky location in Cartesian coordinate
        k=zeros(1,3);  % unit vector pointing from SSB to source
        k(1)=cos(delta)*cos(alpha);
        k(2)=cos(delta)*sin(alpha);
        k(3)=sin(delta);
        theta=acos(k*kp(i,:)');
        %phiI(j)=mod(srcParams.phi0-0.5*srcParams.omega*psrParams.distP(j)*(1-cos(theta)), pi);  % modulus after division, YW 04/30/14 check original def. of phiI


        tmp = FullResiduals(alpha,delta,omega,phi0,phiI(i),alphaP(i),...
            deltaP(i),Amp,iota,thetaN,theta,yr{i});

        timingResiduals{i} = tmp';
        snr_chr2_tmp(i,1) = sum((tmp'.*tmp')./(sd{i}.^2)); % SNR^2 for each pulsar

    end

    snr_chr=sqrt(sum(snr_chr2_tmp,1));  % sum of elements in each column
    snr_est(src) = snr_chr; % save for later use
end

%% Plot
simParams = load([DataDir, filesep,'GWB_Srlz1.mat']);
sim_omega = simParams.omega;
sim_snr = simParams.snr_chr;
cnst = 2*pi*365*24*3600; % constant to convert rad yr^-1 to Hz

figure(1)
plot(sim_snr, sim_omega/cnst,'ob', snr_est, omega_est/cnst, 'r*')
xlabel('SNR')
ylabel('Frequency [Hz]')
title('Simulate vs. Estimate')
legend('simulated', 'estimated')
savefig([ResultsDir,filesep,'plot'])
saveas(gcf, [ResultsDir, filesep, 'plot'], 'png')

% EOS