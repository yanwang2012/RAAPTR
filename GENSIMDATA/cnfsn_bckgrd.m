% Confusion background by removing sources above certain SNR threshold

% QYQ 2022/01/13
%%
FileDir = '/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Final/realizations/2bands/simData/Band_opt_diff/GWBsimDataSKASrlz1Nrlz1.mat';
load(FileDir);
SNR_threshold = 10;
index = snr_chr > SNR_threshold;
Ns = sum(index);

delta_tmp = delta(index);
alpha_tmp = alpha(index);
phi0_tmp = phi0(index);
omega_tmp = omega(index);
Amp_tmp = Amp(index);
iota_tmp = iota(index);
thetaN_tmp = thetaN(index);

timingResiduals_tmp1 = zeros(simParams.Np,simParams.N);
for i=1:1:simParams.Np  % number of pulsar
    for j=1:1:Ns  % number of GW source

        % GW sky location in Cartesian coordinate
        k=zeros(1,3);  % unit vector pointing from SSB to source
        k(1)=cos(delta_tmp(j))*cos(alpha_tmp(j));
        k(2)=cos(delta_tmp(j))*sin(alpha_tmp(j));
        k(3)=sin(delta_tmp(j));
        theta=acos(k*simParams.kp(i,:)');
        %sprintf('%d pulsar theta=%g',i,theta)
        %phiI(i)=mod(phi0-omega*distP(i)*(1-cos(theta)), 2*pi);  % modulus after division
        %phiI(i)=mod(2*phi0-omega_tmp(l)*distP(i)*(1-cos(theta)), pi);  % modulus after division, YW 09/10/13
        phiI(i)=mod(phi0_tmp(j)-0.5*omega_tmp(j)*simParams.distP(i)*(1-cos(theta)), pi);  % modulus after division, YW 04/30/14 check original def. of phiI

        %disp(['pulsar = ', num2str(i), ' ', num2str(phiI(i))])

        tmp = FullResiduals(alpha_tmp(j),delta_tmp(j),omega_tmp(j),phi0_tmp(j),phiI(i),simParams.alphaP(i),simParams.deltaP(i),...
            Amp_tmp(j),iota_tmp(j),thetaN_tmp(j),theta,yr);

        timingResiduals_tmp1(i,:) = timingResiduals_tmp1(i,:) + tmp';

        %fftsignal(i,:)=fft(timingResiduals_tmp(i,:));

        % calculate the perfect fitness value

    end

end

timingResiduals_cnfus = timingResiduals_tmp - timingResiduals_tmp1; % remnant background.
obs = 1:simParams.N; % create observations vector.

%% Plot
figure
subplot(2,1,1);
plot(timingResiduals_cnfus(10,:));
subplot(2,1,2);
plot(timingResiduals_tmp(10,:));
saveas(gcf,'confusion_background.png');

figure
histogram(timingResiduals_cnfus);
saveas(gcf,'confusion_background_hist.png');

figure 
subplot(3,2,1)
plot(timingResiduals_cnfus(1,:));
% plot(obs,timingResiduals_cnfus(1,:),'-r',obs,noise(1,:),'-b');
xlabel('Observations')
ylabel('Residuals[s]')
legend('Signal','Noise')
subplot(3,2,3);
plot(timingResiduals_cnfus(5,:));
% plot(obs,timingResiduals_cnfus(5,:),'-r',obs,noise(5,:),'-b');
xlabel('Observations')
ylabel('Residuals[s]')
subplot(3,2,5)
plot(timingResiduals_cnfus(10,:));
% plot(obs,timingResiduals_cnfus(10,:),'-r',obs,noise(10,:),'-b')
xlabel('Observations')
ylabel('Residuals[s]')
subplot(3,2,[2,4,6])
histogram(timingResiduals_cnfus);
ylabel('Counts')
xlabel('Residuals')
saveas(gcf,'cnfus_bckgrd_tgthr.png');

