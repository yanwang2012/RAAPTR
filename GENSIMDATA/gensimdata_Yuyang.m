% Adopted from gensimdata.m to reproduce Yuyang's simulation results.
% Authors: QYQ
% 1/6/2021
clear;
tic
%%
%Get constants
dy2yr = genptaconsts('dy2yr');

% master directory for simulated data
simDataDir = '/Users/qyq/Research/PulsarTiming/YuYang_data';
NNs = 100; % 100 GW sources
Nrlz = 3; % # of noise realizations
Nnis = 3; % # of H0 data
rng(1) % fix seed for reproducibility

%% set source parameters using Yuyang's data
sim_source = load('/Users/qyq/Research/PulsarTiming/YuYang_data/sim_source.txt');
Amp_tmp = 10.^sim_source(:,6);
alpha_tmp = sim_source(:,1);
delta_tmp = asin(sim_source(:,2));
omega_tmp = sim_source(:,7);
iota_tmp = acos(sim_source(:,3));
thetaN_tmp = sim_source(:,4);
phi0_tmp = sim_source(:,5);
%% ==== Constructing a pulsar timing array using Np pulsars ====
% read in the pulsar catalogue simulated for SKA
skamsp=load('survey_ska.mat');% load input data
[~,I]=sort(skamsp.D);

Np = 100; % first 100 msp
N = 130; % 5 yrs biweekly
start=53187;  % Modified Julian Day, 'July 1, 2004'
deltaT=14;  % observation cadence, in days, set biweekly
load('searchParams_GWBsimDataSKA'); % load xmaxmin

% sky location of the pulsars in the equatorial coordinate
% we need to transfer from hr angle and degree to radian
alphaP=zeros(Np,1);  % right ascension, in radian
deltaP=zeros(Np,1);  % declination, in radian
distP=zeros(Np,1);  % (parallax) distance from SSB to pulsars, from mas/pc to ly
kp=zeros(Np,3);  % unit vector pointing from SSB to pulsars,
sd=zeros(Np,1);  % standard deviation of noise for different pulsar
dy=zeros(N,1);  % observation epoch, in day
yr=zeros(N,1);  % observation epoch, in year


%%
InList=zeros(1,2);
OutList=zeros(1,2);
for i=1:1:Np
    %     if OutList(i,1)>0
    %         alphaP(i)=OutList(i,1);  % in rad
    %     else
    %         alphaP(i)=OutList(i,1)+2*pi;  % in rad
    %     end
    
    InList(1,1)=skamsp.l(I(i));
    InList(1,2)=skamsp.b(I(i));
    
    % transfer the galactic coord. to equatorial coord.
    %[OutList,TotRot]=coco(InList,'g','j2000.0','d','r');
    [OutList,~]=coco(InList,'g','j2000.0','d','r');
    
    alphaP(i)=OutList(1,1);  % in rad
    deltaP(i)=OutList(1,2);  % in rad
    distP(i)=skamsp.D(I(i));  %0.28*kilo*pc2ly;  % in ly
    sd(i)=1.0*10^(-7);  % 100 ns
    
end

% % -------------------------------
% % transfer (hr,min,sec) and (degree,min,sec) to radian; mas/kpc to ly
% tmp1='J0030+0451';
% alphaP(1)=(0*15+30*15/60)*pi/180;
% deltaP(1)=(4+51/60)*pi/180;
% distP(1)=0.28*kilo*pc2ly;  % in ly
% sd(1)=79.2*10^(-8); %0.148;
%
% pname={tmp1 tmp2 tmp3 tmp4 tmp5 tmp6 tmp7 tmp8 tmp9};c
tmp1='place_holder';
pname={tmp1};

% sky location of pulsars in Cartesian coordinate
for i=1:1:Np
    kp(i,1)=cos(deltaP(i))*cos(alphaP(i));
    kp(i,2)=cos(deltaP(i))*sin(alphaP(i));
    kp(i,3)=sin(deltaP(i));
end

for i=1:1:N
    dy(i)=start+(i-1)*deltaT;  % dates conducting observations, in MJD
    yr(i)=2004.5+(i-1)*deltaT*dy2yr;
end

simParams = struct('Np',Np,'N',N,'sd',sd,...
    'alphaP',alphaP,'deltaP',deltaP,'distP',distP,'kp',kp);

% -------------------------------
% calculate timing residuals induced by GW for each pulsar
%Amp=1.3*10^(-7);  % overall amplitude of timing residuals, the same for all pulsars, sec
timingResiduals=zeros(Np,N);  % signal, i.e. GW induced timing residuals, Np pulsars, N observations
timingResiduals_tmp=zeros(Np,N);   % signal without noise
phiI=zeros(Np,1);  % arbitrary phase for each pulsar, relative distance

noise=zeros(Np,N);  % noise
%sd=0.1*10^(-7);  % standard deviation of the normal distribtion (sec)

CA4filenames=cell(Nrlz+Nnis,1);  % cell array for file names of simulated data

perfect_fitness=zeros(NNs,1);  % fitness value for the true parameters

%Standardized true parameter values
%stdTrueCoord = zeros(1,12);  % 4+Np
stdTrueCoord = zeros(NNs,7); % 7 parameters other than pulsar phases

% calculate SNR
snr_chr2_tmp=zeros(Np,NNs);  % squared characteristic snr for each pulsar and source
%snr_chr=zeros(Np,1);  % characteristic srn = <\rho> = sqrt[(h,h)] srn of signal -- the strenghth of signal
%snr_mf=0;  % snr of matched filter = \rho = x/\sigma, normalized matched filter

genHypothesis='H1 data';

tmp=zeros(1,N); % noiseless timing residuals from a source
for i=1:1:Np  % number of pulsar
    for j=1:1:NNs  % number of GW source
        
        % GW sky location in Cartesian coordinate
        k=zeros(1,3);  % unit vector pointing from SSB to source
        k(1)=cos(delta_tmp(j))*cos(alpha_tmp(j));
        k(2)=cos(delta_tmp(j))*sin(alpha_tmp(j));
        k(3)=sin(delta_tmp(j));
        theta=acos(k*kp(i,:)');
        %sprintf('%d pulsar theta=%g',i,theta)
        %phiI(i)=mod(phi0-omega*distP(i)*(1-cos(theta)), 2*pi);  % modulus after division
        %phiI(i)=mod(2*phi0-omega_tmp(l)*distP(i)*(1-cos(theta)), pi);  % modulus after division, YW 09/10/13
        phiI(i)=mod(phi0_tmp(j)-0.5*omega_tmp(j)*distP(i)*(1-cos(theta)), pi);  % modulus after division, YW 04/30/14 check original def. of phiI
        
        %disp(['pulsar = ', num2str(i), ' ', num2str(phiI(i))])
        
        tmp = FullResiduals(alpha_tmp(j),delta_tmp(j),omega_tmp(j),phi0_tmp(j),phiI(i),alphaP(i),deltaP(i),...
            Amp_tmp(j),iota_tmp(j),thetaN_tmp(j),theta,yr);
        
        timingResiduals_tmp(i,:) = timingResiduals_tmp(i,:) + tmp';
        
        %fftsignal(i,:)=fft(timingResiduals_tmp(i,:));
        
        % calculate the perfect fitness value
        
        snr_chr2_tmp(i,j) = dot(tmp,tmp)/sd(i)^2;
        
        % standardization of the true coordinates
        stdTrueCoord(j,1)=(alpha_tmp(j)-xmaxmin(1,2))/(xmaxmin(1,1)-xmaxmin(1,2));  % [0, 2*pi]
        stdTrueCoord(j,2)=(delta_tmp(j)-xmaxmin(2,2))/(xmaxmin(2,1)-xmaxmin(2,2));  % [-pi/2, pi/2]
        stdTrueCoord(j,3)=(omega_tmp(j)-xmaxmin(3,2))/(xmaxmin(3,1)-xmaxmin(3,2));  % [2, 20]
        stdTrueCoord(j,4)= mod(phi0_tmp(j),pi)/pi;  % [0, pi]
        stdTrueCoord(j,5)=(log10(Amp_tmp(j))-xmaxmin(5,2))/(xmaxmin(5,1)-xmaxmin(5,2));
        stdTrueCoord(j,6)=(iota_tmp(j)-xmaxmin(6,2))/(xmaxmin(6,1)-xmaxmin(6,2));
        stdTrueCoord(j,7)=(thetaN_tmp(j)-xmaxmin(7,2))/(xmaxmin(7,1)-xmaxmin(7,2));
        
    end
    
end

snr_chr=sqrt(sum(snr_chr2_tmp,1));  % sum of elements in each column

% signal + noise realizations
nf=0;  % counter for number of files
for jj=1:1:Nrlz
    
    nf=nf+1;
    
    rlz_id=jj; %num2str(jj);
    
    % structure to store the id tag for each metadata file
    snr_id=0;  % not used for GWB
    loc_id=0;
    omg_id=0;
    id=struct('snr_id',snr_id,'loc_id',loc_id,'omg_id',omg_id,'rlz_id',rlz_id);
    
    for i=1:1:Np
        
        % generating a realization of noise
        noise(i,:)=sd(i)*randn(1,N);  % Gaussian noise
        % calculate the actual snr
        %fftnoise(i,:)=fft(noise(i,:));
        
        %snr_tmp=snr_tmp+ fftsignal(i,:)*fftsignal(i,:)/fftnoise(i,:);
        timingResiduals(i,:)=timingResiduals_tmp(i,:)+noise(i,:);  % add noise on signal
        
    end
    
    inParams = struct('Np',Np,'N',N,'s',timingResiduals,'sd',sd,...
        'alphaP',alphaP,'deltaP',deltaP,'kp',kp,'yr',yr,...
        'xmaxmin',xmaxmin);
    
    for j=1:1:NNs
        perfect_fitness(j) = LLR_PSOmpp(stdTrueCoord(j,:),inParams);  % - LogLikelihoodRatio, minimization
        %true_fitness = fitnessTrue_ie(alpha_tmp(j),delta_tmp(j),omega_tmp(l),Amp,iota,thetaN,phi0,phiI,inParams);
    end
    
    %disp(['In simulator2: perfect_fitness: ', num2str(perfect_fitness)]);
    
    % save metadata into a file for each realization (file name rule)
    %filename=strcat('snr',num2str(ii),'loc',num2str(j),'omg',num2str(l),'rlz',num2str(jj),'.mat');
    filename=strcat('GWBsimDataSKA','Srlz',num2str(1),'Nrlz',num2str(jj),'.mat'); % change Srlz when changing rng seed.
    snr=0;  % this variable is not useful here
    alpha=alpha_tmp;
    delta=delta_tmp;
    omega=omega_tmp;
    Amp = Amp_tmp;
    iota = iota_tmp;
    thetaN = thetaN_tmp;
    phi0 = phi0_tmp;
    CA4filenames{nf}=filename;
    
    save([simDataDir,filesep,filename],'genHypothesis','snr','alpha','delta','omega','iota','thetaN','phi0',...
        'timingResiduals','noise','yr', 'pname','id','simParams','perfect_fitness','snr_chr','timingResiduals_tmp','Amp');
    
end

% save the GW source parameters into a file, have duplication with data files
% filename=strcat('GWB_','Srlz',num2str(1),'.mat');
% save([simDataDir,filesep,filename],'AmpOut','alphaOut','deltaOut',...
%     'omega_tmp','fgwOut','iotaOut','thetaNOut','phi0Out','r','stdTrueCoord');

% simulating noise without signal
%Nnis=10; %100;  % number of realization of noise
perfect_fitness = 0;
snr_chr = 0;
timingResiduals_tmp=0;
for ii=1:1:Nnis
    
    genHypothesis='H0 data';
    snr_id=0; %'0';
    loc_id=0; %'0';
    omg_id=0; %'0';
    rlz_id=ii; %num2str(ii);
    
    id=struct('snr_id',snr_id,'loc_id',loc_id,'omg_id',omg_id,'rlz_id',rlz_id);
    
    for i=1:1:Np
        
        % generating a realization of noise
        noise(i,:)=sd(i)*randn(1,N);  % Gaussian noise
        timingResiduals(i,:)=noise(i,:);
        %timingResiduals(i,:)=timingResiduals_tmp(i,:)+noise(i,:);  % add noise on signal
        
    end
    
    % save metadata into a file for each realization (file name rule)
    filename=strcat('noise',num2str(ii),'.mat');
    snr=0.0;  %snr_tmp(ii);
    alpha=0.0;  %alpha_tmp(j);
    delta=0.0;  %delta_tmp(j);
    omega=0.0;  %omega_tmp(l);
    iota=0.0;
    thetaN=0.0;
    phi0=0.0;
    save([simDataDir,filesep,filename],'genHypothesis','snr','alpha','delta','omega','iota','thetaN','phi0',...
        'timingResiduals','noise','yr','pname','id','simParams','perfect_fitness','snr_chr','timingResiduals_tmp','Amp');
    CA4filenames{nf+ii}=filename;
    
end
toc
% end of script
