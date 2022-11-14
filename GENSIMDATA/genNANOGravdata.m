%function []=genNANOGravdata(path_to_parameters,path_to_pulsar_catalog,path_to_output,frqRng)
%Generate PTA data realizations with multiple SMBHB sources using NANOGrav
%pulsar catalog. Adopted from gensimdata.

%Nov 2022
clear;
%Get constants
dy2yr = genptaconsts('dy2yr');

outDir = '/Users/yiqianqian/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/simDATA/SimNANOGrav';
if isfolder(outDir)== 0
    mkdir(outDir);
end

Ns = 100; % number of sources.
Np = 47; % number of pulsars
Nrlz = 1; % number of H1 realization with different noise realizations
Nsrlz = 1; % number of H1 realization with different src realizations
Nnis = 1; % number of H0 realization
parameters(Ns,Np,Nrlz,Nsrlz,Nnis,outDir);
path_to_parameters = [outDir,filesep,'GWBsimDataSKAParams.mat'];
psr_params = load(path_to_parameters);% load all the constants and parameters

start = psr_params.start;
deltaT = psr_params.deltaT;
xmaxmin = psr_params.xmaxmin;

% master directory for simulated data
simDataDir = outDir;

%rng('shuffle')  % initialize the random number generator using a different seed
rng(1105) % for repeatable works.

%% ==== Generate random GW sources ====
for srlz = 1:Nsrlz
    [AmpOut,alphaOut,deltaOut,fgwOut,iotaOut,thetaNOut,phi0Out,r]=GenerateRandomGWSource(Ns);
    % put check for frequency above Nyquist freq
    cadence = 14; % biweekly;
    Nyquist_freq = 1/2 * 1/(cadence*24*3600); % Nyquist frequency (Hz)

    %     % setting up the freq filter
    %     if nargin == 3
    disp("Searching whole band.");
    if max(fgwOut) > Nyquist_freq
        disp("Warning: Frequency above Nyquist frequency:" +Nyquist_freq)
        disp("Will exclude all sources above Nyquist frequency")
    end
    Index = (fgwOut <= Nyquist_freq);
    % disp(['Index is: ', num2str(Index)]);
    if sum(Index) == 0
        disp('There is no source');
        continue
    elseif sum(Index) > 0
        Amp_tmp = AmpOut(Index);
        alpha_tmp = alphaOut(Index);
        delta_tmp = deltaOut(Index);
        iota_tmp = iotaOut(Index);
        thetaN_tmp = thetaNOut(Index);
        phi0_tmp = phi0Out(Index);
        omega_tmp = fgwOut(Index)*2*pi*24*365*3600; % convert sec^-1 (Hz) to rad yr^-1
        NNs = sum(Index);
        %     disp("The number of sources is:" + NNs)
    end
    %     elseif nargin == 4
    %         fgwl = frqRng(1);
    %         fgwh = frqRng(2);
    %         disp("The lower limit is:" + fgwl);
    %         disp("The upper limit is:" + fgwh);
    %         Index = (fgwOut >= fgwl & fgwOut <= fgwh);% select the sources within bands
    %         Amp_tmp = AmpOut(Index);
    %         omega_tmp = 2*pi * fgwOut(Index) * 3.156*10^7; % unit change to yr^-1
    %         alpha_tmp = alphaOut(Index);
    %         delta_tmp = deltaOut(Index);
    %         iota_tmp = iotaOut(Index);
    %         thetaN_tmp = thetaNOut(Index);
    %         phi0_tmp = phi0Out(Index);
    %         NNs = sum(Index);
    %     end
    %
    %     disp("The number of sources is: " + NNs);
    %     disp("The size of Amp_tmp is: " + length(Amp_tmp));
    %disp("m is:" + m);
    %% ==== Constructing a pulsar timing array using Np pulsars ====
    % read in the pulsar catalogue simulated for SKA
    path_to_pulsar_catalog = 'NANOGrav12.5yv4.mat';
    nngpsr=load(path_to_pulsar_catalog);% load input data

    % asign NANOGrav pusar data
    alphaP = nngpsr.alphaP;
    deltaP = nngpsr.deltaP;
    distP = nngpsr.distP;
    yr = nngpsr.yr;
    NN = nngpsr.N;
    sd = nngpsr.sd;
    psr_catalog = nngpsr.Names;

    kp=zeros(Np,3);  % unit vector pointing from SSB to pulsars,

    % convert kpc to ly
    kilo=1000;
    pc2ly=3.26;
    distP = distP*kilo*pc2ly;

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

    % -------------------------------
    % calculate timing residuals induced by GW for each pulsar
    %Amp=1.3*10^(-7);  % overall amplitude of timing residuals, the same for all pulsars, sec


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

    phiI=zeros(Np,1);  % arbitrary phase for each pulsar, relative distance

    psr_data = {}; % struct cell to store different pulsar data.
    noises = {}; % cell array for different pulsar noises.

    % signal + noise realizations
    nf = 0;

    for jj = 1:1:Nrlz
        
        nf = nf + 1;

        for i=1:1:Np  % number of pulsar
            N = NN(i);
            noise=zeros(1,N);  % noise
            tmp=zeros(1,N); % noiseless timing residuals from a source
            timingResiduals=zeros(1,N);  % signal, i.e. GW induced timing residuals, Np pulsars, N observations
            timingResiduals_tmp=zeros(1,N);   % signal without noise
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
                    Amp_tmp(j),iota_tmp(j),thetaN_tmp(j),theta,yr{i});

                timingResiduals_tmp = timingResiduals_tmp + tmp';

                %fftsignal(i,:)=fft(timingResiduals_tmp(i,:));

                % calculate the perfect fitness value

                snr_chr2_tmp(i,j) = sum((tmp'.*tmp')./(sd{i}.^2));

                % standardization of the true coordinates
                stdTrueCoord(j,1)=(alpha_tmp(j)-xmaxmin(1,2))/(xmaxmin(1,1)-xmaxmin(1,2));  % [0, 2*pi]
                stdTrueCoord(j,2)=(delta_tmp(j)-xmaxmin(2,2))/(xmaxmin(2,1)-xmaxmin(2,2));  % [-pi/2, pi/2]
                stdTrueCoord(j,3)=(omega_tmp(j)-xmaxmin(3,2))/(xmaxmin(3,1)-xmaxmin(3,2));  % [2, 20]
                stdTrueCoord(j,4)= mod(phi0_tmp(j),pi)/pi;  % [0, pi]
                stdTrueCoord(j,5)=(log10(Amp_tmp(j))-xmaxmin(5,2))/(xmaxmin(5,1)-xmaxmin(5,2));
                stdTrueCoord(j,6)=(iota_tmp(j)-xmaxmin(6,2))/(xmaxmin(6,1)-xmaxmin(6,2));
                stdTrueCoord(j,7)=(thetaN_tmp(j)-xmaxmin(7,2))/(xmaxmin(7,1)-xmaxmin(7,2));

            end

            snr_chr=sqrt(sum(snr_chr2_tmp,1));  % sum of elements in each column

            % generate noise
            for j = 1:1:N
                noise(j) = sd{i}(j) * randn(1);
            end
            timingResiduals = timingResiduals_tmp + noise;
            noises{i} = noise;
            psr_data{i} = struct('psr_name', psr_catalog{i}, 'alphaP', alphaP(i), 'deltaP',deltaP(i), 'yr', yr{i}, 'sd', sd{i}, 'timingResiduals', timingResiduals,...
                                    'timingResiduals_tmp', timingResiduals_tmp, 'N', N);
        end

        filename=strcat('GWBsimDataSKA','Srlz',num2str(srlz),'Nrlz',num2str(jj),'.mat'); % change Srlz when changing rng seed.
        snr=0;  % this variable is not useful here
        alpha=alpha_tmp;
        delta=delta_tmp;
        omega=omega_tmp;
        Amp = Amp_tmp;
        iota = iota_tmp;
        thetaN = thetaN_tmp;
        phi0 = phi0_tmp;
        CA4filenames{nf}=filename;

%         save([simDataDir,filesep,filename],'genHypothesis','snr','alpha','delta','omega','iota','thetaN','phi0',...
%             'timingResiduals','noise','yr', 'pname','id','simParams','perfect_fitness','snr_chr','timingResiduals_tmp','Amp');
        save([simDataDir, filesep, filename], 'psr_data', 'noises');
    end

    % signal + noise realizations
    %     nf=0;  % counter for number of files
    %     for jj=1:1:Nrlz
    %
    %         nf=nf+1;
    %
    %         rlz_id=jj; %num2str(jj);
    %
    %         % structure to store the id tag for each metadata file
    %         snr_id=0;  % not used for GWB
    %         loc_id=0;
    %         omg_id=0;
    %         id=struct('snr_id',snr_id,'loc_id',loc_id,'omg_id',omg_id,'rlz_id',rlz_id);
    %
    %         % generating a realization of noise
    %         % noise(i,:)=sd(i)*randn(1,N);  % Gaussian noise
    %
    %         % now sd is different for different time
    %         for j = 1:1:N
    %             noise(j) = sd{i}(j) * randn(1);
    %         end
    %         % calculate the actual snr
    %         %fftnoise(i,:)=fft(noise(i,:));
    %
    %         %snr_tmp=snr_tmp+ fftsignal(i,:)*fftsignal(i,:)/fftnoise(i,:);
    %         timingResiduals=timingResiduals_tmp+noise;  % add noise on signal
    %
    %         inParams = struct('Np',Np,'N',N,'s',timingResiduals,'sd',sd,...
    %             'alphaP',alphaP,'deltaP',deltaP,'kp',kp,'yr',yr,...
    %             'xmaxmin',xmaxmin);
    %
    %         for j=1:1:NNs
    %             perfect_fitness(j) = LLR_PSOmpp(stdTrueCoord(j,:),inParams);  % - LogLikelihoodRatio, minimization
    %             %true_fitness = fitnessTrue_ie(alpha_tmp(j),delta_tmp(j),omega_tmp(l),Amp,iota,thetaN,phi0,phiI,inParams);
    %         end

    %disp(['In simulator2: perfect_fitness: ', num2str(perfect_fitness)]);

    % save metadata into a file for each realization (file name rule)
    %filename=strcat('snr',num2str(ii),'loc',num2str(j),'omg',num2str(l),'rlz',num2str(jj),'.mat');
    %     filename=strcat('GWBsimDataSKA','Srlz',num2str(srlz),'Nrlz',num2str(jj),'.mat'); % change Srlz when changing rng seed.
    %     snr=0;  % this variable is not useful here
    %     alpha=alpha_tmp;
    %     delta=delta_tmp;
    %     omega=omega_tmp;
    %     Amp = Amp_tmp;
    %     iota = iota_tmp;
    %     thetaN = thetaN_tmp;
    %     phi0 = phi0_tmp;
    %     CA4filenames{nf}=filename;
    %
    %     save([simDataDir,filesep,filename],'genHypothesis','snr','alpha','delta','omega','iota','thetaN','phi0',...
    %         'timingResiduals','noise','yr', 'pname','id','simParams','perfect_fitness','snr_chr','timingResiduals_tmp','Amp');

    %     end

    % save the GW source parameters into a file, have duplication with data files
    filename=strcat('GWB_','Srlz',num2str(srlz),'.mat');
    save([simDataDir,filesep,filename],'Amp','alpha','delta',...
        'omega','iota','thetaN','phi0','r','stdTrueCoord');

    % simulating noise without signal
    %Nnis=10; %100;  % number of realization of noise
%     perfect_fitness = 0;
%     snr_chr = 0;
%     timingResiduals_tmp=0;
%     for ii=1:1:Nnis
% 
%         genHypothesis='H0 data';
%         snr_id=0; %'0';
%         loc_id=0; %'0';
%         omg_id=0; %'0';
%         rlz_id=ii; %num2str(ii);
% 
%         id=struct('snr_id',snr_id,'loc_id',loc_id,'omg_id',omg_id,'rlz_id',rlz_id);
% 
%         for i=1:1:Np
% 
%             % generating a realization of noise
%             noise(i,:)=sd(i)*randn(1,N);  % Gaussian noise
%             timingResiduals(i,:)=noise(i,:);
%             %timingResiduals(i,:)=timingResiduals_tmp(i,:)+noise(i,:);  % add noise on signal
% 
%         end
% 
%         % save metadata into a file for each realization (file name rule)
%         filename=strcat('srlz',num2str(srlz),'noise',num2str(ii),'.mat');
%         snr=0.0;  %snr_tmp(ii);
%         alpha=0.0;  %alpha_tmp(j);
%         delta=0.0;  %delta_tmp(j);
%         omega=0.0;  %omega_tmp(l);
%         iota=0.0;
%         thetaN=0.0;
%         phi0=0.0;
%         save([simDataDir,filesep,filename],'genHypothesis','snr','alpha','delta','omega','iota','thetaN','phi0',...
%             'timingResiduals','noise','yr','pname','id','simParams','perfect_fitness','snr_chr','timingResiduals_tmp','Amp');
%         CA4filenames{nf+ii}=filename;

%     end
end

% end of script
