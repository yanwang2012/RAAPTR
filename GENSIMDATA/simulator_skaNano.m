% It is a simulator based on simulator2.m, we make it to simulate the data 
% for the SKA pulsars. RA, dec, distP are from..., sd is calculated based
% on the assumptions used in R. Smits ...
% Aug. 19, 2015

clear;
tic;


% master directory for simulated data
simDataDir = 'simDataSKA_loc_J334_snr6';
%simDataDir = 'simDataSKA_loc_PG1302';
searchParamsFile = 'searchParams_simDataSKA_X';
mkdir(simDataDir);

rng('shuffle')  % initialize the random number generator using a different seed

% ======  Useful constants  ======
pc2ly=3.261563777;  % 1 pc=3.26 ly (Julian)
dy2yr=1.0/365.25;  % 1 day=365.25 yr (Julian)  ??
kilo=1.0*10^3;  % kilo 1000

% ==== GW source parameters ====
Nsnr=6;  % number of snr of matched filter  100, 30, 8
snr_net=zeros(Nsnr,1);  % network snr
% snr_net(1)=8.0*12.7; % loc=3
% snr_net(2)=30.0*12.7;
% snr_net(3)=100.0*12.7;

% snr_net(1)=8.0*12.7/1.34/85.48*1.14;  % loc=3
% snr_net(2)=30.0*12.7/1.34/85.48*1.14;
% snr_net(3)=100.0*12.7/1.34/85.48*1.14;

% snr_net(1)=30.0*12.7/1.34/85.48;  % loc=4
% snr_net(2)=60.0*12.7/1.34/85.48;
% snr_net(3)=100.0*12.7/1.34/85.48;

% snr_net(1)=30.0*12.7/1.34/85.48*1.14/1.28;  % loc=5
% snr_net(2)=60.0*12.7/1.34/85.48*1.14/1.28;
% snr_net(3)=100.0*12.7/1.34/85.48*1.14/1.28;

% snr_net(1)=30.0*12.7/1.34/85.48*1.14/1.58/0.873/1.16/0.94;  % loc=6
% snr_net(2)=60.0*12.7/1.34/85.48*1.14/1.58/0.873/1.16/0.94;
% snr_net(3)=100.0*12.7/1.34/85.48*1.14/1.58/0.873/1.16/0.94;

snr_net(1)=5.0/8.4/1.5;  % loc=6
snr_net(2)=8.0/8.4/1.5;
snr_net(3)=10.0/8.4/1.5;
snr_net(4)=12.0/8.4/1.5;
snr_net(5)=15.0/8.4;  %/1.5;
snr_net(6)=20.0/8.4;  %/1.5;

%snr_net(1)=8.0*12.7/1.34/1.54/85.48;  % loc=5
%snr_net(2)=30.0*12.7/1.34/1.54/85.48;
%snr_net(3)=100.0*12.7/1.34/1.54/85.48;

snr_tmp=zeros(Nsnr,1);
snr_tmp(1)=1.0/13.2242*snr_net(1); %10.0; %0.5; %8.0;
snr_tmp(2)=1.0/13.2242*snr_net(2); %10.0;
snr_tmp(3)=1.0/13.2242*snr_net(3);  %12.0;
snr_tmp(4)=1.0/13.2242*snr_net(4);
snr_tmp(5)=1.0/13.2242*snr_net(5);
snr_tmp(6)=1.0/13.2242*snr_net(6);

Nloc=6;  % number of source sky location
alpha_tmp=zeros(Nloc,1);
delta_tmp=zeros(Nloc,1);
alpha_tmp(1)=pi/4+1.2;
delta_tmp(1)=pi/4-0.16;
alpha_tmp(2)=4.367;  % cond(A)=1.017
delta_tmp(2)=0.8796;
% location used in Taylor et al.
alpha_tmp(3)=0.9999;
delta_tmp(3)=0.5007;  % pi/2 - acos(0.48)
% loc 4, NOT Taylor et al.
alpha_tmp(4)=0.0;
delta_tmp(4)=1.2;
% loc 4, NOT Taylor et al.
alpha_tmp(5)=0.0;
delta_tmp(5)=0.0;
alpha_tmp(6)=3.5;
delta_tmp(6)=0.3;

alpha_tmp(7)=3.5;
delta_tmp(7)=-0.2;
alpha_tmp(8)=3.5;
delta_tmp(8)=-0.7;
alpha_tmp(9)=3.5;
delta_tmp(9)=-1.2;

alpha_tmp(10)=0.9338; % J334
delta_tmp(10)=0.0246;
alpha_tmp(11)=3.4252; % PG1302
delta_tmp(11)=-0.1841;

Nomg=5;  % number of GW frequency
omega_tmp=zeros(Nomg,1);
omega_tmp(1)=2*pi/0.3925; %16.0081;  % 0.3925 yr  high freq
omega_tmp(2)=4.2743;  % 1.47 yr  low freq
omega_tmp(3)=2*pi/(1.6);  % orbital freq is 10^(-8) Hz

omega_tmp(4)=2*4.2454;  % Loc 10
omega_tmp(5)=2*1.21;  % Loc 11


% optimized value
iota=0.4949;  %0.0;  %pi/4+0.6;  % inclination between orbital plane and plane of the sky
thetaN=0.5;  %pi/2;  % angle to the line of nodes
phi0=2.89;  %0.0;  % initial orbital phase

% ==== Constructing a pulsar timing array using Np pulsars ====
% read in the pulsar catalogue simulated for SKA
skamsp=load('/Your/path/to/the/file/survey_ska.mat');
[~,I]=sort(skamsp.D);
Np=1000;  % number of pulsars in the timing array

% sky location of the pulsars in the equatorial coordinate
% we need to transfer from hr angle and degree to radian
alphaP=zeros(Np,1);  % right ascension, in radian
deltaP=zeros(Np,1);  % declination, in radian
distP=zeros(Np,1);  % (parallax) distance from SSB to pulsars, from mas/pc to ly
kp=zeros(Np,3);  % unit vector pointing from SSB to pulsars,
sd=zeros(Np,1);  % standard deviation of noise for different pulsar

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
    sd(i)=1.0*10^(-8);  %79.2*10^(-8); %0.148;
    
end


% -------------------------------
% transfer (hr,min,sec) and (degree,min,sec) to radian; mas/kpc to ly
%tmp1='J0030+0451';
alphaP(1)=(0*15+30*15/60)*pi/180;
deltaP(1)=(4+51/60)*pi/180;
distP(1)=1.376*kilo*pc2ly;  % in ly
sd(1)=1.0*10^(-8); %0.148;

%tmp2='J0613-0200';
alphaP(2)=(6*15+13*15/60)*pi/180;
deltaP(2)=-(2+0/60)*pi/180;
distP(2)=6.318*kilo*pc2ly;
sd(2)=1.0*10^(-8); %0.178;

%tmp3='J1713+0747';
alphaP(3)=(17*15+13*15/60)*pi/180;
deltaP(3)=(7+47/60)*pi/180;
distP(3)=7.524*kilo*pc2ly;
sd(3)=1.0*10^(-8); %0.03;

%tmp4='J1909-3744';
alphaP(4)=(19*15+9*15/60)*pi/180;
deltaP(4)=-(37+44/60)*pi/180;
distP(4)=3.532*kilo*pc2ly;
sd(4)=1.0*10^(-8); %0.038;

%-----------
%tmp5='J1012+5307';
alphaP(5)=(10*15+12*15/60)*pi/180;
deltaP(5)=(53+7/60)*pi/180;
distP(5)=1.045*kilo*pc2ly;
sd(5)=1.0*10^(-8); %0.276;

%tmp6='J1455-3330';
alphaP(6)=(14*15+55*15/60)*pi/180;
deltaP(6)=-(33+30/60)*pi/180;
distP(6)=6.593*kilo*pc2ly;
sd(6)=1.0*10^(-8); %0.787;

%tmp7='J1600-3053';
alphaP(7)=(16*15+0*15/60)*pi/180;
deltaP(7)=-(30+53/60)*pi/180;
distP(7)=13.532*kilo*pc2ly;
sd(7)=1.0*10^(-8); %0.163;

%tmp8='J1640+2224';
alphaP(8)=(16*15+40*15/60)*pi/180;
deltaP(8)=(22+24/60)*pi/180;
distP(8)=3.675*kilo*pc2ly;
sd(8)=1.0*10^(-8); %0.409;

% ---- add other nanograv pulsars  09/25/2014  YW
%tmp9='J1643-1224';
alphaP(9)=(16*15+43*15/60)*pi/180;
deltaP(9)=-(12+24/60)*pi/180;
distP(9)=2.735*kilo*pc2ly;
sd(9)=1.0*10^(-8);

%tmp10='J1744-1134';
alphaP(10)=(17*15+44*15/60)*pi/180;
deltaP(10)=-(11+34/60)*pi/180;
distP(10)=0.453*kilo*pc2ly;
sd(10)=1.0*10^(-8);

%tmp11='J1853+1308';
alphaP(11)=(18*15+53*15/60)*pi/180;
deltaP(11)=(13+8/60)*pi/180;
distP(11)=7.243*kilo*pc2ly;
sd(11)=1.0*10^(-8);

%tmp12='B1855+09';  % J1857+0943
alphaP(12)=(18*15+57*15/60)*pi/180;
deltaP(12)=(9+43/60)*pi/180;
distP(12)=3.239*kilo*pc2ly;
sd(12)=1.0*10^(-8);

%tmp13='J1910+1256';
alphaP(13)=(19*15+10*15/60)*pi/180;
deltaP(13)=(12+56/60)*pi/180;
distP(13)=13.542*kilo*pc2ly;
sd(13)=1.0*10^(-8);

%tmp14='J1918-0642';
alphaP(14)=(19*15+18*15/60)*pi/180;
deltaP(14)=-(6+42/60)*pi/180;
distP(14)=6.437*kilo*pc2ly;
sd(14)=1.0*10^(-8);

%tmp15='B1953+29';  % J1955+2908
alphaP(15)=(19*15+55*15/60)*pi/180;
deltaP(15)=(29+8/60)*pi/180;
distP(15)=0.875*kilo*pc2ly;
sd(15)=1.0*10^(-8);

%tmp16='J2145-0750';
alphaP(16)=(21*15+45*15/60)*pi/180;
deltaP(16)=-(7+50/60)*pi/180;
distP(16)=3.982*kilo*pc2ly;
sd(16)=1.0*10^(-8);

%tmp17='J2317+1439';
alphaP(17)=(23*15+17*15/60)*pi/180;
deltaP(17)=(14+39/60)*pi/180;
distP(17)=5.281*kilo*pc2ly;
sd(17)=1.0*10^(-8);

tmp1='place_holder';
pname={tmp1};

%%
% sky location of pulsars in Cartesian coordinate
for i=1:1:Np
    kp(i,1)=cos(deltaP(i))*cos(alphaP(i));
    kp(i,2)=cos(deltaP(i))*sin(alphaP(i));
    kp(i,3)=sin(deltaP(i));
end

% -------------------------------
% starting epoch of the observations
start=53187;  % Modified Julian Day, 'July 1, 2004'
%finish=start+ceil(365.25*5);  % approximately, set 5 yrs
deltaT=14;  % observation cadence, in days, set biweekly
%N=389;  % 14.9 yr, 128;  % number of biweekly observations for all pulsars, fft prefer 2^n
N=130;  % 5 yrs biweekly
%N=260;  % 10 yrs biweekly
dy=zeros(N,1);  % observation epoch, in day
yr=zeros(N,1);  % observation epoch, in year
for i=1:1:N
    dy(i)=start+(i-1)*deltaT;  % dates conducting observations, in MJD
    yr(i)=2004.5+(i-1)*deltaT*dy2yr;
end

simParams = struct('Np',Np,'N',N,'sd',sd,...
                  'alphaP',alphaP,'deltaP',deltaP,'kp',kp);


% -------------------------------
% calculate timing residuals induced by GW for each pulsar
%Amp=1.3*10^(-7);  % overall amplitude of timing residuals, the same for all pulsars, sec
timingResiduals=zeros(Np,N);  % signal, i.e. GW induced timing residuals, Np pulsars, N observations
timingResiduals_tmp=zeros(Np,N);   % signal without noise
phiI=zeros(Np,1);  % arbitrary phase for each pulsar, relative distance

Nrlz=50;  %500;  %200;  %25;  % number of noise realizations H1
noise=zeros(Np,N);  % noise
%sd=0.1*10^(-7);  % standard deviation of the normal distribtion (sec)

Nnis=2; %500; %100;  % number of realization of noise only cases H0

CA4filenames=cell(Nrlz*Nsnr+Nnis,1);  % cell array for file names of simulated data

perfect_fitness=0.0;  % fitness value for the true parameters

% set the range of the parameters
xmaxmin=zeros(7,2);  % x_max, x_min for each parameter x
xmaxmin(1,1)=2*pi;  % alpha
xmaxmin(1,2)=0.0;
xmaxmin(2,1)=pi/2;  % delta 
xmaxmin(2,2)=-pi/2;
xmaxmin(3,1)=20.0;  % angular velocity -- omega for GW
xmaxmin(3,2)=2.0;
xmaxmin(4,1)=pi;  % initial phase
xmaxmin(4,2)=0;
xmaxmin(5,1)=-8.0;  %-5.0;  %10^(-6);  % amplitude, in sec
xmaxmin(5,2)=-12.0;  % -10.0;  %10^(-8);
xmaxmin(6,1)=pi;  % inclination
xmaxmin(6,2)=0;
xmaxmin(7,1)=pi;  % polarization
xmaxmin(7,2)=0;

save(searchParamsFile,'xmaxmin');

%Standardized true parameter values
%stdTrueCoord = zeros(1,12);  % 4+Np 
stdTrueCoord = zeros(1,7); % 7 parameters other than pulsar phases

% calculate SNR
%snr_chr=0;  % characteristic srn = <\rho> = sqrt[(h,h)] srn of signal -- the strenghth of signal
%snr_mf=0;  % snr of matched filter = \rho = x/\sigma, normalized matched filter

nf=0;  % counter for number of files
for ii=6:1:Nsnr
    
    genHypothesis='H1 data';
    Amp=snr_tmp(ii)*2*10^(-9);
    
    for j=10:1:10  %Nloc
        
        % GW sky location in Cartesian coordinate
        k=zeros(1,3);  % unit vector pointing from SSB to source
        k(1)=cos(delta_tmp(j))*cos(alpha_tmp(j));
        k(2)=cos(delta_tmp(j))*sin(alpha_tmp(j));
        k(3)=sin(delta_tmp(j));
        
        for l=4:1:4  %Nomg
            
            snr_id=ii; %num2str(ii);
            loc_id=j; %num2str(j);
            omg_id=l; %num2str(l);
            
            snr_chr=0.0;  % initialize to zero
            
            % simulate the signals
            for i=1:1:Np
                
                theta=acos(k*kp(i,:)');
                %sprintf('%d pulsar theta=%g',i,theta)
                %phiI(i)=mod(phi0-omega*distP(i)*(1-cos(theta)), 2*pi);  % modulus after division
                %phiI(i)=mod(2*phi0-omega_tmp(l)*distP(i)*(1-cos(theta)), pi);  % modulus after division, YW 09/10/13
                phiI(i)=mod(phi0-0.5*omega_tmp(l)*distP(i)*(1-cos(theta)), pi);  % modulus after division, YW 04/30/14 check original def. of phiI
                
                %disp(['pulsar = ', num2str(i), ' ', num2str(phiI(i))])
                
                timingResiduals_tmp(i,:)=FullResiduals(alpha_tmp(j),delta_tmp(j),omega_tmp(l),phi0,phiI(i),alphaP(i),deltaP(i),...
                    Amp,iota,thetaN,theta,yr);
                
                %fftsignal(i,:)=fft(timingResiduals_tmp(i,:));

                % calculate the perfect fitness value
                
                snr_chr = snr_chr + dot(timingResiduals_tmp(i,:),timingResiduals_tmp(i,:)) / sd(i)^2;
                
                % standardization of the true coordinates
                stdTrueCoord(1)=(alpha_tmp(j)-xmaxmin(1,2))/(xmaxmin(1,1)-xmaxmin(1,2));  % [0, 2*pi]
                stdTrueCoord(2)=(delta_tmp(j)-xmaxmin(2,2))/(xmaxmin(2,1)-xmaxmin(2,2));  % [-pi/2, pi/2]
                stdTrueCoord(3)=(omega_tmp(l)-xmaxmin(3,2))/(xmaxmin(3,1)-xmaxmin(3,2));  % [2, 20]
                stdTrueCoord(4)= mod(phi0,pi)/pi;  % [0, pi]
                stdTrueCoord(5)=(log10(Amp)-xmaxmin(5,2))/(xmaxmin(5,1)-xmaxmin(5,2));
                stdTrueCoord(6)=(iota-xmaxmin(6,2))/(xmaxmin(6,1)-xmaxmin(6,2));
                stdTrueCoord(7)=(thetaN-xmaxmin(7,2))/(xmaxmin(7,1)-xmaxmin(7,2));
                
            end
            
            
            %snr_chr=sqrt(snr_chr/Np);  % averaged snr--root of mean square of individual snr
            snr_chr=sqrt(snr_chr);
            
            % signal + noise realizations
            for jj=1:1:Nrlz
                
                nf=nf+1;
                
                rlz_id=jj; %num2str(jj);
                
                % structure to store the id tag for each metadata file
                id=struct('snr_id',snr_id,'loc_id',loc_id,'omg_id',omg_id,'rlz_id',rlz_id);
                
                %snr_tmp=0;
                
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
                
                perfect_fitness = LLR_PSOmpp(stdTrueCoord,inParams);  % - LogLikelihoodRatio, minimization
                %true_fitness = fitnessTrue_ie(alpha_tmp(j),delta_tmp(j),omega_tmp(l),Amp,iota,thetaN,phi0,phiI,inParams);
                
                %disp(['In simulator2: perfect_fitness: ', num2str(perfect_fitness)]);
                
                % save metadata into a file for each realization (file name rule)
                filename=strcat('snr',num2str(ii),'loc',num2str(j),'omg',num2str(l),'rlz',num2str(jj),'.mat');
                snr=snr_tmp(ii);
                alpha=alpha_tmp(j);
                delta=delta_tmp(j);
                omega=omega_tmp(l);
                CA4filenames{nf}=filename;
                
                save([simDataDir,filesep,filename],'genHypothesis','snr','alpha','delta','omega','iota','thetaN','phi0',...
                    'timingResiduals','noise','yr', 'pname','id','simParams','perfect_fitness','snr_chr','timingResiduals_tmp','Amp');
                             
            end
            
        end
        
    end
    
end


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
        'timingResiduals','noise','yr','pname','id','simParams','perfect_fitness','snr_chr','timingResiduals_tmp');
    CA4filenames{nf+ii}=filename;
    
end


toc; % stop the stopwatch timer
% end of script
