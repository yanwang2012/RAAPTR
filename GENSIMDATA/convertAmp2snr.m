% Convert Amp back to SNR
% Yiqian Qian 2nd April, 2019
%% Load the information of pulsar
function [snr_chr]=convertAmp2snr(path_to_estimatedData,path_to_pulsar_catalog)
% path_to_simulationData = '~/Research/PulsarTiming/SimDATA/Mauritius/GWBsimDataSKA/GWBsimDataSKASrlz1Nrlz9.mat';
% path_to_estimatedData = '~/Research/PulsarTiming/SimDATA/Mauritius/Mauritius_results_mat/6_GWBsimDataSKASrlz1Nrlz9.mat';
% path_to_pulsar_catalog = '~/Research/PulsarTiming/GENSIMDATA/survey_ska.mat';
%% Load the information of estimated GW source
source = load(path_to_estimatedData,'bestRealLoc');

alpha = source.bestRealLoc(1);
delta = source.bestRealLoc(2);
omega = source.bestRealLoc(3);
phi0 = source.bestRealLoc(4);
Amp = source.bestRealLoc(5);
iota = source.bestRealLoc(6);
thetaN = source.bestRealLoc(7);


%% ==== Constructing a pulsar timing array using Np pulsars ====
% read in the pulsar catalogue simulated for SKA
%Get constants
dy2yr = genptaconsts('dy2yr');
Np = 1000; % number of pulsar
N = 130; % 5 yrs biweekly
% starting epoch of the observations
start=53187;  % Modified Julian Day, 'July 1, 2004'
%finish=start+ceil(365.25*5);  % approximately, set 5 yrs
deltaT=14;  % observation cadence, in days, set biweekly

skamsp=load(path_to_pulsar_catalog);% load input data
[~,I]=sort(skamsp.D);
distP=zeros(Np,1);  % (parallax) distance from SSB to pulsars, from mas/pc to ly
% sky location of the pulsars in the equatorial coordinate
% we need to transfer from hr angle and degree to radian
alphaP=zeros(Np,1);  % right ascension, in radian
deltaP=zeros(Np,1);  % declination, in radian
kp=zeros(Np,3);  % unit vector pointing from SSB to pulsars,
sd=zeros(Np,1);  % standard deviation of noise for different pulsar
%dy=zeros(N,1);  % observation epoch, in day
yr=zeros(N,1);  % observation epoch, in year
phiI=zeros(Np,1);  % arbitrary phase for each pulsar, relative distance

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

%% sky location of pulsars in Cartesian coordinate
for i=1:1:Np
    kp(i,1)=cos(deltaP(i))*cos(alphaP(i));
    kp(i,2)=cos(deltaP(i))*sin(alphaP(i));
    kp(i,3)=sin(deltaP(i));
end

for i=1:1:N
    %dy(i)=start+(i-1)*deltaT;  % dates conducting observations, in MJD
    yr(i)=2004.5+(i-1)*deltaT*dy2yr;
end

%% calculate SNR
snr_chr2_tmp=zeros(Np,1);  % squared characteristic srn for each pulsar and source
tmp=zeros(1,N); % noiseless timing residuals from a source
for j=1:1:Np  % number of pulsar
    %for j=1:1:Ns  % number of GW source
        
        % GW sky location in Cartesian coordinate
        k=zeros(1,3);  % unit vector pointing from SSB to source
        k(1)=cos(delta)*cos(alpha);
        k(2)=cos(delta)*sin(alpha);
        k(3)=sin(delta);
        theta=acos(k*kp(j,:)');
        phiI(j)=mod(phi0-0.5*omega*distP(j)*(1-cos(theta)), pi);  % modulus after division, YW 04/30/14 check original def. of phiI
        
        
        tmp = FullResiduals(alpha,delta,omega,phi0,phiI(j),alphaP(j),deltaP(j),...
            Amp,iota,thetaN,theta,yr);

        
        snr_chr2_tmp(j,1) = dot(tmp,tmp)/sd(j)^2;

end

snr_chr=sqrt(sum(snr_chr2_tmp,1));  % sum of elements in each column