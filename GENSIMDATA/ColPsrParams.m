function [pulsarParams]=ColPsrParams(path_to_pulsar_catalog)
% Collect pulsar paramters for function Amp2Snr
% [PulsarParams]=ColPsrParams(PATH_TO_PULSAR_CATALOG)

% QYQ 9th April, 2019
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

pulsarParams = struct('Np',Np,'N',N,'phiI',phiI,'alphaP',alphaP,'deltaP',deltaP,...
                      'kp',kp,'distP',distP,'yr',yr,'sd',sd);