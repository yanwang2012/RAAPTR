function []=parameters(NumGWsources,NumPulsar,NumNoiseReali,NumRealiNoise)
%% A script to generate the all the parameters and save it to .mat
% for the function FSGWB
% Yi-Qian, Sep 16, 2018
%% parameters
% ==== Generate random GW sources ====
Ns = NumGWsources;  % number of GW sources
[Amp,alpha_tmp,delta_tmp,fgw,iota,thetaN,phi0,r]=GenerateRandomGWSource(Ns);
omega_tmp = 2*pi* fgw * 3.156*10^7;  % convert sec^-1 (Hz) to yr^-1
Np=NumPulsar;  % number of pulsars in the timing array

% starting epoch of the observations
start=53187;  % Modified Julian Day, 'July 1, 2004'
%finish=start+ceil(365.25*5);  % approximately, set 5 yrs
deltaT=14;  % observation cadence, in days, set biweekly
%N=389;  % 14.9 yr, 128;  % number of biweekly observations for all pulsars, fft prefer 2^n
N=130;  % 5 yrs biweekly
%N=260;  % 10 yrs biweekly
dy=zeros(N,1);  % observation epoch, in day
yr=zeros(N,1);  % observation epoch, in year

% noise
Nrlz=NumNoiseReali;  % number of noise realizations H1
noise=zeros(Np,N);  % noise
%sd=0.1*10^(-7);  % standard deviation of the normal distribtion (sec)

Nnis=NumRealiNoise;  % number of realization of noise only cases H0

% set the range of the parameters
xmaxmin=zeros(7,2);  % x_max, x_min for each parameter x
xmaxmin(1,1)=2*pi;  % alpha
xmaxmin(1,2)=0.0;
xmaxmin(2,1)=pi/2;  % delta 
xmaxmin(2,2)=-pi/2;
xmaxmin(3,1)=200.0;  % angular velocity -- omega for GW
xmaxmin(3,2)=1.0;
xmaxmin(4,1)=pi;  % initial phase
xmaxmin(4,2)=0;
xmaxmin(5,1)=-6.0;  %-5.0;  %10^(-6);  % amplitude, in sec
xmaxmin(5,2)=-15.0;  % -10.0;  %10^(-8);
xmaxmin(6,1)=pi;  % inclination
xmaxmin(6,2)=0;
xmaxmin(7,1)=pi;  % polarization
xmaxmin(7,2)=0;

save('parameter.mat')
% end of function