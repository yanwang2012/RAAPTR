%% to generate random sources of GW
%% use Monte-Carlo simulation method
%% 2017.2.22 QYQ
function [Amp,alpha,delta,fgw,iota,Psi,Phi0,r]=GenerateRandomGWSource(Ns)
%Ns = 1000; % number of sourses
%% parameters to locate the source
log10Mc = random('uniform',6,10,1,Ns);
Mc = 10.^(log10Mc);
[alpha,delta]=SpherePointPicking(Ns);
log10fgw = random('uniform',-8,-6,1,Ns); % generate uniform distribution of log10(fgw)
%log10fgw = random('uniform',-9,-6,1,Ns); % generate uniform distribution of log10(fgw)
fgw = 10.^(log10fgw);% get the real fgw
% iota = random('uniform',0,pi,1,Ns);
iota = acos(random('uniform',-1,1,1,Ns));
Psi = random('uniform',0,pi,1,Ns);
Phi0 = random('uniform',0,pi,1,Ns);
%Mpc = 3.08567758*10^22;% define million parsec
R = 1000; % the max r of observe universe use the unit Mpc
pr0 = random('uniform',0,1,1,Ns);
pr1 = pr0.^(1/3);
%histogram(pr1);
%figure;
r = R*pr1;
%Amp = 1./r .* Mc .^ (5/3) .* (2 .* pi .* fgw) .^ (-1/3);%
%Amp = 1.43e-23 .* 1./r .* Mc .^ (5/3) .* (2 .* pi .* fgw) .^ (-1/3);% the amplitude of GW
Amp = 9.7876e-24 .* 1./r .* Mc .^ (5/3) .* (2 .* pi .* fgw) .^ (-1/3);% the amplitude of GW
% x = r.*cos(alpha).*sin(delta);
% y = r.*sin(alpha).*sin(delta);
% z = r.*cos(delta);
% plot3(x,y,z,'.r');