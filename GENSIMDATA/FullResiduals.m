% Function to calculate the timing residuals for a full set of parameters 
% (initial phase, arbitrary phase ...)
% Yan Wang, March 4, 2013.
% /07/25/13, Vectorization of code

function r=FullResiduals(alpha,delta,omega,phi0,phiI,alphaP,deltaP,Amp,iota,thetaN,theta,t)

N=length(t);

%r=zeros(N,1);
c=zeros(8,1);  % coefficient function
A=zeros(N,8);  % base functions

%sprintf('len N=%d',N)
%theta=acos(k.*kp);
%theta=acos(dot(k,kp));
%Delta=distP*(1-cos(theta));

alphatilde=alpha-alphaP;

Pp=-cos(deltaP)^2*(1-2*cos(alphatilde)^2+cos(alphatilde)^2*cos(delta)^2)...
   +sin(deltaP)^2*cos(delta)^2-0.5*sin(2*deltaP)*cos(alphatilde)*sin(2*delta);

Pc=2*cos(deltaP)*sin(alphatilde)*(cos(deltaP)*cos(alphatilde)*sin(delta)...
   -sin(deltaP)*cos(delta));

Fp=Pp/(1-cos(theta));
Fc=Pc/(1-cos(theta));

%disp(['in FullRes..','Pp',num2str(Pp),'Pc',num2str(Pc),'theta',num2str(theta)])

CosIota=cos(iota);
TwoThetaN=2*thetaN;
tmp1=Amp*(1+CosIota^2);
tmp2=Amp*2*CosIota;
tmpC=cos(TwoThetaN);
tmpS=sin(TwoThetaN);
% here Amp=\zeta*omega^(-1/3) that defined in my notes
c(1,1)=-tmp1*tmpC;
c(2,1)=-tmp2*tmpS;
c(3,1)=-c(1,1);
c(4,1)=c(2,1);
c(5,1)=tmp1*tmpS;
c(6,1)=-tmp2*tmpC;
c(7,1)=-c(5,1);
c(8,1)=c(6,1);
%
% c(1,1)=-Amp*(1+cos(iota)^2)*cos(2*thetaN);
% c(2,1)=-Amp*2*cos(iota)*sin(2*thetaN);
% c(3,1)= Amp*(1+cos(iota)^2)*cos(2*thetaN);
% c(4,1)=-Amp*2*cos(iota)*sin(2*thetaN);
% c(5,1)= Amp*(1+cos(iota)^2)*sin(2*thetaN);
% c(6,1)=-Amp*2*cos(iota)*cos(2*thetaN);
% c(7,1)=-Amp*(1+cos(iota)^2)*sin(2*thetaN);    
% c(8,1)=-Amp*2*cos(iota)*cos(2*thetaN);

%phiI=phi0-omega*distP*(1-cos(theta));
tmpC2=cos(2*phi0)-cos(2*phiI);
tmpS2=sin(2*phi0)-sin(2*phiI);
FpC=Fp*tmpC2; % (cos(2*phi0)-cos(2*phiI));
FpS=Fp*tmpS2; % (sin(2*phi0)-sin(2*phiI));
FcC=Fc*tmpC2; % (cos(2*phi0)-cos(2*phiI));
FcS=Fc*tmpS2; % (sin(2*phi0)-sin(2*phiI));
% FpC=Fp*(cos(2*phi0)-cos(2*phiI));
% FpS=Fp*(sin(2*phi0)-sin(2*phiI));
% FcC=Fc*(cos(2*phi0)-cos(2*phiI));
% FcS=Fc*(sin(2*phi0)-sin(2*phiI));
%OmegaT=zeros(N,1);
OmegaT=omega*t;

A(:,1)=FpC*cos(OmegaT);
A(:,2)=FpS*cos(OmegaT);
A(:,3)=FpS*sin(OmegaT);
A(:,4)=FpC*sin(OmegaT);
A(:,5)=FcC*cos(OmegaT);
A(:,6)=FcS*cos(OmegaT);
A(:,7)=FcS*sin(OmegaT);
A(:,8)=FcC*sin(OmegaT);

r=A*c;

% for i=1:1:N
%     
%     A(i,1)=Fp*(cos(2*phi0)-cos(2*phiI))*cos(omega*t(i));
%     A(i,2)=Fp*(sin(2*phi0)-sin(2*phiI))*cos(omega*t(i));
%     A(i,3)=Fp*(sin(2*phi0)-sin(2*phiI))*sin(omega*t(i));
%     A(i,4)=Fp*(cos(2*phi0)-cos(2*phiI))*sin(omega*t(i));
%     A(i,5)=Fc*(cos(2*phi0)-cos(2*phiI))*cos(omega*t(i));
%     A(i,6)=Fc*(sin(2*phi0)-sin(2*phiI))*cos(omega*t(i));
%     A(i,7)=Fc*(sin(2*phi0)-sin(2*phiI))*sin(omega*t(i));
%     A(i,8)=Fc*(cos(2*phi0)-cos(2*phiI))*sin(omega*t(i));
%     
%     for j=1:1:8 
%         r(i)=r(i)+c(j)*A(i,j);
%     end
%     
% end

% END of function