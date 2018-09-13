% function calculating the c coefficients
% optimized the calculation for x, y, z.  11/20/2014

function [c,varargout]=cfunc(alpha,delta,alphaP,deltaP,theta,Amp,omega,iota,thetaN,phi0,Phi,s,sdi)
%function c = cfunc(alpha,delta,alphaP,deltaP,theta,Amp,omega,iota,thetaN,phi0,Phi,s,sdi)

alphatilde=alpha-alphaP;

Pp=-cos(deltaP)^2*(1-2*cos(alphatilde)^2+cos(alphatilde)^2*cos(delta)^2)...
   +sin(deltaP)^2*cos(delta)^2-0.5*sin(2*deltaP)*cos(alphatilde)*sin(2*delta);

Pc=2*cos(deltaP)*sin(alphatilde)*(cos(deltaP)*cos(alphatilde)*sin(delta)...
   -sin(deltaP)*cos(delta));

Fp=Pp/(1-cos(theta));
Fc=Pc/(1-cos(theta));

A=2*Amp*sqrt( (1+cos(iota)^2)^2*(Fp*cos(2*thetaN)-Fc*sin(2*thetaN))^2 ...
        + 4*cos(iota)^2*(Fp*sin(2*thetaN)+Fc*cos(2*thetaN))^2 );

%tmp=-2*cos(iota)/(1+cos(iota)^2)*(Fp*sin(2*thetaN)+Fc*cos(2*thetaN))/(Fp*cos(2*thetaN)-Fc*sin(2*thetaN));
% solve psi, atan or atan2 ?
%psi=atan(tmp);
psi=atan2( -2*cos(iota)*(Fp*sin(2*thetaN)+Fc*cos(2*thetaN)), ...
           (1+cos(iota)^2)*(Fp*cos(2*thetaN)-Fc*sin(2*thetaN)) );

% sin_ts=sin(phi0+psi+Phi);  % N by 1 matrix
% cos_ts=cos(phi0+psi+Phi);
% 
% B=0.5*A*sin(phi0)*sin_ts;
% C=-0.5*A*cos(phi0)*sin_ts;
% D=0.5*A*sin(phi0)*cos_ts;
% E=-0.5*A*cos(phi0)*cos_ts;

% x=B-E;
% y=C+D;  % 5/5/2014
% z=B+E;

x=0.5*A*cos(psi+Phi);
y=-0.5*A*sin(psi+Phi);
z=-0.5*A*cos(2*phi0+psi+Phi);  % 11/20/2014

% c are combination of inner weighted product of s with X,Y,Z
sx=InnProduct(s,x,sdi);  % scalar
sy=InnProduct(s,y,sdi);
sz=InnProduct(s,z,sdi);
xx=InnProduct(x,x,sdi);
xy=InnProduct(x,y,sdi);
xz=InnProduct(x,z,sdi);
yy=InnProduct(y,y,sdi);
yz=InnProduct(y,z,sdi);
zz=InnProduct(z,z,sdi);

c(1)=-sx+xz;
c(2)=sy-yz;
c(3)=0.5*(xx-yy);
c(4)=-xy;

if nargout > 1
    varargout{1}=[sx,sy,sz,xx,xy,xz,yy,yz,zz];
end

% end of function