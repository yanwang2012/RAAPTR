function [theta,phi]=SpherePointPicking(n)
%Uniform random Sphere Point Picking

% r = 1;
NN = n;
u = random('uniform',0,1,1,NN);
v = random('uniform',0,1,1,NN);
theta = 2*pi*u;
phi = asin(2*v-1);

% x = zeros(3,1);
% y = zeros(3,1);
% z = zeros(3,1);
% for i=1:1:NN
% x(i) = r*sin(phi(i))*cos(theta(i));
% y(i) = r*sin(phi(i))*sin(theta(i));
% z(i) = r*cos(phi(i));
% end
% plot3(x,y,z,'.r');