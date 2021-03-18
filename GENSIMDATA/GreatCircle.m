function [v]=GreatCircle(center,pend1,pend2,r)
% Connect two points on a sphere with great circle.
% [V] = GREATCIRCLE(CENTER,ENDPOINT1,ENDPOINT2,RADIUS)
% the output v vector contains the coordinates of the great circle
% v(x,y,z). CENTER is the center of the sphere CENTER(x0,y0,z0), and two
% endpoints, ENDPOINT1(x1,y1,z1) and ENDPOINT2(x2,y2,z2). RADIUS is the
% radius of the sphere.

% QYQ 2021/3/18

v1 = [pend1(1)-center(1);pend1(2)-center(2);pend1(3)-center(3)]; % Vector from center to 1st point
% r = norm(v1); % The radius
v2 = [pend2(1)-center(1);pend2(2)-center(2);pend2(3)-center(3)]; % Vector from center to 2nd point
v3 = cross(cross(v1,v2),v1); % v3 lies in plane of v1 & v2 and is orthog. to v1
v3 = r*v3/norm(v3); % Make v3 of length r
% Let t range through the inner angle between v1 and v2
t = linspace(0,atan2(norm(cross(v1,v2)),dot(v1,v2)));
v = v1*cos(t)+v3*sin(t); % v traces great circle path, relative to center