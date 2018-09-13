% log likelihood ratio for a single pulsar
% 05/26/14 YW
function [phiI,lhI]=likelihood(y,inn)

% there are two solution (+-x,y) for each y=cos(2\phi)
x(1)=sqrt(1-y^2);
x(2)=-sqrt(1-y^2);
lh=zeros(2,1);

for i=1:1:2
    lh(i) = inn(1)*y + inn(2)*x(i) + inn(3) - 0.5*(inn(4)*y^2 + inn(7)*x(i)^2 ...
        + 2*inn(5)*x(i)*y + 2*inn(6)*y + 2*inn(8)*x(i) + inn(9));  % Eq 22
end

if lh(1)>lh(2)
    tmp=x(1);
    lhI=lh(1);
else
    tmp=x(2);
    lhI=lh(2);
end
    
%if y>=0
if tmp>=0
    phiI=atan2(tmp,y)/2;  % y=cos(2*phiI)==X in Matlab atan2 func
else
    phiI=(atan2(tmp,y)+2*pi)/2;
end

% disp('In likelihood: ')
% disp([y,tmp]);
% disp([cos(2*phiI),sin(2*phiI)]);

% End of function