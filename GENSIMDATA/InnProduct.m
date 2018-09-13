% Function to calculate the noise weightred inner product (ip)
% white Gaussian noise with standard variance sd
% Yan Wang, March 4, 2013.

function ip=InnProduct(X,Y,sd)

ip=dot(X,Y)/sd^2;  % vector dot product

%ip=(X.*Y)/sd^2;


% END of function