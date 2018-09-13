% Function to calculate the likelihood ratio for particle swarm 
% optimazition (PSO) code developed by Prof. Soumya Mohanty.
% Maximize pulsar phases. Yan Wang, Jan 21, 2014.

function [fitVal,varargout]=LLR_PSOmpp(xVec,inParams)
%function fitVal=LLR_PSO(xVec,inParams)

%rows: points
%columns: coordinates of a point
[nrows,ncols]=size(xVec);

%storage for fitness values
fitVal = zeros(nrows,1);

validPts = chkstdsrchrng(xVec);
% %Set fitness for invalid points to infty
fitVal(~validPts)=inf;
% if isempty(inParams)
%     %use default range of coordinates
%     %(only the numerical values below should be changed for different
%     %fitness functions)
%     xVec(validPts,:) = s2rscalar(xVec(validPts,:),-5,5);
% else
%     xVec(validPts,:) = s2rvector(xVec(validPts,:),inParams);
% end

%x=zeros(1,ncols);
% ===============================
%realCoord = zeros(size(xVec));
Np=inParams.Np;
realCoord = zeros(1,7+Np);
for lpc = 1:nrows
    if validPts(lpc)
    % Only the body of this block should be replaced for different fitness
    % functions
        x = xVec(lpc,:);
        % x(1) alpha; x(2) delta; x(3) omega; x(4) phi0
        % x(5) phi1; x(6) phi2; x(7) phi3; x(8) phi4 ...
        % 0<=x<=1 is convert to physical unit in 'LogLikelihoodRatio'
            %[ft,dummy]=LogLikelihoodRatioMP(x,inParams);  % maximum phase
            %[ft,dummy]=LogLikelihoodRatioMP4(x,inParams);  % maximum phase, both fzero and roots used
            [ft,dummy]=LogLikelihoodRatioMP5(x,inParams);  % maximum phase
            %[ft,dummy]=LogLikelihoodRatioNonc(x,inParams);
            %[ft,dummy1,dummy]=LogLikelihoodRatioNonc2(x,inParams);
            fitVal(lpc)= ft;
            realCoord(lpc,:)=dummy;
            %a=dummy1;
    end
end
% ===============================


%Return real coordinates if requested
if nargout > 1
    varargout{1}=realCoord;
end
% END of function