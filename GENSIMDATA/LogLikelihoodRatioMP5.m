% log likelihood ratio function which maximizing over pulsar phase-MP,
% this function should be implemented as efficient as possible, since it
% would be called for a large amount of times, single source
% 01/21/14, Yan Wang: adopt and modify from LogLikelihoodRatio()
% 11/11/14, add fzero() after roots(), defaul TolX for fzero() is eps=2.2204e-16
% 11/13/14, compare fitness values at boundaries and at the stationary points
% 12/19/14, cancel the changes made on 11/11/14, but still check boundary

function [LLR,varargout]=LogLikelihoodRatioMP5(x,inParams)

%% transfer parameters from structure inParams
Np = inParams.Np;
%N = inParams.N;
s = inParams.s;
sd = inParams.sd;
alphaP = inParams.alphaP;
deltaP = inParams.deltaP;
kp = inParams.kp;
yr = inParams.yr;
xmaxmin = inParams.xmaxmin;

% transform x from standard coordinates [0,1] (requested by PSO) to physical 
% coordinates for 7 intrinsic parameters
alpha=x(1)*(xmaxmin(1,1)-xmaxmin(1,2))+xmaxmin(1,2);  % [0, 2*pi]
delta=x(2)*(xmaxmin(2,1)-xmaxmin(2,2))+xmaxmin(2,2);  % [-pi/2, pi/2]
omega=x(3)*(xmaxmin(3,1)-xmaxmin(3,2))+xmaxmin(3,2);
phi0=x(4)*(xmaxmin(4,1)-xmaxmin(4,2))+xmaxmin(4,2);
Amp=10^(x(5)*(xmaxmin(5,1)-xmaxmin(5,2))+xmaxmin(5,2));
iota=x(6)*(xmaxmin(6,1)-xmaxmin(6,2))+xmaxmin(6,2);
thetaN=x(7)*(xmaxmin(7,1)-xmaxmin(7,2))+xmaxmin(7,2);

%% calculate c for each pulsar
Phi=omega*yr;  % Phi=omega*t, N by 1 matrix
c=zeros(Np,4);  % Eq. 25
% sin_ts=zeros(N,1);
% cos_ts=zeros(N,1);
% B=zeros(N,1);
% C=zeros(N,1);
% D=zeros(N,1);
% E=zeros(N,1);

% sky location of source in Cartesian coordinate
k=zeros(1,3);  % unit vector pointing from SSB to source
%for i=1:1:Ns
k(1)=cos(delta)*cos(alpha);
k(2)=cos(delta)*sin(alpha);
k(3)=sin(delta);
%end

% quartic equation coefficients, closed form solution
e=zeros(Np,5);
lh=zeros(Np,6);  % likelihood for each pulsar  zeros(Np,4); two additions for boundaries
phiItmp=zeros(Np,6);  % tmp solution to phiI
phiI=zeros(Np,1);
LLR=0.0;  % log likelihood ratio

%rtmp=0.0;

for i=1:1:Np
    
    theta=acos( dot(k,kp(i,:)) );
    [c(i,:),inn]= cfunc( alpha,delta,alphaP(i),deltaP(i),theta,...
                           Amp,omega,iota,thetaN,phi0,Phi,s(i,:),sd(i) );
    
    e(i,1) = 4*(c(i,3)^2+c(i,4)^2);  % y^4
    e(i,2) = 4*(c(i,1)*c(i,3)+c(i,2)*c(i,4));
    e(i,3) = c(i,1)^2+c(i,2)^2-4*(c(i,3)^2+c(i,4)^2);
    e(i,4) = -2*c(i,2)*c(i,4)-4*c(i,1)*c(i,3);
    e(i,5) = c(i,4)^2-c(i,1)^2;  % y^0
    
    % calculate phiI(i) which maximize the likelihood ratio function
    % calculate roots of a quartic equation (there may be complex solutions)
    %r=zeros(4,1);
    r=roots([e(i,1) e(i,2) e(i,3) e(i,4) e(i,5)]);
    
%    fun=@(y) e(i,1)*y^4 + e(i,2)*y^3 + e(i,3)*y^2 + e(i,4)*y + e(i,5);  % quartic eq.
    
    nr=0;  % number of EFFECTIVE roots (real number && abs(r)<1), may be not the true roots for part(Lambda)
    for j=1:1:4
        if imag(r(j))==0 && abs(r(j))<=1.0  % real solution for y=cos(2\phi)
            
%            rtmp=fzero(fun,r(j));  % look for a refined solution close to r(j)
            
%             if rtmp>1.0
%                disp('Warning ... rtmp>1, need a break');
%                break          
%             end

%            if abs(rtmp)<=1.0  % fzero may make r(j) out of range
%                r(j)=rtmp;   % r(j) is refined
                [phiItmp(i,j),lh(i,j)] = likelihood(r(j),inn);  % phiI = [0,pi]
                nr=nr+1;
%            else  % no refined solution obtained, still use the ones from roots
                %disp(['rtmp>1, rtmp = ', num2str(rtmp), ', for PSR ', num2str(i)]);
%                 save('test.mat');
%                 pause;
                
%                [phiItmp(i,j),lh(i,j)] = likelihood(r(j),inn);  % phiI = [0,pi]
%                nr=nr+1;
%                disp(['In LogLikelihoodRatioMP4: fzero makes rtmp = ',num2str(rtmp),' for PSR: ',num2str(i)]);
%            end
            
        else
            phiItmp(i,j)=NaN('double');  % for complex root or abs > 1.0
            lh(i,j) = -Inf;  % ??
            %disp(['In LogLikelihoodRatioMP: There is a problem for psr: ',num2str(i)]);
        end
    end
    
    % find fitness values at boundaries even with nr!=0
    [~,tmp1]=likelihood(-1.0,inn); 
    [~,tmp2]=likelihood(1.0,inn);
    lh(i,5)=tmp1;
    lh(i,6)=tmp2;
    phiItmp(i,5)=pi/2;  % the same as the nr==0 part
    phiItmp(i,6)=0.0;
    
    if nr>0
        [C I]=max(lh(i,:));  % C is the largest value, I is index
        
%         if abs(C)==Inf
%             disp(['In LogLikelihoodRatioMP: C=Inf for pulsar: ',num2str(i)]);
%         end
        
        if I==5 || I==6
            disp(['Using fitness at boundary for PSR: ',num2str(i)]);
        end

        phiI(i)=phiItmp(i,I);
        LLR = LLR + C;
        
%         if abs(LLR)==Inf
%             disp(['In LogLikelihoodRatioMP: LLR=Inf for pulsar: ',num2str(i)]);
%         end
        
    elseif nr==0
        % there is no turning point for the log-likelihood function of this
        % pulsar, however, the extremum shoule exist at the boundary of the
        % allowed range of the variable
        
%         phiI(i)=NaN('double');
%         LLR = LLR + 0.0;  % this pulsar has no contribution to overal LLR

        disp(['In LogLikelihoodRatioMP5: NO effective root (nr=0) for PSR: ',num2str(i),', use boundary fitness']);
%         disp([alpha,delta,omega,phi0,Amp,iota,thetaN]);
        
%         [~,tmp1]=likelihood(-1.0,inn); 
%         [~,tmp2]=likelihood(1.0,inn);
        
%         if abs(tmp1)==Inf || abs(tmp2)==Inf
%             disp(['In LogLikelihoodRatioMP: ERROR for pulsar: ',num2str(i)]);
%         end
        
        if tmp1>tmp2
            C=tmp1;
            phiI(i)=pi/2;
        else
            C=tmp2;
            phiI(i)=0.0;
        end
        LLR = LLR + C;
        
    end
    
%     disp('In LogLikelihoodRatioMP: ')
%     %disp([alpha,delta,omega,phi0,Amp,iota,thetaN,phiI']');
%     disp([i,cos(2*phiI(i)),sin(2*phiI(i))]);
    
end

LLR = -LLR;  % return -log likelihood ratio

if nargout > 1
    varargout{1}=[alpha,delta,omega,phi0,Amp,iota,thetaN,phiI'];  % in real coord.
end

%disp()


% END of function