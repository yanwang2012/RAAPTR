% antenna pattern of a PTA, response matrix and its condition number 
% ref 'NetRespMat.m' for a previous version
% YW, May 24, 2016

clear;

% set path to the data
%simDataDir='/Users/ywang/Research/PULSARTIMING/MultiCW/simData17_snr123_loc3/';
simDataDir='/Users/qyq/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/BANDEDGE/Band_opt/simData/Acond for SKA';
%resDataDir='/Users/ywang/Research/PULSARTIMING/MultiCW/simData17_snr123_loc3/results_nelderMead/';

% load data for the PTA configeration
filename='snr3loc3omg3rlz11.mat';
simData=load([simDataDir,filesep,filename]);
%resData=load([resDataDir,filesep,filename]);

% original
Na=300;
Nd=151;  % count both north and sourth poles

% in order to generate Acond for 3D skymap
% Na = 301; % for Na-1 faces sphere
% Nd = 301;

alpha=zeros(Na,1);  % ra of source
delta=zeros(Nd,1);  % dec of source
ks=zeros(1,3);

Acond=zeros(Nd,Na);  % condition number of A at a sky location
Inten=zeros(Nd,Na);  % sqrt(Fp^2+Fc^2)

Np=simData.simParams.Np;
alphaP=simData.simParams.alphaP;
deltaP=simData.simParams.deltaP;
kp=simData.simParams.kp;

A=zeros(Np,2);  % network response matrix: Np pulsars, 2 GW polarizations +,x
Fp=zeros(Np,1);  % antenna pattern func for +,  range vector of A
Fc=zeros(Np,1);  % antenna pattern func for x


for k=1:1:Nd
    delta(k,1)=pi/2-(k-1)*pi/(Nd-1);
    
    for j=1:1:Na
        alpha(j,1)=(j-1)*2*pi/(Na-1);
        
        ks(1)=cos(delta(k))*cos(alpha(j));
        ks(2)=cos(delta(k))*sin(alpha(j));
        ks(3)=sin(delta(k));
        
        for i=1:1:Np
            alphatilde=alpha(j)-alphaP(i);
            theta=acos(ks*kp(i,:)');
                 
            Pp=-cos(deltaP(i))^2*(1-2*cos(alphatilde)^2+cos(alphatilde)^2*cos(delta(k))^2)...
                +sin(deltaP(i))^2*cos(delta(k))^2-0.5*sin(2*deltaP(i))*cos(alphatilde)*sin(2*delta(k));
            
            Pc=2*cos(deltaP(i))*sin(alphatilde)*(cos(deltaP(i))*cos(alphatilde)*sin(delta(k))...
                -sin(deltaP(i))*cos(delta(k)));
            
            % + polarization
            Fp(i)=Pp/(1-cos(theta));
            A(i,1)=Fp(i);
            
            % x polarization
            Fc(i)=Pc/(1-cos(theta));
            A(i,2)=Fc(i);
            
            %Inten(j,k)=Inten(j,k) + sqrt(Fp(i)^2 + Fc(i)^2);
                 
        end     
        
        Inten(k,j)=sqrt(Fp(15)^2 + Fc(15)^2);
        % calculate the condition number
        Acond(k,j)=cond(A,2);  % cond(A,inf)
        
    end
    
end


[v,ind]=max(Acond);
[v1,ind1]=max(max(Acond));
disp(['The largest element in Acond is: ' num2str(v1) ' at (' num2str(ind(ind1)) ',' num2str(ind1) ')']);
disp(['The values for ra and dec are: ', num2str(delta(ind(ind1))),' and ',num2str(alpha(ind1))]);

% skymap of the condition number
figure
%surf(alpha,delta,Acond);
imagesc(alpha,delta,Acond);
xlim([0 2*pi])
ylim([-pi/2 pi/2])
xlabel('\alpha')
ylabel('\delta')
%title(['Skymap for the condition number of A for ', num2str(Np), ' pulsars']);
%surf(alpha,delta,log10(Acond)');  % log10 plot

figure
surf(alpha,delta,Inten);
xlim([0 2*pi])
ylim([-pi/2 pi/2])
xlabel('\alpha')
ylabel('\delta')
title(['Skymap for Inten for ', num2str(Np), ' pulsars']);
% imagesc

%% 3D skymap
[x,y,z] = sph2cart(alpha',delta,ones(Nd,Na));
% [x,y,z] = sphere(Na-1);
figure
surf(x,y,z,Acond,'EdgeColor','Interp')
axis equal
c = colorbar;
c.Label.String = 'Condition Number';
saveas(gcf,[simDataDir,filesep,'3DSkymap.png'])
savefig([simDataDir,filesep,'3DSkymap'])

% End