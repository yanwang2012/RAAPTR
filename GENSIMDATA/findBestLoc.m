% post analysis for multiple gw sources
% comparing the bestRealLoc with the injected gw sources
%% Yi-Qian 23th, Mar, 2019
clear;
ana_loc = load('summary.mat');
inj_gws = load('~/Research/PulsarTiming/SimDATA/INPUTDATA/SimSrlz1/GWBsimDataSKASrlz1Nrlz1.mat');
Nl = length(ana_loc.bestRealLoc);
Ns = length(inj_gws.alpha);
detec = zeros(3,1); % detected location alpha, delta and frequency omega
for i =1:Nl
    for j = 1:Ns
        if abs(ana_loc.bestRealLoc(3,i)-inj_gws.omega(j)) <= 0.1 % match the frequency
            detec(3,i) = ana_loc.bestRealLoc(3,i);
            display(['The number of source is ',num2str(j)])
            if abs(ana_loc.bestRealLoc(1,i)-inj_gws.alpha(j)) <= 0.1 && ana_loc.bestRealLoc(2,i)-inj_gws.delta(j) <= 0.1 % match alpha and delta
                detec(1,i) = ana_loc.bestRealLoc(i,i);
                detec(2,i) = ana_loc.bestRealLoc(2,i);
            end
        end
    end
end