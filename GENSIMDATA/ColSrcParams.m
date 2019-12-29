function [sourceParams]=ColSrcParams(path_to_estimatedData)
% Collect the source parameters for function Amp2Snr
% [sourceParams]=ColSrcParams(path_to_estimatedData)

% Yiqian Qian 2nd April, 2019
%% Load the information of pulsar
% path_to_simulationData = '~/Research/PulsarTiming/SimDATA/Mauritius/GWBsimDataSKA/GWBsimDataSKASrlz1Nrlz9.mat';
% path_to_estimatedData = '~/Research/PulsarTiming/SimDATA/Mauritius/Mauritius_results_mat/6_GWBsimDataSKASrlz1Nrlz9.mat';
% path_to_pulsar_catalog = '~/Research/PulsarTiming/GENSIMDATA/survey_ska.mat';
%% Load the information of estimated GW source
source = load(path_to_estimatedData,'bestRealLoc');
alpha = source.bestRealLoc(1);
delta = source.bestRealLoc(2);
omega = source.bestRealLoc(3);
phi0 = source.bestRealLoc(4);
Amp = source.bestRealLoc(5);
iota = source.bestRealLoc(6);
thetaN = source.bestRealLoc(7);
phiI = source.bestRealLoc(8:1007);

sourceParams = struct('alpha',alpha,'delta',delta,'omega',omega,'phi0',phi0,'Amp',Amp,...
                        'iota',iota,'thetaN',thetaN,'phiI',phiI);
% END of function