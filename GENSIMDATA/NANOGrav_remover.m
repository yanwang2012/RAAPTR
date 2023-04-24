%% PTA data preprocessing for SHAPES ("Block Averaging")
%Preprocessing of PTA data prior to running SHAPES for trend estimation.
%The steps in this script are as follows. Block average: average clusters
%of observations falling within a specified short period (e.g., 1 day).
%Interpolation: interpolate the block averaged time series to a uniform
%grid in time
%Scaling: Scale the data values to make the noise standard deviation close
%to unity.
%Trend estimation: run SHAPES to estimate the trend.
%De-scale: Remove the scaling used above from the estimated trend.

clear;
% Simulated data generated by YQ with 1 source injected. 
DataDir = '/Users/yiqianqian/Library/Mobile Documents/com~apple~CloudDocs/Research/PulsarTiming/SimDATA/Lomb-Scargle';
load([DataDir,filesep,'NANOGravSimDataNrlz1.mat']);

%Pulsar number
psrNum = 2; % B1937+21

%Pulsar Name
p_name = psr_data{psrNum}.psr_name;

%Duration of block in which data points will be averaged
mrgBlck = 1/365;%year

%Sampling period for interpolation after block averaging
intrpSmplIntrvl = 10/365;%year

%SHAPES parameters
rGain = 0.01; %regulator gain
nKnts = 5:5:30; %Knot numbers
nRuns = 4; %PSO runs

%----------DO NOT CHANGE BELOW-----------
timeVec = psr_data{psrNum}.yr;
dataVec = psr_data{psrNum}.timingResiduals;

mrgTimeVec = -ones(size(timeVec)); %unused elements will be discarded
mrgDataVec = -ones(size(dataVec));

nSmpls = length(dataVec);

counterNotAv = 1; %Counter for earliest observation not averaged
counterAv = 0; %Counter for latest observation from averaging
while counterNotAv <= nSmpls
    %Indices of samples belonging to the next block to be averaged
    mrgIndx = timeVec >= timeVec(counterNotAv) &   timeVec  <= timeVec(counterNotAv) + mrgBlck;
    if ~isempty(mrgIndx)
        counterAv = counterAv + 1;
        mrgDataVec(counterAv) = mean(dataVec(mrgIndx));
        mrgTimeVec(counterAv) = mean(timeVec(mrgIndx));
        %disp(diff(timeVec(mrgIndx)));
    end
    %Move past the block that was averaged: number of elements in this
    %block are the number of 1's in mrgIndx
    counterNotAv = counterNotAv + sum(mrgIndx);
end
mrgTimeVec((counterAv+1):end) = [];
mrgDataVec((counterAv+1):end) = [];

%Interpolate
intrpTimeVec = mrgTimeVec(1):intrpSmplIntrvl:mrgTimeVec(end);
intrpData = interp1(mrgTimeVec, mrgDataVec, intrpTimeVec);

%Scale
intrpDataScl = intrpData*1e7;

%Trend estimation
[outStrct,bestModelOutStrct] = shps(struct('dataX',intrpTimeVec,'dataY',intrpDataScl,'nBrks',nKnts,'rGain',0.01),...
               struct('nRuns',nRuns,'psoParams',[]));
           
%De-scale
trndEst = bestModelOutStrct.bestModelSig/1e7;


%Evaluate at original time values
intrpTrndEst = interp1(intrpTimeVec,trndEst,timeVec);

%Residual after subtracting out the trend estimate
resOfres = dataVec - intrpTrndEst;

save([DataDir,filesep,num2str(psrNum),'_detrend.mat'],'resOfres');

%Plots
close all;
figure;
hold on
plot(timeVec,dataVec,'o');
plot(mrgTimeVec, mrgDataVec,'+');
legend('PTA Data','Block-averaged');
axis tight;

figure;
hold on;
plot(timeVec,dataVec,'o');
plot(mrgTimeVec, mrgDataVec,'+');
plot(intrpTimeVec,intrpData,'LineWidth',2.0);
legend('PTA Data','Block-averaged','Block-averaged Interpolated');
axis tight;

figure;
hold on;
plot(timeVec,dataVec,'o');
plot(intrpTimeVec, trndEst,'LineWidth',2.0);
plot(timeVec, resOfres,'.','MarkerSize',4.0);
legend('PTA Data','SHAPES estimate','Residual');
axis tight;





