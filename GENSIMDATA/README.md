Updated 09/13/2018 YW


This folder contains the code for generating simulated PTA data.


'simulator_ska_GWB.m' is the simulator that generates the data for the SKA era PTA. 'H1 data' (such as 'GWBsimDataSKASrlz1Nrlz1.mat') contains the timing residuals of SKA observed pulsars (Np=1000) including both the GW induced residuals from a number of (Ns=100) SMBHBs across the Universe and the measurement noise (sd=100 ns, contributed by radiometer and jitter noise); 'H0 data' (such as 'noise1.mat') contains noise only timing residuals.  

This simulation has a few outputs: 

(1) In simDateDir folder ./GWBsimDataSKA_xxx: GWB_Srlz1.mat contains the information of the simulated SMBHBs, this info may have overlaps with the ones in the H1/H0 data for the ease of subsequent analysis. Each data may take 20 sec or so on a single processor (2.8 GHz Intel Core i7).

(2) ./GWBsimDataSKA_xxx: H1 data GWBsimDataSKASrlz1Nrlz1.mat


(3) ./GWBsimDataSKA_xxx: H0 data noise1.mat (We try to maintain the same data structure of H0 data as H1 data.)


(4) In home folder: search parameter file (searchParamsFile='searchParams_GWBsimDataSKA_xxx') contains variable 'xmaxmin' which provides the maximum and minimum values of the searched parameters (parameter ranges). The detection algorithm 'MaxPhase' and 'AvPhase' (in C) may use a different parameter file which separate the full range into several bins for the purpose of multiple source detection and estimation. 


Inputs needed: 

(1) Pulsar catalog: skamsp=load('/Users/XXX/YYY/survey_ska.mat');

(2) Subfunctions called: 

    * GenerateRandomGWSource, SpherePointPicking: generates a random realization of Ns SMBHBs with parameters, such as chirp mass (Mc), location (alpha, delta), distance (r), GW frequency (fgw), angles (iota, Psi, Phi0), sampled either uniformly or log uniformly. These parameters are saved in H1 data file. The parameter values are set to 0 for H0 data. 

    * coco, cosine, rotm_coo: coordinates transformation

    * likelihood, InnProduct,LLR_PSO, etc:  calculating the time series of timing residual and the likelihood. 

Note: simulator_skaNano.m is an earlier version of the code, which is not used in the current simulation. 

