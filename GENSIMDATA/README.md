Updated 09/13/2018 YW


This folder contains the code for generating simulated Pulsar Timing Array (PTA) data. A PTA data set consists of the timing residuals from Np pulsars. The timing residuals for all the pulsars are generated with the same start and end times and the same sampling rate. Henceforth, we refer to the collection of timing residuals as a data realization, or just the data.


The main code is 'simulator_ska_GWB.m', which generates multiple data realizations of two types, namely,  signal+noise ('H1 data')   and noise-only ('H0 data'). The name of the file containing a data realization indicates what type of data it is. For example, 'GWBsimDataSKASrlz1Nrlz1.mat' contains H1 data while  'noise1.mat' contains H0 data.

The noise is assumed to consist of only radiometer and jitter noise, which is well approximated as Gaussian white noise, because these are expected to be the dominant sources of noise for timing residuals obtained with the Square Kilometer Array (SKA). At present, the noise standard deviation (sd) is assumed to be 100 nano-seconds (ns) for each pulsar.

The GW signal in each timing residual arises from a set of Ns=100 Supermassive Black Hole Binaries (SMBHBs) scattered across the Universe. 

The simulated data files are organized as follows: 

(1) In the simDataDir folder named as ./GWBsimDataSKA_xxx: GWB_Srlz1.mat contains information about the particular simulation. This info may have some overlap with that cointained in the H1/H0 data files for ease of subsequent analysis. Each data file may take 20 sec or so to generate on a single processor (2.8 GHz Intel Core i7).

(2) ./GWBsimDataSKA_xxx: H1 data (e.g., GWBsimDataSKASrlz1Nrlz1.mat)


(3) ./GWBsimDataSKA_xxx: H0 data (e.g., noise1.mat). We try to maintain the same variable names for the contents in H0 and H1 data files.


(4) In the home folder: search parameter file (searchParamsFile='searchParams_GWBsimDataSKA_xxx') contains variable 'xmaxmin' which provides the maximum and minimum values of the searched parameters (parameter ranges). The detection algorithm 'MaxPhase' and 'AvPhase' (in C) may use a different parameter file which separate the full range into several bins for the purpose of multiple source detection and estimation. 


Inputs needed: 

(1) Pulsar catalog: skamsp=load('/Users/XXX/YYY/survey_ska.mat');

(2) Subfunctions called: 

    * GenerateRandomGWSource, SpherePointPicking: generates a random realization of Ns SMBHBs with parameters, such as chirp mass (Mc), location (alpha, delta), distance (r), GW frequency (fgw), angles (iota, Psi, Phi0), sampled either uniformly or log uniformly. These parameters are saved in H1 data file. The parameter values are set to 0 for H0 data. 

    * coco, cosine, rotm_coo: coordinates transformation

    * likelihood, InnProduct,LLR_PSO, etc:  calculating the time series of timing residual and the likelihood. 

Note: simulator_skaNano.m is an earlier version of the code, which is not used in the current simulation. 

