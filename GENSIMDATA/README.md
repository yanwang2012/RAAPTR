This folder contains the code for generating simulated Pulsar Timing Array (PTA) data. The data from a PTA is a collection of the timing residuals from Np pulsars. The code generates statistically independent realizations of PTA data and stores them in files.

The main code is **gensimdata.m**, which generates data realizations of two types, namely,  signal+noise ('H1 data')   and noise-only ('H0 data'). The name of the file containing a data realization indicates what type of data it is. For example, 'GWBsimDataSKASrlz1Nrlz1.mat' contains H1 data while  'noise1.mat' contains H0 data. In each file, the timing residuals for all the pulsars are generated with the same start and end times, and the same sampling rate. Each data file may take 20 sec or so to generate on a single processor (2.8 GHz Intel Core i7).

The noise in each timing residual is assumed to consist of only radiometer and jitter noise, which is well approximated as Gaussian white noise, because these are expected to be the dominant sources of noise for timing residuals obtained with the Square Kilometer Array (SKA). At present, the noise standard deviation (sd) is assumed to be 100 nano-seconds (ns) for each pulsar. The GW signal in each timing residual arises from a set of Ns Supermassive Black Hole Binaries (SMBHBs) scattered across the Universe. (Currently, Ns=100 in the codes.)

The simulated data files are stored in an output_folder. Its contents are organized as follows.

(1) The file GWB_Srlz1.mat contains information about the particular simulation. This info may have some overlap with that cointained in the H1/H0 data files for ease of subsequent analysis.

(2) H1 data files (e.g., GWBsimDataSKASrlz1Nrlz1.mat).


(3) H0 data files (e.g., noise1.mat). We try to maintain the same variable names for the contents in H0 and H1 data files.

Example and more information

you can find an example and more information about the functions at
[user-guide](https://github.com/yanwang2012/RAAPTR/blob/master/GENSIMDATA/gensimdata_Guide.md)
