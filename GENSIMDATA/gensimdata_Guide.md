# Function 'gensimdata.m'
This is the main function to use for generating simulation data files.
## How to use
(See the 'test_gensimdata.m' script for a worked out example of how to generate simulated data files.)
The input argument for gensimdata are: The name of a .mat input parameter file containing the variables listed below, the name of a .mat file containing information about the pulsars in the PTA, ...
[SDM: test_gensimdata.m is an example. In this place, you are describing gensimdata.m. So, you should explain the arguments of gensimdata.m directly here and describe the contents of the files that it reads without asking the reader to edit test_gensimdata.m.]
### Variables required in input parameter file
**NumGwsources**

```matlab
number of GW sources.
```

**NumPulsar**

```matlab
number of pulsars.
```

**NumNoiseReali**

```matlab
number of noise realizations H1
```

**NumRealiNoise**

```matlab
number of realization of noise only cases H0
```
### Variables required in file containing PTA information
...

You choose these Parameters in the function **test_gensimdata.m** then it will generate a parameter file named **parameter.mat**. In this file, you have all the parameters you need to simulate the data set. To change other parameter in searching area see next step.
### Function 'Paremeters.m'
In this function you can set up the search area for PSO and the starting epoch of the observation.

The **xmaxmin** variable is the search range of PSO, you can set it in this function.
