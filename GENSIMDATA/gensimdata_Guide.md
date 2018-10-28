# Function 'gensimdata.m'
This is the main function to use for generating simulation data files.
## How to use
The input arguments for **gensimdata.m** are: The path to a .mat input parameter file containing the variables listed below, the path to a .mat file containing information about the pulsars in the PTA, and the output folder where you want to store the output files.

### Function 'Paremeters.m'
#### Call the Function
By calling the function **parameters** to generate your own parameter file with the  input arguments listed below:

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
#### Edit the function itself to set the range of PSO
In this function you can set up the search area for PSO and the starting epoch of the observation.

The **xmaxmin** variable is the search range of PSO, you can set it in this function.
### Variables required in file containing PTA information
This file is a pulsar catalog called **survey_ska.mat** contains the information of pulsars.

### Example
See the **test_gensimdata.m** file for a worked out example about how to generate the simulate data files.
