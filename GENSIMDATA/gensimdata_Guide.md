# Function 'gensimdata.m'
This is the main function to generate simulation data about multiple source GW signals.
## How to use
Use the 'test_gensimdata.m' file to generate data files.

### Parameters in 'test_gensimdata.m'
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
You choose these Parameters in the function **test_gensimdata.m** then it will generate a parameter file named **parameter.mat**. In this file, you have all the parameters you need to simulate the data set. To change other parameter in searching area see next step.
### Function 'Paremeters.m'
In this function you can set up the search area for PSO and the starting epoch of the observation.

The **xmaxmin** variable is the search range of PSO, you can set it in this function.
