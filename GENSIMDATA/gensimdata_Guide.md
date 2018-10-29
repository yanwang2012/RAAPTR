# Function 'gensimdata.m'
This is the main function to use for generating simulation data files.
## How to use
The input arguments for **gensimdata.m** are: The path to a .mat **input_parameter_file** containing the variables listed below, the path to a .mat **pulsar_catalog_file** containing information about the pulsars in the PTA, and the **output_folder** name where the data realizations will be stored. 

### Variables required in the input_parameter_file
[Describe.]

### Variables required pulsar_catalog_file
This file is a pulsar catalog called **survey_ska.mat** contains the information of pulsars.
[Describe the variables required in this file. Then point the reader to survey_ska.mat as an example pulsar_catalog_file.]

#### Edit the function itself to set the range of PSO 
[The user should not have to change any function... defeats the point of having functions. The code should be modified so that this is not required.]
In this function you can set up the search range for each parameter and the starting epoch of the observation.

The **xmaxmin** variable is the search range of PSO, you can set it in this function.
[In which file is xmaxmin stored and how does gensimdata.m use this file? This file is not mentioned above in the list of inputs to gensimdata. This variable is needed by the AvPhase/MaxPhase codes. I don't think gensimdata.m needs to know xmaxmin. Its description should be moved to a different place where the AvPhase/MaxPhase usage is described.] 


### Example
See **test_gensimdata.m** for a worked out example of how to generate simulated data files.

The input_parameter_file can be generated using the function **paremeters.m** by calling it with with the  input arguments listed below:

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
