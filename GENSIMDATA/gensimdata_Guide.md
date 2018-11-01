This is a brief description of the main function, **gensimdata.m**, used for generating simulation data files.
## Usage
GENSIMDATA(P,C,O)
P is the path to a .mat **input_parameter_file** containing the variables listed below, C is the path to a .mat **pulsar_catalog_file** containing information about the pulsars in the PTA, and O is the **output_folder** name where the data realizations will be stored.

### input_parameter_file
The Function **parameters.m** can be used to generate the input_parameter_file. This function generates two files, **parameter.mat** and **searchParams_simDataSKA_X.m**. The .mat file **parameter.mat** contains the variables: the location of random GW sources, number of pulsars in the timing array, starting epoch of the observation, information about the noise realizations. The search range for Particle Swarm Optimization (PSO) is stored separately in the file **searchParams_simDataSKA_X.mat** for running on LS5(LoneStar 5).

**Note**: Though these variables are already in the input_parameter_file but will still be individually stored in a file **searchParams_simDataSKA_X.mat** for the future use of multiple sources.

The input arguments for parameters.m are listed below:

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

### pulsar_catalog_file
This file contains information about the pulsars in a Pulsar Timing Array (PTA). The file **survey_ska.mat** provides an example of a pulsar_catalog_file.

### Subfunctions needed by GENSIMDATA

* GenerateRandomGWSource, SpherePointPicking: generates a random realization of Ns SMBHBs with parameters, such as chirp mass (Mc), location (alpha, delta), distance (r), GW frequency (fgw), angles (iota, Psi, Phi0), sampled either uniformly or log uniformly. These parameters are saved in H1 data file. The parameter values are set to 0 for H0 data.

* coco, cosine, rotm_coo: coordinates transformation

* likelihood, InnProduct,LLR_PSOmpp, [Do not use etc.... all functions needed by gensimdata should be listed] etc:  calculating the time series of timing residual and the likelihood.

### Example
See **test_gensimdata.m** for a worked out example of how to generate simulated data files.
