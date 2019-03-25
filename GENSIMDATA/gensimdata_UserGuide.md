This is a brief description of the main function, **gensimdata.m**, used for generating simulation data files. Use 'help gensimdata' in Matlab to see more details.

## Usage
GENSIMDATA(P,C,O)

P is the path to a .mat **input_parameter_file** containing the variables listed below, C is the path to a .mat **pulsar_catalog_file** containing information about the pulsars in the PTA, and O is the **output_folder** name where the data realizations will be stored.

### input_parameter_file
This is a .mat file containing the following variables: the locations of the GW sources, number of pulsars in the timing array, starting epoch of the observation, information about the noise realizations.

### pulsar_catalog_file
This file contains information about the pulsars in a Pulsar Timing Array (PTA). The file **survey_ska.mat** provides an example of a pulsar_catalog_file.

### Generating the parameter files
The Function **parameters.m** provides a template that can be modified to generate the input_parameter_file. This function generates two files, **Sim_Params_X.mat** and **searchParams_simDataSKA_X.mat**. The latter file contains parameters for Particle Swarm Optimization (PSO) that is used in the analysis codes (see the documentation in the **MxAvPhase** folder).

**Note**: Though these variables are already in the input_parameter_file but will still be individually stored in a file **searchParams_simDataSKA_X.mat** for the future use of multiple sources.

The input arguments for **parameters.m** are listed below:

**Ns**

```matlab
number of GW sources.
```

**Np**

```matlab
number of pulsars.
```

**Nrlz**

```matlab
number of noise realizations H1
```

**Nnis**

```matlab
number of realization of noise only cases H0
```

### Subfunctions needed by GENSIMDATA

* GenerateRandomGWSource, SpherePointPicking: generates a random realization of Ns SMBHBs with parameters, such as chirp mass (Mc), location (alpha, delta), distance (r), GW frequency (fgw), angles (iota, Psi, Phi0), sampled either uniformly or log uniformly. These parameters are saved in H1 data file. The parameter values are set to 0 for H0 data.

* coco, cosine, rotm_coo: coordinates transformation

* likelihood, InnProduct,LLR_PSOmpp, [Do not use etc.... all functions needed by gensimdata should be listed] etc:  calculating the time series of timing residual and the likelihood.

### Example
See **test_gensimdata.m** for a worked out example of how to generate simulated data files.
