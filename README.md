# Thermal limits of marine species based on occurrence data

Experimenting methods to estimate species' thermal limits based on occurrence data.

## Method tests

We are testing the proposed method with 3 distinct conditions:

1. Simulated data
   - Simulated suitability based on a virtual surface
   - Simulated suitability based on a real surface (i.e. SST data from a location)
2. True occurrence data with absence information
3. True occurrence data with inferred absence information

The functions to simulate the species are on the script `simulate_species.R`. You can also generate multiple datasets by using `simulate_datasets.R`

Codes for model testing are named [type]_tests.R. For example, the code for the simulated data is `simulation_tests.R`.

The current [Stan](https://mc-stan.org/) model in use is **model 4** (`model4.stan`)

## Action plan
Last modified in September 2024

1. Adjust simulation codes
2. Adjust models based on simulated data (trying skew normal, boundaries, etc.)
3. Test again the model with the simulated data
	1. Ideally, try the gaussian vs skew normal (?)
4. Select a subset of species for which we have good absence data
	1. That would be fish species, but maybe other benthic species covered by the Reef Life Survey
5. Select a subset of species for which we have experimental data
	1. Use GlobTherm database
6. Apply the model to this subset of species
7. Analyse effectiveness
8. Apply the model to the same species but using the month/depth matched data (OBISTherm dataset, preparing)
9. Analyse effectiveness and differences
10. Next steps
	1. Apply this to all OBIS data and make available in the website -> this is our main target, to support some of the projects (e.g. eDNA expeditions, PacMAN, MPA Europe)
	2. Try other things
		1. Models combining phylogenetic information (multiple species within the same group or from related groups)
  		2. Other ideas? 
	3. Publication of the method (?)
