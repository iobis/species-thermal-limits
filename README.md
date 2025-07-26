# Thermal limits of marine species based on occurrence data

Experimenting methods to estimate species' thermal limits based on occurrence data.

- [Overview](#overview)
- [Models](#models)
- [Auxiliary code](#auxiliary-code)
- [Deployment](#deployment)
- [Contributors](#contributors)
- [Action plan](#action-plan)

## Overview

Temperature is a key factor influencing the distribution and abundance of marine species, as it directly affects physiological processes. Organism activity (growth, reproduction, and survival) is optimum at a specific temperature, decreasing towards the extremes, up to a certain point where it becomes lethal (thermal limits). Understanding these limits is crucial for predicting species responses to climate change and managing marine ecosystems.

Although not universal, performance curves are often bell-shaped, following a Gaussian distribution. It is possible to estimate the parameters of these curves using occurrence data, which can be used to predict species' thermal limits. This is usually done by considering the 95th or 99th quantile of the temperatures at which the species was recorded. However, this approach is specially challenging for species that occur close to the maximum (or minimum) temperatures available in the environment (this is the case of tropical species), since the upper temperature of occurrence no longer reflects the species' thermal limit, but rather the maximum temperature available in the environment. This project aims to develop a method to estimate species' thermal limits based on occurrence data, using occupancy models. The method is developed using hierarchical models under a Bayesian framework.

## Models

Models were developed with increasing complexity, to enable fully testing the method with simulations. Each model is named according to its order of development (starting with `codes/model1.stan`), but some models might have been removed (e.g. `model3.stan`). Then, each model is accompanied by a test script, named accordingly (e.g. `codes/model1_tests.R`).

### Model 1

**Description**

The first model is a simple occupancy model.

Let:

- $y_i \in \{0, 1\}$ be the observed detection at survey $i$  
- $\text{sid}[i] \in \{1, \dots, N_{\text{spp}} \}$ be the species ID at survey $i$  
- $\text{sst}[i]$ be the sea surface temperature at survey $i$

Then we have parameters:

- $\mu_j \sim \mathcal{N}(20, 3)$ — thermal optimum (mean) for species $j$  
- $\sigma_j \sim \mathcal{N}(5, 1)$ — thermal niche breadth (standard deviation) for species $j$  
- $p \sim \text{Beta}(2, 2)$ — detection probability (shared across species)

The model is defined as:

For each survey $i$, define the suitability-based occupancy probability:

$$
q_i = \exp\left( -\frac{1}{2} \left( \frac{\text{sst}_i - \mu_j}{\sigma_j} \right)^2 \right)
\quad \text{where } j = \text{sid}[i]
$$

This defines the probability of occupancy using a Gaussian distribution.

The detection model accounts for both occupancy and imperfect detection, that is, models the probability that the species is either truly absent, or present but undetected.

- If $y_i = 1$:

$$
\log P(y_i = 1) = \log(p) + \log(q_i)
$$

- If $y_i = 0$:

$$
\log P(y_i = 0) = \log\left( (1 - q_i) + q_i (1 - p) \right)
$$

**What we tested**

- Model functionality with simulated data under normal conditions
- Model functionality with simulated data under extreme conditions (i.e. truncated distribution)
- Varying number of occurrences and absences
- Real species (fish data)

### Model 2

**Description**

This model is written exactly like Model 1, but it enables to pass different priors for the $\mu_j$ and $\sigma_j$ parameters.

**What we tested**

The effect of different priors on the model performance for those two previous tests:

- Model functionality with simulated data under normal conditions
- Model functionality with simulated data under extreme conditions (i.e. truncated distribution)

We also do the test with real species, but is not the focus of this model.

### Model 4

**Description**

Model 4 modified model 1 to now express the occupancy suitability in the log scale:

$$
\log q_i = -\frac{1}{2} \left( \frac{\text{sst}_i - \mu_j}{\sigma_j} \right)^2
\quad \text{where } j = \text{sid}[i]
$$

With consequent change in the detection model for $y_i = 1$:

$$
\log P(y_i = 1) = \log(p) + q_i
$$

It improves numerical stability. 

**What we tested**

Same tests as Model 1.

### Model 5

It extends Model 4 by adding a new parameter to the detection model, which defines a maximum occupancy probability for each species. This is to accommodate the possibility that the species might be rare.

**Description**

It adds the parameter $\theta_j \sim \mathcal{Beta}(5, 1)$.

Modifying the occupancy suitability formula to:

$$
\log q_i = -\frac{1}{2} \left( \frac{\text{sst}_i - \mu_j}{\sigma_j} \right)^2 + \log(\theta_j)
\quad \text{where } j = \text{sid}[i]
$$

**What we tested**

Same tests as Model 1, but now the simulation accounts for different maximum occupancy probabilities per species.

### Model 6

Model 6 adds partial pooling to Model 5 (a varying effects model). That way we estimate the variation between species through a global mean and standard deviation for the $\mu_j$ and $\sigma_j$ parameters, which are then adjusted for each species to obtain their unique optimum and tolerance. More specifically, this is the new model:

We model species-level thermal preferences hierarchically:

$\text{log tmu} \_j \sim \mathcal{N}(\mu_{\text{tmu}}, \sigma_{\text{tmu}})$ [Temperature optimum for species $j$]  
$\text{log tsd} \_j \sim \mathcal{N}(\mu_{\text{tsd}}, \sigma_{\text{tsd}})$ [Thermal tolerance (SD) for species $j$]  
$\theta \_j \sim \text{Beta}(5, 1)$ [Maximum occupancy probability for species $j$]  
$p \sim \text{Beta}(2, 2)$ [Detection probability]  

Where:

- $\mu_{\text{tmu}} \sim \text{LogNormal}(20, 3)$ [Global mean thermal optimum]  
- $\mu_{\text{tsd}} \sim \text{LogNormal}(5, 1)$ [Global mean thermal tolerance (SD)]  
- $\sigma_{\text{tmu}}, \sigma_{\text{tsd}} \sim \text{Exponential}(1)$ [Standard deviations for the global means]  

We interpret:

- $\exp(\text{log tmu} \_j)$ as the thermal optimum for species $j$
- $\exp(\text{log tsd} \_j)$ as the thermal tolerance (SD) for species $j$
- $\theta_j = \text{tomax}_j$ as the maximum occupancy probability

For survey $i$ and species $j = \text{sid}[i]$, define:

``` math
\log q_i = -\frac{1}{2} \left( \frac{\text{sst}_i - \exp(\text{log\_tmu}_j)}{\exp(\text{log\_tsd}_j)} \right)^2 + \log(\theta_j)
```

**What we tested**

- Model functionality with simulated data under normal conditions
- Model functionality with simulated data under extreme conditions (i.e. truncated distribution)
- Varying number of occurrences and absences, and the effect of low absences on pushing the species optimum towards the global mean
- Real species (fish data)

Note that here we explicitly compare the results of Model 6 with Model 5, to check the effect of partial pooling on the model performance.

### Model 7

Model 7 also uses partial pooling, but extends the previous one by modeling the species-level parameters jointly using a multivariate normal distribution with a full correlation structure.

More formally we have this structure:

Considering the following parameters:

- $\mathbf{z}_j \sim \mathcal{N}(\mathbf{0}, \mathbf{I}_3)$ — standard normal latent variables  
- $\boldsymbol{\Sigma} = \text{diag}(\boldsymbol{\sigma}) \cdot \mathbf{R} \cdot \text{diag}(\boldsymbol{\sigma})$ — full covariance matrix  
- $\mathbf{L}_{\mathbf{R}} \sim \text{LKJ}(4)$ — Cholesky factor of correlation matrix $\mathbf{R}$ 
- $\boldsymbol{\mu} = \text{spp}\_{means} \in \mathbb{R}^3$

For each species $j$, we compute:

``` math
\begin{aligned}
\mathbf{v}_j &= \boldsymbol{\mu} + \mathbf{L}_{\boldsymbol{\Sigma}} \mathbf{z}_j \\
\log \mu_j &= v_{j,1} \\
\log \sigma_j &= v_{j,2} \\
\theta_j &= \text{logit}^{-1}(v_{j,3})
\end{aligned}
```

Where:

- $\mu_j$ is the species mean temperature (optimum)
- $\sigma_j$ is the species temperature standard deviation (tolerance)
- $\theta_j$ is the maximum occupancy probability

**What we tested**

- Model functionality with simulated data under normal conditions
- Model functionality with simulated data under extreme conditions (i.e. truncated distribution)
- Varying number of occurrences and absences, and the effect of low absences on pushing the species optimum towards the global mean
- Real species (fish data)

Note that, again, here we explicitly compare the results of Model 7 with Model 5, to check the effect of partial pooling on the model performance.

## Auxiliary code

Codes used for preparing non-simulated data are on `codes/data-processing`.

## Deployment

If you want to run the models, you need to install `Stan` and `rstan`. You can found instructions [here](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).

You will also need the [`rethinking`](https://github.com/rmcelreath/rethinking) R package, which can be installed from GitHub:

```r
devtools::install_github("rmcelreath/rethinking")
```

Check all dependencies in the `r-dependencies.R` file.

## Contributors

- Richard McElreath - model development, testing and documentation
- Silas Principe - testing and documentation
- Pieter Provoost - testing and documentation

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
