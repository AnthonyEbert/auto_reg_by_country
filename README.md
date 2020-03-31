# Bayesian analysis of the COVID-19 pandemic

Bayesian analysis of COVID-19 case time-series for individual regions using a simulation based model and Approximate Bayesian Computation (ABC).

## Authors

1. Christopher Drovandi (c.drovandi@qut.edu.au),
                School of Mathematical Sciences, 
                Science and Engineering Faculty, 
                Queensland Univeristy of Technology 
Google Scholar: (https://scholar.google.com.au/citations?user=jc0284YAAAAJ&hl=en)

2. David J. Warne (david.warne@qut.edu.au),
                School of Mathematical Sciences, 
                Science and Engineering Faculty, 
                Queensland Univeristy of Technology 
Google Scholar: (https://scholar.google.com.au/citations?user=t8l-kuoAAAAJ&hl=en)

## Model summary
The model is a discrete state continuous-time Markov Process model that is based on an SIR model with important extentions. The form of the interactions are:
$$S + I \overset{\alpha_0 + \alpha f(C)}{\rightarrow} 2I,\quad C \overset{\beta}{\rightarrow} R,\quad I\overset{\gamma}{\rightarrow} C, \quad C\overset{\delta}{\rightarrow} D,\quad\text{and}\quad I\overset{\eta\beta}{\rightarrow} R^u$$
where $S$ is susceptible, $I$ is infected (latent variable), $C$ are confirmed cases, $R$ are case recoveries, $D$ are case fatalities, and $R^u$ are removed (recoveries of deaths of infected latent variable). Parameters are the infection rates $\alpha_0$ and $\alpha$, case recovery rate $\beta$, case fatality rate $\delta$, case detection rate $\gamma$, and overal removal rate relative to case recovery $\eta$. The function $f(C)$ represents a regulatory response of the community (e.g., social distancing etc...) we assume the form $f(C) = \dfrac{1}{1+C^{n}}$ where $n \in [0,1]$. Inital contion for $I$ is given by $I_0 = \kappa C_0$. 

## Functions and scripts

The followinf files are provided:
* `run_smc_intensity_reg.m` runs analysis pipeline (computes posteriors and samples prior predictive) given a country ID (number from 1 to 252).

* `simuldata_reg.m` forwards stochastic simulation of model given a vector of parameters.

* `smry.m` summary statistic function for usage in the ABC method.

* `smc_abc_rw.m` Adaptive sequential Monte Carlo sampler for ABC. This should not need to be modified even if the model completely changes.

* `TauLeapingMethod.m` routine for approximate stochastic simulation of discrete-state continous-time Markov process (used in model simulation).

* `import_covid19data.m` function for importing COVID-19 data.

* `import_JHUtoA3.m` function to import conversion table from Johns-Hopkins University country/region names to ISO-3166 alpha-3 codes.

* `import_country_full_list.m` function to import country information by ISO-3166 alpha-3 codes. For example, population data.

## Usage

1. Edit `run_smc_intensity_reg.m` to ensure  the variable `DATA_DIR` is pointing to a directory that contains clones of the two github repos `COVID16data` and `COVID-19_ISO-3166`.
2. Define a country/region id in the MATLAB workspace (i.e., the row number in the population table for the region of interest), e.g., for China
use `country_id = 44` or for Italy use `country_id = 114`.
3. If required edit the smc parameters, defaults are quite reasonable for the moment. Please contact authors if there are any difficulties editing these.
4. Run the main script `run_smc_intensity_reg`. This may take a while to run. Currently, the code is not parallelised to ensure reproducibility (RNG in parallel are not reproducible). If Parallel SMC is desired, then edit `run_smc_intensity_reg.m` to utilise `smc_abc_rw_par.m`.


