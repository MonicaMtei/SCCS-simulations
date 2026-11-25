# SCCS simulations

# Project Title: When do pre-exposure periods work in the self-controlled case series (SCCS)? A simulation study investigating event-dependent exposure bias

## Overview
This project investigates the performance of different analytical approaches to address event-dependent exposure bias in self-controlled case series when evaluating the safety outcomes

## Simulation design
The simulation comprises of 8 base case scenarions (a combination of 2 relative incidences 1 and 2, 2 observation period lengths of 1 and 4 years and 2 vaccine uptake distribution patterns namely uniform and left skewed)  in which there is no event-dependent exposure. Each base case was varied tp introduced increasing extent of event-dependent exposures from temporary delays (variation_id 1-4) and permanent vaccine differals (variation_id 5&6)
Each base case and its variation was simulated 5000 times and the following statistical approaches were applied
1)The standard SCCS model
2)Pre-exposure period of 14 days included in the standard SCCS model
3)Pre-exposure period of 28 days included in the standard SCCS model
4)Pre-exposure period of 60 days included in the standard SCCS model
5)Modified SCCS model for event-dependent exposures

## Performance measures
The approaches were evaluated using bias, empirical standard error, average-model based standard error, mean square error and coverage


### Repository structure

The R/ directory contains the following scripts:

| Script | Description |
|---|---|
| `01_libraries.R` | Required R packages |
| `02_seeds.R` | Seed generation for reproducibility |
| `03_simulate_base_skewed.R` | Simulation of base case with left-skewed vaccine uptake |
| `04_simulate_base_uniform.R` | Simulation of base case with uniform vaccine uptake |
| `05_variations_to_each_base.R` | Introduce variations to each base case |
| `06_fit_models_and_extract_estimates.R` | Fitting SCCS models and extracting estimates |
| `07_run_replicates.R` | Running simulation replicates |
| `08_performance_measures.R` | Calculation of simulation performance measures |


### Simulation results

The `results/` directory contains the individual simulation
replicate results and summarised performance measures.

### Replicate results

The `results/replicate_results/` directory contains 56 CSV
files, one for each simulation scenario, with results from
5,000 replicates per scenario and their variations.

These files allow users to reproduce the performance
measures without rerunning the simulations.

### Performance measures

The `results/performance_measures/` directory contains
the summarised performance measures across the simulation
scenarios, doses and analytical approaches.

Performance measures include:

- Bias
- Empirical standard error
- Average model-based standard error
- Mean squared error
- Coverage
- Monte Carlo standard errors for each performance measure
