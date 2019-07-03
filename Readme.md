# Deductive and Computerizable Estimation for Double-Sampling Designs

The code produces deductive estimate for survival probability in a double-sampling design. For methodology details, see working paper:

Qian T., and Frangakis C. Deductive Semiparametric Estimation in Double-Sampling Designs. https://arxiv.org/abs/1902.11147

## Note

Note that all directory names in the R code need to be modified accordingly if running on another computer.

It will take about 13,000 hours on a single core CPU to run the entire data analysis.

It will take about 34,000 hours on a single core CPU to run the entire simulation.


## Code Structure

Note that all directory names in the R code need to be modified accordingly if running on another computer.

### simulation/: code for simulation

simulation/dgm.R: generative model

simulation/fcn_DE.R: wrapper functions for the deductive estimator

simulation/fcn_double-sampling_survival-prob.R: core functions for the deductive estimator

simulation/fcn_weighted-KM.R: functions for the Kaplan-Meier estimators

simulation/make_tables.R: make tables for the simulation section in the paper

simulation/parallel_seeds.csv: random seeds for parallel jobs

simulation/simulation_parallel_jhpce.R: main file to conduct simulation (R code)

simulation/simulation_parallel_jhpce.sh: main file to conduct simulation (shell script)

simulation/table_generation/: utility code helping the table generation in the paper

### data analysis/: code for data analysis

data analysis/data_analysis_for_paper_harvardrc.R: main file to compute deductive estimators for data analysis (R code)

data analysis/data_analysis_for_paper_harvardrc.slurm: main file to compute deductive estimators for data analysis (slurm file)

data analysis/data_analysis_for_paper_KmAllMethods.R: main file to compute Kaplan-Meier estimators for data analysis

data analysis/fcn_double-sampling_survival-prob.R: core functions for the deductive estimator

data analysis/fcn_weighted-KM.R: functions for the Kaplan-Meier estimators

data analysis/make_plot_and_table.R: make tables and figures for the data analysis section in the paper

data analysis/parallel_seeds.csv: random seeds for parallel jobs

data analysis/table_generation/: utility code helping the table generation in the paper
