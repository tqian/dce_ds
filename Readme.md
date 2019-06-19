# Deductive and Computerizable Estimation for Double-Sampling Designs

The code produces deductive estimate for survival probability in a double-sampling design. For methodology details, see working paper:

Qian T., and Frangakis C. Deductive Semiparametric Estimation in Double-Sampling Designs. https://arxiv.org/abs/1902.11147




## Code Structure

fcn_double-sampling_survival-prob.R: function definition.

main_double-sampling_survival-prob.R: code for estimating survival probability for PEPFAR study on HIV. Corresponding data set is omitted due to data security.

main_simulation.R: code for simulation.

main_simulation.sh: shell script for "qsub -t" to the cluster for parallelization.

misc/: random seeds for parallelization.
