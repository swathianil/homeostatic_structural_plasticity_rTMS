# Repetitive transcranial magnetic stimulation (rTMS) triggers dose-dependent homeostatic rewiring in recurrent neuronal networks
Swathi Anil, Han Lu, Stefan Rotter, Andreas Vlachos

This repository consists of HPC-ready code used for the simulations in the paper.

The project aims to apply tetanic electrical current injection-like stimulation onto a recurrent neuron network, that follows structural plasticity rules based on firing-rate homeostasis. 
Simulations were carried out using bwForCluster NEMO and job scripts have been implemented to suit the MOAB workload manager. 

## Instructions
Each experiment (stimulation condition) consists of a package of scripts for data collection and analysis. These include:

### MAIN:
 - jobcall.sh shell script to run the condition through the network stages of growth, stimulation, post-stimulation and decay. (1)
 - analysis.sh shell script to run analyses on the four stages, excecuted via four different *.moab job calls 

### PERIPHERAL: 
 - master.py data collection script 
 - analysis.py analyses scripts 
 - *.yaml files for importing parameters: simulation.yaml, Experiment.yaml, network.yaml
 - jobcall.moab for job submission via jobcall.sh
 - analysis.moab for job submission via analysis.sh 

 In order to run a condition, populate and run the jobcall.sh script.
 If you need to customise parameters, you can do so in the Experiment.yaml file. Stimulation times can be accessed in simulation.yaml.
 Values in shell script override values in *.yaml

 ## Resource and requirements
 - Python v3.x
 - NEST v 2.2
 - yaml (available with Nemo python module)
 -Requested resources (relevant to bwForCluster NEMO HPC users):
  - 5 nodes* 20 procs
  - Walltime ~ 20h for data collection
  
## Data
Data supporting this preprint can be found at: 
https://doi.org/10.5281/zenodo.8374484
The jupyter notebook 'GraphPlotter' in this repository can be used to interactively regenerate the figures available in this preprint.

