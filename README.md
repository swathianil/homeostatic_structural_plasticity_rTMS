# Repetitive transcranial magnetic stimulation (rTMS) triggers dose-dependent homeostatic rewiring in recurrent neuronal networks

This repository consists of HPC-ready version of codes for Project Network. The project aims to apply tetanic electrical current injection-like stimulation onto a recurrent neuron network, that follows structural plasticity rules based on firing-rate homeostasis. 

## Instructions
Every condition consists of a package of data collection and analysis scripts. These include:

MAIN:
 - *.sh shell script to run the condition through stages of growth, stimulation, post-stimulation and decay. (1)
 - *an.sh shell script to run analyses on the four stages, split into four different *.moab job calls (2)

PERIPHERAL: 
 - *.py data collection script (eg.,100Hz80kpA.py) (3)
 - Four *.py analyses scripts  (eg., 100gr.py, 100stim.py, 100post.py, 100decay.py) (4)
 - *.yaml files for importing parameters: simulation.yaml, Experiment.yaml, network.yaml
 - *.moab script for jobcall assistance 

 In order to run a condition, run the respective (1) *.sh file.
 If you need to customise parameters (except amplitude), you can do so in the (3) *.py script, after the *.yaml files are unpacked. Be sure to then update it on all the analyses (4) *.py files as well.
 The amplitude is passed in via the shell script and it can be modified there. Be sure to change the job name as well.
 Parameters with absolute values are in the *.yaml files. Others that are deriviatives are within (1)
 
 ## Resource and requirements
 - Python v3.x
 - NEST v 2.2
 - yaml (available with Nemo python module)
 -Requested resources (relevant to HPC users):
  - 5 nodes* 20 procs
  - Walltime ~ 20h for data collection
  
