: '
Shell script for job submission to HPC (MOAB compatible)

        Task: 
            Create folders for data and results
            Run experiment through 4 phases: growth, stim, post, decay

        Parameters:
            amplitude (float)  : rTMS stimulation strength in pA
            whereami  (string) : fetch current location (all scripts must be stored here)
            stim      (string) : descriptive folder name with nomenclature {Frequency_in_Hz}_{Aplitude_in_nA}_{number_of_pulses}p_seed{experiment_seed}
            seed      (int)    : Experiment seed. Other relevant seeds will be derived from thisvalue
            experiment(string) : Name of *.yaml file contanining experiment parameters
            script    (string) : Name of master script, *.py file used to run the experiment

        Prerequisite:
	    "logfiles" and "data" folders in current location. data folder must have following subfolders: growth, stim, post, decay
        
        Output:
            *.npy files distributed in "data" directory according to phase of experiment (growth, stim, post, decay, rTMS, post-rTMS) 
            log error file (*.e*) and output file (*.o*) in logfiles folder in current location

        Example:
            msub -e logfiles/ -o logfiles/ jobcall.moab -v amplitude=30000. -v whereami="$PWD" -v stim="10Hz_30nA_10pc_9000p_seed1" -v seed=1 -v experiment="Experiment.yaml" -v script="master.py" -N 10Hz_9000p

'
amplitude=34000.
whereami=$PWD
stim="test_seed1"
seed=1
experiment="Experiment.yaml"
resultsfolder="test"
name="testjob"

cp -r ./template ./data/$stim # make datadirectory
mkdir ./results/$resultsfolder # make resultsdirectory

# Submit job
msub -e logfiles/ -o logfiles/ jobcall.moab -v amplitude=$amplitude -v whereami="$PWD" -v stim=$stim -v seed=$seed -v experiment=$experiment -v script="master.py" -N  $name
