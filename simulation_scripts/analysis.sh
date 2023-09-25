: '
Shell script for job submission to HPC (MOAB compatible)

        Task: 
            Analyse data obtained via jobcall.sh through 4 phases: growth, stim, decay, post

        Parameters:
            amplitude (float)    : rTMS stimulation strength in pA
            whereami  (string)   : fetch current location (all scripts must be stored here)
            stim      (string)   : descriptive folder name with nomenclature {Frequency_in_Hz}_{Aplitude_in_nA}_{number_of_pulses}p_seed{experiment_seed}
            seed      (int)      : Experiment seed. Other relevant seeds will be derived from thisvalue
            experiment(string)   : Name of *.yaml file contanining experiment parameters
            anscript  (string)   : Name of analysis script, *.py file 
            phase     (string)   : analysis phase (growth/stim/post/decay)
        
        Output:
            *.npy files distributed in 'results' directory, results from all phases are consolidated

        Example:
            msub analysis.moab -v amplitude=30000. -v whereami="$PWD" -v stim='10Hz_30nA_10pc_9000p_seed1' -v seed=1 -v experiment='Experiment.yaml' -v anscript='analysis.py' -v phase='growth'  -N 10Hz_30nA_10pc_9000p_seed1_gr

'
amplitude=34000.
stim="10Hz_30nA_10pc_9000p_seed1"
seed=1
experiment="Experiment.yaml"

msub analysis.moab -v amplitude=$amplitude -v whereami="$PWD" -v stim=$stim -v seed=$seed -v experiment=$experiment -v anscript='analysis.py' -v phase='growth'  -N  $stim:gr

msub analysis.moab -v amplitude=$amplitude -v whereami="$PWD" -v stim=$stim -v seed=$seed -v experiment=$experiment -v anscript='analysis.py' -v phase='stim'  -N $stim:stim
 
msub analysis.moab -v amplitude=$amplitude -v whereami="$PWD" -v stim=$stim -v seed=$seed -v experiment=$experiment -v anscript='analysis.py' -v phase='post'  -N $stim:post
 
msub analysis.moab -v amplitude=$amplitude -v whereami="$PWD" -v stim=$stim -v seed=$seed -v experiment=$experiment -v anscript='analysis.py' -v phase='decay'  -N  $stim:decay

