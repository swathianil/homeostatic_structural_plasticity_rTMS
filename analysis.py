# Dependencies
import numpy as np
import yaml
import nest
import sys

# Accept and assign variables
scriptpath = sys.argv[1]  # Path to Scripts
scriptpath += "/"
Experiment = sys.argv[2]  # Experiment defining file

# Load Parameters and unpack
yaml_files = ["network.yaml", "simulation.yaml", Experiment]
for file in yaml_files:
    E = scriptpath + file
    with open(E) as file:
        parameters = yaml.full_load(file)
    locals().update(parameters)

# Accept and assign remaining variables
amp = sys.argv[3]  # Raw value of amplitude
stim = str(sys.argv[4])  # stimulation name
arrayindex = int(sys.argv[5])  # seed
phase = str(sys.argv[6])
datapath = scriptpath + "data/" + (stim) + "/" + phase + "/"
TMS_amplitude = float(amp)

# Construct resultpath from stim nomenclature
pos = stim.find("seed") - 1
resname = str(stim[:pos])
resultspath = scriptpath + "results/" + resname + "/" + "seed_" + str(arrayindex) + "_"

# Neuron & network parameters
n_rank = total_num_virtual_procs * 1
NE = order * 4
NI = order * 1
N = NE + NI
stimulated_pop = int(stimulated_fraction * NE)

# Stimulation settings: TMS parameters
if TMS_protocol == "simple_rTMS":
    stimulation_time = 1000 * (n_pulses / TMS_stim_frequency)
    stimulation_cicles = stimulation_time / stimulation_cicle_dur

# Define simulation step size
growth_step = growth_time / cicles
stimulation_step = stimulation_time / stimulation_cicles
post_stimulation_step = post_stimulation_time / post_stimulation_cicles
decay_step = decay_time / decay_cicles

# Create simulation time steps
stimulation_end = growth_time + stimulation_time
post_stimulation_end = stimulation_end + post_stimulation_time
decay_end = post_stimulation_end + decay_time

growth_steps = np.arange(growth_step, growth_time + 1, growth_step)
stimulation_steps = np.arange(
    growth_time + stimulation_step, stimulation_end + 1, stimulation_step
)
post_stimulation_steps = np.arange(
    stimulation_end + post_stimulation_step,
    post_stimulation_end + 1,
    post_stimulation_step,
)
decay_steps = np.arange(post_stimulation_end + decay_step, decay_end + 1, decay_step)

all_steps = np.concatenate(
    ([0], growth_steps, stimulation_steps, post_stimulation_steps, decay_steps)
)

# Initiate containers for results
firing_rate = np.zeros([3, len(all_steps)])
connectivity = np.zeros((2, 2, len(all_steps)))
rate_stim = []
rate_exc = []
rate_inh = []
rates_all = []

# Functions
# Rate calculator
def analysis(senders, times, simstep_dur, N):
    rates = []
    for i in np.arange(0, N, 1):
        idx = np.where(senders == i + 1)
        time = times[idx]
        rate = len(time) / (simstep_dur / 1000.0)
        rates.append(rate)
    return rates


# Store phase relevant data: Positional index for Tensor, simulations steps ,step size
phase_data = {
    "growth": [1, growth_steps, growth_step],
    "stim": [cicles + 1, stimulation_steps, stimulation_step],
    "post": [
        cicles + len(stimulation_steps) + 1,
        post_stimulation_steps,
        post_stimulation_step,
    ],
    "decay": [
        (cicles + len(stimulation_steps) + post_stimulation_cicles + 1),
        decay_steps,
        decay_step,
    ],
}

# Extract phase relevant data
aux = phase_data[phase][0]
steplist = phase_data[phase][1]
stepsize = phase_data[phase][2]

# Analysis
for step in steplist:
    if step == 0:
        continue
    senders = np.array([])
    times = np.array([])
    sources = np.array([])
    targets = np.array([])
    events = np.array([])
    local_connections = np.array([])

    for rank in np.arange(n_rank):
        extension = str(step) + "_rank_" + str(rank) + ".npy"
        # 
        senders = np.concatenate(
            (senders, np.load(datapath + "senders_" + extension)), axis=0
        )
        times = np.concatenate(
            (times, np.load(datapath + "times_" + extension)), axis=0
        )
        sources = np.concatenate(
            (sources, np.load(datapath + "sources_" + extension)), axis=0
        )
        targets = np.concatenate(
            (targets, np.load(datapath + "targets_" + extension)), axis=0
        )

    matrix = np.zeros((NE, NE)) # Initiate matrix for plastic connections (E-E)
    for ii in np.arange(sources.shape[0]):
        matrix[int(targets[ii] - 1), int(sources[ii] - 1)] += 1
    connectivity[0, 0, aux] = np.mean(matrix[:stimulated_pop, :stimulated_pop])
    connectivity[0, 1, aux] = np.mean(matrix[:stimulated_pop, stimulated_pop:])
    connectivity[1, 0, aux] = np.mean(matrix[stimulated_pop:, :stimulated_pop])
    connectivity[1, 1, aux] = np.mean(matrix[stimulated_pop:, stimulated_pop:])

    rates = analysis(senders, times, stepsize, N)
    rate_stim.append(rates[0:stimulated_pop])
    rate_exc.append(rates[stimulated_pop:NE])
    rate_inh.append(rates[NE:])
    rates_all.append(rates)

    spike_count = np.histogram(senders, bins=np.array([0, stimulated_pop, NE, N]))[0]
    firing_rate[:, aux] = (
        spike_count
        / np.array([stimulated_pop, NE - stimulated_pop, NI]).astype(float)
        / stepsize
        * 1000.0
    )

    aux += 1


np.save(resultspath + phase + "_Con", connectivity)

np.save(resultspath + phase + "_FR", firing_rate)

np.save(resultspath + "All_Steps", all_steps)

np.save(resultspath + phase + "_Steps", steplist)

np.save(resultspath + phase + "_rate_stim", rate_stim)

np.save(resultspath + phase + "_rate_exc", rate_exc)

np.save(resultspath + phase + "_rate_inh", rate_inh)

np.save(resultspath + phase + "_rate_all", rates_all)


def checkpoint(resultspath, phase):
    filename = resultspath + phase + "_AnalysisComplete.npy"
    x = 1
    np.save(filename, x)


checkpoint(resultspath, phase)
