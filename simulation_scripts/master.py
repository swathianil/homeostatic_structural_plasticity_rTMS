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
datapath = scriptpath + "data/" + (stim) + "/"
TMS_amplitude = float(amp)

# Seed assignment
rank = nest.Rank()
seed = 7928 + (1000 * arrayindex)  # hardwired arbitrarily (!)
grng_seed = seed
rng_seeds = range(seed + 1, seed + total_num_virtual_procs + 1)

# Neuron & network parameters
n_rank = total_num_virtual_procs * 1  # n_rank will be used in job parsing by NEST MPI
NE = 4 * order  # number of excitatory neurons
NI = 1 * order  # number of inhibitory neurons
N = NE + NI  # total number of neurons
CE = int(eps * NE)  # number of incoming excitatory synapses per inhibitory neuron
CI = int(eps * NI)  # number of incoming inhibitory synapses per neuron
beta_Ca = 1.0 / tau_Ca  # increment on calcium trace per spike (1/ms)

neuron_params = {
    "C_m": CMem,
    "tau_m": tauMem,
    "t_ref": t_ref,
    "E_L": E_L,
    "V_reset": V_reset,
    "V_m": V_m,
    "beta_Ca": beta_Ca,
    "tau_Ca": tau_Ca,
    "V_th": theta,
}

nu_th = theta / (J * CE * tauMem)
nu_ex = eta * nu_th
rate = 1000.0 * nu_ex * CE  # Poisson background input
weight = J
synapse_params = {
    "min_delay": min_delay,
    "max_delay": max_delay,
}
stimulated_pop = int(stimulated_fraction * NE) # Subset of excitatory population receiving rTMS

# Stimulation settings: TMS parameters
if TMS_protocol == "simple_rTMS":
    stimulation_time = 1000 * (n_pulses / TMS_stim_frequency)  # in ms
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

# Setup Base network
# Configure kernel
nest.ResetKernel()
nest.EnableStructuralPlasticity()
nest.SetKernelStatus({"resolution": dt, "print_time": False})
nest.SetKernelStatus(
    {
        "structural_plasticity_update_interval": int(MSP_update_interval / dt),
        "grng_seed": grng_seed,
        "rng_seeds": rng_seeds,
        "total_num_virtual_procs": int(total_num_virtual_procs),
    }
)

# Create generic neuron with Axon and Dendrite
nest.SetDefaults(neuron_model, neuron_params)
nest.CopyModel(neuron_model, "excitatory")
nest.CopyModel(neuron_model, "inhibitory")
nest.CopyModel(
    "static_synapse", "inhibitory_synapse", {"weight": -g * weight, "delay": max_delay}
)
nest.CopyModel("static_synapse", "EI_synapse", {"weight": weight, "delay": max_delay})
nest.CopyModel(synapse_model, "msp_excitatory")
nest.CopyModel("static_synapse", "device", {"weight": weight, "delay": max_delay})
nest.SetDefaults("msp_excitatory", {"weight": weight, "delay": max_delay})

# Assign synaptic elements with growth curve to excitatory neuron model
gc_den = {
    "growth_curve": growth_curve_d,
    "z": z0_mean,
    "growth_rate": -slope * target_rate,
    "eps": target_rate,
    "continuous": False,
}
gc_axon = {
    "growth_curve": growth_curve_a,
    "z": z0_mean,
    "growth_rate": -slope * target_rate,
    "eps": target_rate,
    "continuous": False,
}
nest.SetDefaults(
    "excitatory", "synaptic_elements", {"Axon_exc": gc_axon, "Den_exc": gc_den}
)

# Use SetKernelStatus to activate the plastic synapses
nest.SetKernelStatus(
    {
        "min_delay": min_delay,
        "max_delay": max_delay,
        "structural_plasticity_synapses": {
            "syn1": {
                "model": "msp_excitatory",
                "weight": weight,
                "post_synaptic_element": "Den_exc",
                "pre_synaptic_element": "Axon_exc",
            }
        },
        "autapses": False,
    }
)

# Create & configure network nodes
pop_exc = nest.Create("excitatory", NE)
pop_inh = nest.Create("inhibitory", NI)
poisson_generator_inh = nest.Create("poisson_generator")
poisson_generator_ex = nest.Create("poisson_generator")
step_current_generator = nest.Create("step_current_generator")
spike_detector = nest.Create("spike_detector")
nest.SetStatus(poisson_generator_ex, {"rate": rate})
nest.SetStatus(poisson_generator_inh, {"rate": rate})
nest.SetStatus(spike_detector, {"withtime": True, "withgid": True})

# Connect devices
nest.Connect(poisson_generator_inh, pop_inh, "all_to_all", model="device")
nest.Connect(poisson_generator_ex, pop_exc, "all_to_all", model="device")
nest.Connect(pop_exc + pop_inh, spike_detector, "all_to_all", model="device")
nest.Connect(
    step_current_generator, pop_exc[:stimulated_pop], "all_to_all", model="device"
)

# Establish static connections
nest.Connect(pop_exc, pop_inh, {"rule": "fixed_indegree", "indegree": CE}, "EI_synapse")
nest.Connect(
    pop_inh,
    pop_exc + pop_inh,
    {"rule": "fixed_indegree", "indegree": CI},
    "inhibitory_synapse",
)

# Configure stimulator: step_current_generator
stimstart = growth_steps[-1]
timestamps = np.linspace(
    stimstart,
    stimstart + stimulation_time - (1000.0 / TMS_stim_frequency),
    n_pulses + 1,
).astype(int)
start_list = timestamps
stop_list = start_list + TMS_pulse_duration
steplist = np.append(start_list, stop_list)
steplist.sort()
amplist = []
for ind, val in enumerate(steplist):
    if ind % 2 == 0:
        amplist.insert(ind, TMS_amplitude)
    else:
        amplist.insert(ind, 0.0)
nest.SetStatus(
    step_current_generator,
    {"amplitude_times": list(steplist), "amplitude_values": amplist},
)


# Functions
def storage(  # Save data from simulations
    events, sources, targets, simulation_time, datapath, phase,
):
    times = events["times"]
    senders = events["senders"]
    extension = str(simulation_time) + "_rank_" + str(rank) + ".npy"
    datapath = datapath + phase + "/"
    np.save(datapath + "times_" + extension, times)
    np.save(datapath + "senders_" + extension, senders)
    np.save(datapath + "sources_" + extension, sources)
    np.save(datapath + "targets_" + extension, targets)


def storage_time(  # Save time points/steps
    stepsize, timesteps, index, datapath, phase,
):
    extension = ".npy"
    datapath = datapath + phase + "/"
    np.save(datapath + str(phase) + "_stepsize" + extension, stepsize)
    np.save(datapath + str(phase) + "_timesteps" + extension, timesteps)
    np.save(datapath + str(phase) + "_index" + extension, index)


def simulate_cicle(  # Simulation cycle and data recording
    time_steps, stimulation_time, save, datapath, phase, index
):
    step = np.mean(np.diff(time_steps))  # extract step size
    local_connections = []
    for simulation_time in time_steps:
        if np.isnan(step) == True:  # single_TMS
            step = stimulation_step

        nest.Simulate(step)

        local_connections = nest.GetConnections(pop_exc, pop_exc)
        sources = np.array(nest.GetStatus(local_connections, "source"))
        targets = np.array(nest.GetStatus(local_connections, "target"))

        events = nest.GetStatus(spike_detector, "events")[0]
        times = events["times"]
        senders = events["senders"]

        if save == 1:
            storage(
                events, sources, targets, simulation_time, datapath, phase,
            )
            storage_time(step, time_steps, index, datapath, phase)
        nest.SetStatus(spike_detector, "n_events", 0)  # flush spike_detector


def checkpoint(
    datapath, phase
):  # print checkpoint to indicate end of simulation of a phase in the phase data folder
    filename = datapath + phase + "/" + "MissionComplete.npy"
    x = 1  # arbitrary filler
    np.save(filename, x)


# Simulation
# Growth
phase = "growth"
index = 1
simulate_cicle(growth_steps, stimulation_time, 1, datapath, phase, index)
checkpoint(datapath, phase)

# Stimulation
phase = "stim"
index += cicles
simulate_cicle(stimulation_steps, stimulation_time, 1, datapath, phase, index)
checkpoint(datapath, phase)

# Post-stimulation
phase = "post"
index += len(stimulation_steps)
simulate_cicle(post_stimulation_steps, stimulation_time, 1, datapath, phase, index)
checkpoint(datapath, phase)

# Decay
phase = "decay"
index += len(post_stimulation_steps)
simulate_cicle(decay_steps, stimulation_time, 1, datapath, phase, index)
checkpoint(datapath, phase)
