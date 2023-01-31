import brian2
from brian2 import *
import numpy as np

gl = 10.0*nsiemens   # Leak conductance
el = -70*mV          # Resting potential

er = -80*mV          # Inhibitory reversal potential
vt = -55.0*mV         # Spiking threshold
memc = 200.0*pfarad  # Membrane capacitance
tau_ampa = 5*ms   # Glutamatergic synaptic time constant
tau_gaba = 10.0*ms  # GABAergic synaptic time constant
res = '''
    v = -80*mV 
    g_ampa = 0*nS
    g_gaba = 0*nS
 '''         # reset potential

c = []
d = np.arange(0, 45, 0.5).tolist()

# %% run if testing exc connections
for k in d:
    eqs = '''
       
       dv/dt=(-gl*(v-el)-(g_ampa*v+g_gaba*(v-er)))/memc : volt 
       dg_ampa/dt = -g_ampa/tau_ampa : siemens (unless refractory)
       dg_gaba/dt = -g_gaba/tau_gaba : siemens (unless refractory)
    
        '''
    
    G = NeuronGroup(2, eqs, threshold = 'v > vt', refractory = 2*ms, reset = res ,method = 'euler')
    S = Synapses(G,G,model='''w : 1''', on_pre='g_ampa += (w)*0.3*nS') 
    S.connect(i=0, j=1)
    S.w = k
    
    M = StateMonitor(G, 'v', record=True)
    G.v = -70*mV
    run(200*ms)
    G.v[0] = -54.9*mV
    run(200*ms)
    plot(M.v[1])
    a = max(M.v[1])
    b = float(a)*1000+70
    c.append(b)

# %% run if testing inh connections
c = []
d = np.arange(0, 4, 0.5).tolist()

for k in d:
    eqs = '''
       
       dv/dt=(-gl*(v-el)-(g_ampa*v+g_gaba*(v-er)))/memc : volt 
       dg_ampa/dt = -g_ampa/tau_ampa : siemens (unless refractory)
       dg_gaba/dt = -g_gaba/tau_gaba : siemens (unless refractory)
    
        '''
    
    G = NeuronGroup(2, eqs, threshold = 'v > vt', refractory = 2*ms, reset = res ,method = 'euler')
    S = Synapses(G,G,model='''w : 1''', on_pre='g_gaba += (w)*3*nS') #change of 0.3 mV
    S.connect(i=0, j=1)
    S.w = k
    
    M = StateMonitor(G, 'v', record=True)
    G.v = -70*mV
    run(200*ms)
    G.v[0] = -54.9*mV
    run(200*ms)
    plot(M.v[1])
    a = min(M.v[1])
    b = float(a)*1000+70 
    c.append(b)
    
    

