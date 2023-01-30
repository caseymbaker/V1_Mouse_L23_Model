# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 22:32:46 2022

@author: cmbak
"""

import brian2
from brian2 import *
import random
import numpy as np
import matplotlib.pyplot as plt
import time
import scipy.io
import networkx as nx
import pandas as pd  
import itertools
from itertools import chain
import h5py
from scipy.io import savemat

defaultclock.dt = 0.1*ms

# %% Load Network Connections, Weights, and Background 
mat = scipy.io.loadmat('connections.mat') #adjacency matrix All
dat = mat.get("data")

matE = scipy.io.loadmat('ExcData.mat') #adjacency matrix Exc Presynaptic
datE = matE.get("dataE")
presynE = datE[:,0]-1 # -1 to convert indicies from matlab to python
presynE = presynE.astype(int)
postsynE = datE[:,1]-1
postsynE = postsynE.astype(int)
weightsE = datE[:,2] 

matI = scipy.io.loadmat('InhData.mat') #adjacency matrix Inh Presyn
datI = matI.get("dataI")
presynI = datI[:,0]-1 # -1 to convert indicies from matlab to python
presynI = presynI.astype(int)
postsynI = datI[:,1]-1
postsynI = postsynI.astype(int)
weightsI = datI[:,2] 

G = nx.DiGraph()
for k in range(len(dat)):
    G.add_edge(dat[k][0],dat[k][1],weight = dat[k][2])
    
neuron_num = max(dat[:,0:1])[0].astype(int)

# load background stim 
mat2 = scipy.io.loadmat('stimfile.mat')
stim = mat2.get("rstim4")
stim = np.transpose(stim)

# %%  run simulation 5 times (once for each neuron pair created in SelectNeuronsToStim.m)
for iz in range(1,6):

    node1stim = scipy.io.loadmat('nods'+str(iz)+'.mat') 
    nodes = node1stim.get("nodes")
    
    neuron_num = max(dat[:,0:1])[0].astype(int)
    
    
    volts = []
    mat2 = scipy.io.loadmat('stimfile.mat')
    stim = mat2.get("rstim4")
    stim = np.transpose(stim)
    gl = 10.0*nsiemens   # Leak conductance
    erest = np.full((1, 3200), -50)
    irest1 = np.full((1,330), -55)
    irest2 = np.full((1,330), -65)
    irest3 = np.full((1,140), -65)
    elres = np.hstack((erest,irest1,irest2,irest3)) #leak voltage
    erest2 = np.full((1, 3200), -68)
    neighbor3 = []
    irest12 = np.full((1,330), -65)
    irest22= np.full((1,330), -60)
    irest32 = np.full((1,140), -65)
    vres = np.hstack((erest2,irest12,irest22,irest32)) #initial voltage
    # reset potential
    res = ''' 
        v = -70*mV 
        g_ampa = 0*nS
        g_gaba = 0*nS
     '''         
    er = -80*mV          # Inhibitory reversal potential
    vt = -55.0*mV         # Spiking threshold
    memc = 200.0*pfarad  # Membrane capacitance
    tau_ampa = 5*ms   # Glutamatergic synaptic time constant
    tau_gaba = 10.0*ms  # GABAergic synaptic time constant
    spiketrain =np.ones((1,20000))*100 #remove this at the end, just initializing 

    
    for j in range(5):
            
        
        #stimulate neurons for 6 frames (0.6 ms) every 333ms (30Hz) for 250 ms (starting after 1s)
        stimlist=numpy.zeros((20000,neuron_num))
        pulselist=numpy.arange(10000,12500,333)
        
        stimlist[pulselist,(nodes-1)]=150 #unit = siemens 
        stimlist[pulselist+1,(nodes-1)]=150#unit = siemens 
        stimlist[pulselist+2,(nodes-1)]=150#unit = siemens 
        stimlist[pulselist+3,(nodes-1)]=150#unit = siemens 
        stimlist[pulselist+4,(nodes-1)]=150#unit = siemens 
        stimlist[pulselist+5,(nodes-1)]=150#unit = siemens 
    
        P=TimedArray(stimlist*nS,dt=defaultclock.dt)
        P2 = TimedArray(stim*pA,dt =defaultclock.dt)
    
       
        eqs = '''
       
        dv/dt=(-gl*(v-el)-(g_ampa*v+g_gaba*(v-er))-P(t,i)*v+P2(t,i))/memc : volt 
        dg_ampa/dt = -g_ampa/tau_ampa : siemens (unless refractory)
        dg_gaba/dt = -g_gaba/tau_gaba : siemens (unless refractory)
        el : volt
    
        '''
       
        group = NeuronGroup(neuron_num, eqs, threshold = 'v > vt', refractory = 2*ms, reset = res ,method = 'euler')
        
        group.el = elres*mV
        group.v = vres*mV       

        S = Synapses(group,group,model='''w : 1''', on_pre='g_ampa += (w)*0.3*nS',delay = 1.5*ms) 
        S.connect(i = presynE, j = postsynE)
        S.w = weightsE
        
        
        dlys = list(np.random.randint(low = 4, high = 6, size = len(postsynI)))
        S2 = Synapses(group,group,model='''w : 1''', on_pre='g_gaba += (w)*3*nS') 
        S2.connect(i = presynI, j = postsynI)
        S2.w = weightsI
        S2.delay = dlys*ms
        
        p1 = PoissonInput(group[0:3199],'g_ampa', N=1, rate = .08*Hz,weight = 30*nS)
        p4 = PoissonInput(group[3530:3859],'g_ampa', N=1, rate = .05*Hz,weight = 30*nS)
        p2 = PoissonInput(group[3200:3529],'g_ampa', N=1, rate = 1.6*Hz,weight = 30*nS)
        p3 = PoissonInput(group[3860:3999],'g_ampa', N=1, rate = 1.35*Hz,weight = 30*nS)
    
    
    
        M = StateMonitor(group,('v','g_ampa','g_gaba'),record=True)
        M2 = SpikeMonitor(group)
        run(2*second)
        plot(M2.t/ms,M2.i,'.')
 
    
        t1 = numpy.zeros((neuron_num,20000))
        trains = M2.spike_trains()
        for z in range(neuron_num):
            spikes = trains[z]
            for k in range(len(spikes)):
                q = spikes[k]*ms
                q = round(q*10/ms/ms);
                t1[z,q]=1
                
        spiketrain = np.vstack((spiketrain,t1)) 
        
        volts.append(M.v)
    
    
       
    spiketrain = spiketrain[1:len(spiketrain),]
    savemat("spikes" +str(iz)+".mat", mdict={'arr': spiketrain})
    savemat("volt" + str(iz)+".mat", mdict={'volts': volts})





