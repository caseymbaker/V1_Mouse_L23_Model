# V1_Mouse_L23_Model
Code to create a spiking model of layer 2/3 of mouse visual cortex and simulate optogenetic perturbation of pattern completion neurons in ensembles

If you use any code from this repository, please cite: 

## Background/Motivation
***Cortical Ensembles*** are groups of densely connected neurons that fire together to encode stimuli. [Recent Work](https://www.cell.com/cell/pdf/S0092-8674(19)30616-6.pdf) found preliminary evidence that ensemble activity in layer 2/3 of the mouse visual cortex can drive behavior. However, the link between ensembles and behavior is underdeveloped. Establishing a causal link between ensembles and behavior requires an efficient way to activate ensembles in the absence of stimuli. 

 ***Pattern Completion Neurons*** are pairs of neurons that, when stimulated, can activate the entire ensemble. However, selection of pattern completion neurons *in vivo* is challenging. We used our model to identify properties of efficient pattern completion neurons. 
 
## Network Information 

**Highlights:**
- 4000 Leaky Integrate-and-Fire Neurons
    - 3200 Excitatory
    - 800 Inhibitory
       - 330 PV
       - 330 SOM
       - 140 VIP
- Connection weights follow a lognormal distribution 
- Excitatory connections are influenced by similarity in preferred orientation

**Connectivity Pattern:**

![plot](./Miscellaneous/schematic.png)

*see [Methods](publication link) for more details*

## Code Information

1. Run CreateNetwork/MakeNetwork.m to generate the connection weights between all neurons and generate background current 
 - change numframs and bc variables based on if you are going to stimulate neurons in the network
2. Run CreateNetwork/SelectNeuronsToStim.m to generate pairs of neurons in the same ensemble to stimulate
3. Run code in RunStimulation folder 
 - output: volts.mat and spikes.mat which are variables that have the voltage and spiking information for each neuron for each trial
