# BOLD_noise
MATLAB toolbox for simulation of noise in BOLD data based on properties observed in real BOLD fMRI datasets

BOLD noise simulator
Aim of this toolbox is simulation of noise in BOLD data based on properties observed in real datasets. The simulated noise consists particularly from:
1. Noise caused by movements of subject's head
2. Noise caused by breathing and heartbeats
3. 1/f noise

Requirements
Toolbox was written in MATLAB 2015b with statistical toolbox. It uses several functions from toolbox SPM12.

Movement simulation
This simulation is performed in function get_movement_noise. First the head movements are simulated as translations and rotations in x, y and z axes . The initial head position is zero for all six of these movement parameters and differences are added in each step (scan) of the simulation. The added differences are based on normal distribution of differences observed in our 6 datasets (251 subjects). Next, the simulator computes differences, squares and squared differences of simulated movement parameters. It follows simulation of explained variability of BOLD signal by each of 24 movement parameters (EVMP). For this purpose, we estimated maps of EVMP using GLM model. The EVMP are voxelwise selected from observed exponential distribution. Simulation is performed inside the mask, where at least 80% of subjects from dataset had valid data. The movement noise is computed as a sum of multiplicated EVMPs by particular movement parameter. Simulated 4D movement noise is stored in the output variable ‘M_noise’.
 
Physiological noise simulation
Simulation is performed in function get_physio_noise. First, the RR intervals of cardiac activity and time varying breathing periods are simulated. The simulation is based on distribution of values observed in 2 datasets (52 subjects). We used RETROICOR to extract typical impulse response functions of artifacts related to breathing and cardiac activity. We also used this tool to get positions where physiological artifacts describe significant part of the measured signal. The simulation of the physiological noise is performed only for these positions. Impulse responses for cardiac and breathing artifacts consist from set of 8 basis functions. For simulation is selected one impulse response for cardiac and one for breathing artifact from observed distribution. In the next step, we are using Cholesky decomposition to keep correlation structure of the cardiac and breathing noise artifact. The correlation structure is selected from distribution of correlation coefficients between extracted physiological artifacts. Simulated 4D physiological noise is stored in the output variable ‘P_noise’.

1/f noise simulation
Simulation is performed in function get_1_f_noise. 1/f noise is characteristic with hyperbolic curve in amplitude spectra of BOLD fMRI data. The hyperbola is defined by parameter 1/f. We extracted BOLD signals from noise related positions in four of our datasets (170 subjects). Then we fitted hyperbola on the amplitude spectra and created distribution of 1/f parameters. For simulation is one 1/f parameter selected from normal distribution of observed parameters. In the next step is in the area of brain mask generated random noise whose amplitude spectra is multiplied by selected hyperbola. 1/f parameter for simulated 4D data has normal distribution with selected mean and standard deviation XX. Simulated 4D 1/f noise is stored in the output variable ‘F_noise’.

