# SL-SP-MRF
A simple implementation of the SL-SP algorithm for MRF reconstruction:

**”High-quality MR Fingerprinting Reconstruction using Structured Low-rank Matrix Completion and Subspace Projection“**
# Start
To run the simulation open: **"main.m"** in Matlab and press Run button.
# Code Ref
Some of the codes are taken from
1. Paper：“Low-rank magnetic resonance fingerprinting”
http://webee.techni-on.ac.il/Sites/People/YoninaEldar/software_det18.php
2. Paper：“A Fast Algorithm for Convolutional Structured Low-Rank Matrix Recovery”
https://github.com/cbig-iowa/giraf
3. Brian Hargreaves webpage: 
http://mrsrl.stanford.edu/~brian/mritools.html
4. Michael Lustig webpage:
https://people.eecs.berkeley.edu/~mlustig/

# Code details
* **Bloch**: FISP simulations
* **data**: store data for simulation, including FA, TR, Mask, etc..
* **methods**: SL-SP methods
* **simulations**: dictionary generation, dictionary matching, etc..
*  showRseults: plot results

The code is tested in the matlab2021 on a Windows workstation with an Intel Xeon 3.80 GHz
CPU and an Nvidia Quadro GV100 GPU.
