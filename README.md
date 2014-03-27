harlemcmc-shake
===============

This code generates two short animations illustrate the differences between a Metropolis-Hastings (MH) sampler and a Hamiltonian Monte Carlo (HMC) sampler. 

The animations can be seen at 
https://www.youtube.com/watch?v=Vv3f0QNWvWQ

We test the samplers on nine mixtures of Gaussians, ranging from a mixture of four elongated Gaussians in the upper right-hand corner, to a single spherical Gaussian in the middle, to a mixture of 100 spherical Gaussians in the lower lefthand corner. We plot the contours of these distributions in white and the last 10 samples in red.

In the MH case, the samples are connected in sequential order by yellow lines. In the HMC case, the samples are connected (in yellow) by the Hamiltonian evolution sequence generated during the proposal of the next (accepted) sample. The MH proposals are set to spherical Gaussians with relatively small variance in the outer eight target Gaussian mixtures and relatively large variance for the central spherical Gaussian target distribution.

Authors:
Tamara Broderick (UC Berkeley)
http://www.stat.berkeley.edu/~tab/

David Duvenaud (University of Cambridge)
http://mlg.eng.cam.ac.uk/duvenaud/


Enjoy!

