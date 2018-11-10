Author: Cathy Zhang
Date: 2018.11.04

This package is an implementation of deblurring and denoising via Total Variation method.


----------------- Code Descripition ----------------

Run 'main.m' to test the rest the result. You can change the parameters in this file to get desirable results.

NOTE: Besides the required inputs for the main routine, you can also feed in a struct 'opts' into the the function to specify other optional parameters such as mu & tmax & type. See the comments in the 'main.m' function or use 'help main' to view details.

'tv_deblur_cyclic.m' & 'tv_deblur_noncyclic.m' are two major functions, which are the implementation of tv method with periodic boundary condition blur and non-periodic boundary condition blur respectively.

See the comments in the function files for more detail. You can also use the 'help' function in MATLAB to view the detailed description.


--------- Recommended Choice for HyperParameters -------

lambda = 1e-3
opts.mu = 1
opts.tmax = 30

-------------- Image file Specification -------------

All the testing image file is contained in the folder '\image'. There are four testing images: 'classic.jpg', 'kamiya.jpg', 'koala.jpg', 'shizuwo.jpg'.

The blurred and recovered image is contained in the folder '\image\recovered image cyclic' & '\image\recovered image noncyclic'. Different folders correspond to different boundary conditions for the blur.
 