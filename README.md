# Efficient determination of the transfer function via mechanical networks, for the calculation of the $H_2$ norm
## Gareth H. Willetts and Timothy H. Hughes, Department of Earth and Environmental Sciences, University of Exeter, Penryn, Cornwall.
## ghw205@exeter.ac.uk, t.h.hughes@exeter.ac.uk

Here you will find the associated code to reproduce the results of the paper titled "A graph theoretic approach to determining the transfer functions of mechanical networks, towards efficient computation of their H2 norms", submitted for CONTROL2024.

The files are as follows:

ReproducibleTrainScript.mlx

: is an annotated MATLAB livescript, which talks the user through the use of a graphical method to determine a transfer function of a model of a train suspension system, and using this transfer function to calculate the H2 norm (symbolically or numerically). This is used in conjunction with an interior point optimisation scheme to minimise the H2 norm.

inbuilt_optim_H2.m

: is a MATLAB function that uses inline functions to evaluate a symbolic transfer function numerically.

Inbuilt_optim.m

: uses a free-body diagram approach to calculate a symbolic expression of the transfer function for the model of a train suspension system. This is evaluated using "inbuilt_optim_H2.m", and then passed to an interior point optimisation scheme to optimise the H2 norm.

Train_Graph_Functions_Timing.m

: follows the same example as above, using the functions combine_edge.m, combine_vertices.m and sum_admittance_products.m to determine the transfer function only, and times how long MATLAB takes to execute the code, averaged over 1000 runs.
