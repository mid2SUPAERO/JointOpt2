# JointOpt2

In this depository we have the V2 version with new mathematical constraints (simplified):

The simulation code  FGASLJ1DbarTEPS.m

The raw optimization main.m (using finite differences and nl constraints)

the MultiStart optimization main_ms.m (need global optimization toolbox + Parrallel Computing toolbox) *

Need to: 

add a new objective function (such as min variance)
treat constraints as linear : Ax=B


*
https://fr.mathworks.com/help/gads/example-parallel-multistart.html
