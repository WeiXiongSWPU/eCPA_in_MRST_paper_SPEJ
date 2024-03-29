# Code supplement for  "Compositional simulation for carbon storage in porous media using an electrolyte association equation of state"

This repo contains the code for the compostional simulator developed for the
paper "Compositional simulation for carbon storage in porous media using an
electrolyte association equation of state" (submitted to the SPE Journal). This
repository contains a stripped-down version of MRST 2022b that contains the
necessary code to run the examples. The code is planned for integration in an
upcoming MRST release, but this repository will remain as a self-contained
snapshot of the code used to produce the figures and results in the paper.

Before you run the examples, the startup.m file must be run to activate MRST.
Once this is done, you can find the example scripts in the folder code_for_spej.
To run the first example you can run the following in Matlab:

```matlab
cd <path/to/repository>
run startup.m
cd code_for_spej
run case1_single_cell.m
```

Please refer to [MRST website](https://www.sintef.no/projectweb/mrst) on how to run
and explore the examples.
