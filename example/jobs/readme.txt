This shows an example script of how to run a job that creates, runs and visualizes a model.
You will need to adapt the paths to your own system. These examples are purely intended as a 
reference of what you ideally can do with this package. You will probably have to spend some
time installing all the dependencies though :(

Start the jobscript by using the command:

$ sbatch Run_Nr *nr*

With *nr* being a placeholder for an integer in the range of 0 to ((nr_parameters+1) * nr_trajectories)

DEPENDENCIES
------------

-SLURM as workload manager
 ->If your cluster is using SGE, this page might be helpful to convert https://srcc.stanford.edu/sge-slurm-conversion

-Python 3.6 environment (with xarray and imod-python installed)

-GNU Parallel

-iMOD-SEAWAT (ask permission to join the SVN repository to compile it)

-Paraview with OSMESA support (see delta_aquifer/examples/configure with scripts to install this)

