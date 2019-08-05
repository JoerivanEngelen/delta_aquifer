Compiling iMOD seawat

Location SVN branch:
https://repos.deltares.nl/repos/seawat/branches/imod-wq_pks/

Dependencies
------------

-Cmake

-Python 2.7

-mpi (intel flavour)

-fortran (intel flavour)

Example
-------
$ cd /home/jengelen/seawat/compile/
$ module load python
$ python makegen.py
$ cd src
$ module load mpi/impi
$ module load fortran/intel
$ make


