These folders can be used to build paraview. It is not a fully functional superbuild though.
It is adviced to place the repro folder in:
$HOME/visualization/

And then run the command:
mv $HOME/visualization/repro/* $HOME/visualization/build/

Order
-----
1. $ bash config_llvm
2. $ bash config_mesa
3. $ bash config_paraview

Assumptions
-----------
-cmake is installed and available as command (If not install and make sure it is in included in $PATH)
-there is a python 2.7 installation in an conda environment named "py27"

Errors
------
-Paraview is still partly built against the default python installed in /usr/lib/python. Potential solution in Discourse by Ben Boeckel: 
"I’d just go through the CMakeCache.txt and find any mention of Python from /usr/lib and change it to the one in the Anaconda env."
I have not tested whether that works.

