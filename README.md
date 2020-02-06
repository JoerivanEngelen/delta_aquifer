# delta_aquifer

Scripts to create, process and post-process a 3D variable-density groundwater model of a synthetic delta. 
The synthetic delta can be created based on a set of parameters that allow for a wide range of geometries, 
hydrogeological parameterizations and boundary conditions. An explanation of the concepts and the literature study 
these are based on can be found in the accompanying paper that is currently in the works.

As of currently used to conduct a global sensitivity analysis on the Dutch national cluster (Cartesius):


The repository is structured as follows:
* _analyze_input_.  
Scripts to analyze the input distributions are here
* _data_.  
Data is located here.
* _debug_help_.  
Dummy models are located here, to allow for quicker debugging of this package. Poorly supported
* _delta_aquifer_.  
The main workhorse, where scripts to construct models are located.
* _example_.  
Examples where parts of a model are created and plotted. Also contains jobscripts that can be used on clusters. Compiling instructions are also located here.
* _post_.  
Post-processing scripts and Paraview macros.
* _process_.  
Scripts to convert output of one model to input for another model.  
Required for models that run longer than 8000 years, due to a calendar flaw in iMOD-SEAWAT.
* _test_.  
Should contain tests. No support now.
