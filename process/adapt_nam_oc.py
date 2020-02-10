# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 11:53:11 2020

@author: engelen

Adapt namfile and ocfile to write 

"""

import sys, os
from glob import glob

#%%Path management
if len(sys.argv) > 1:
    modelfol  = sys.argv[1]
else:
    ##For Testing
    modelfol = r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\iMOD-SEAWAT\modifications"

path_nam = glob(os.path.join(modelfol, "*.nam"))[0]
path_oc  = glob(os.path.join(modelfol, "*.oc"))[0]

#%%Read namfile
with open(path_nam, 'r') as f:
    namlines = f.readlines()

namlines = [l for l in namlines if l != ""] #Filter optional last empty lines out

lastid = int(namlines[-1].split(" ")[1])

namlines.append("DATA(BINARYPAR) {} heads.hds REPLACE\n".format(lastid+1))
namlines.append("DATA(BINARYPAR) {} budget.cbb REPLACE\n".format(lastid+2))

with open(path_nam, "w") as f:
    f.writelines(namlines)

#%%Read output control
with open(path_oc, 'r') as f:
    oclines = f.readlines()
    
oclines = [l for l in oclines if l != ""] #Filter optional last empty lines out

headline = oclines.pop(1)
headline = headline.split(" ")
headline[-1] = str(lastid + 1) + "\n"
oclines.insert(1, " ".join(headline))

oclines.insert(2, "COMPACT BUDGET\n")

with open(path_oc, "w") as f:
    f.writelines(oclines)