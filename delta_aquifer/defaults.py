# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 15:54:29 2019

@author: engelen
"""

from collections import OrderedDict

def get_seawat_default_runfile(): 
    seawat_default_runfile = OrderedDict(
        [
            # gen
            ("modelname", "results"),
            ("writehelp", True),
            # dis
            ("nstp", 1),
            ("sstr", "tr"),
            ("laycbd", 0),
            # bas6
            ("hnoflo", -9999.0),
            # oc
            ("savehead", True),
            ("saveconclayer", True),
            ("savebudget", False),
            ("saveheadtec", False),
            ("saveconctec", False),
            ("savevxtec", False),
            ("savevytec", False),
            ("savevztec", False),
            # lfp
            ("ilpfcb", 1),
            ("hdry", 1.0e30),
            ("nplpf", 0),
            ("laytyp", 0),
            ("layavg", 0),
            ("chani", 1.0),
            ("layvka", 0),
            ("sto", 1e-5),
            # pcg
            ("mxiter", 100),
            ("iter1", 30),
            ("hclose", 0.0001),
            ("rclose", 1.0),
            ("relax", 0.98),
            ("nbpol", 0),
            ("iprpcg", 1),
            ("mutpcg", 1),
            # pksf
            ("pksf", False),
            ("mxiterpks", 1000),
            ("inneritpks", 30),
            ("hclosepks", 0.0001),
            ("rclosepks", 1.0),
            ("npc", 2),
            ("partopt", 0),
            ("pressakey", False),
            # btn
            ("cinact", -9999.0),
            ("thkmin", 0.01),
            ("nprs", 0),
            ("ifmtcn", -1),
            ("chkmas", True),
            ("nprmas", 10),
            ("nprobs", 1),
            ("tsmult", 1.0),
            ("dt0", 0.0),
            ("mxstrn", 10000.0),
            ("ttsmult", 1.0),
            ("ttsmax", 0.0),
            ("por", 0.35),
            # adv
            ("mixelm", -1),
            ("percel", 0.8),
            # dsp
            ("al", 1),
            ("trpt", 0.1),
            ("trpv", 0.01),
            ("dmcoef", 1.e-9),
            # gcg
            ("mt3d_mxiter", 1000),
            ("mt3d_iter1", 300),
            ("mt3d_isolve", 2),
            # vdf
            ("mtdnconc", 1),
            ("mfnadvfd", 2),
            ("nswtcpl", 1),
            ("iwtable", 0),
            ("densemin", 1000.0),
            ("densemax", 1025.0),
            ("denseref", 1000.0),
            ("denseslp", 0.7143),
            # drn
            ("mxactd", 1.0e6),
            ("idrncb", 0),
            # chd
            ("mxactc", 1.0e6),
            # gbh
            ("mxactb", 1.0e6),
            ("ighbcb", 0),
            # riv
            ("mxactr", 1.0e6),
            ("irivcb", 0),
            # rch
            ("nrchop", 3),
            ("irchcb", 0),
            # wel
            ("mxactw", 1.0e6),
            ("iwelcb", 0),
            # ssm
            ("mxss", 1.0e6),
        ]
    )
    return(seawat_default_runfile)