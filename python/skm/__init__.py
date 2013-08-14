"""
Module for A_FB Skim

"""

__author__ = "Xin Shi <Xin.Shi@cern.ch>"
__copyright__ = "Copyright (c) Xin Shi"

import sys 

def main(args):
    if args[0] == 'mc':
        return skim_mc(args[1:])

    module_name = 'skm.%s' % args[0]
    top_module = __import__(module_name)
    module = sys.modules[module_name]
    module.main(args[1:])


def skim_mc(args):
    label = args[0]

    if label == 'B2KstarMuMu/RECO_100M_v1.6':
        cmd = 'This is a dummy command'
    print cmd 
    sys.exit()

    



    
