"""
Module for A_FB Selection

"""

__author__ = "Xin Shi <Xin.Shi@cern.ch>"
__copyright__ = "Copyright (c) Xin Shi"

import sys 

def main(args):
    module_name = 'sel.%s' % args[0]
    top_module = __import__(module_name)
    module = sys.modules[module_name]
    module.main(args[1:])

