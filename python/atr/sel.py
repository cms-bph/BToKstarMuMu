"""
Attributes for the Selected Files. 

"""

__author__ = "Xin Shi <Xin.Shi@cern.ch>"
__copyright__ = "Copyright (c) Xin Shi"


import os
import sys
from tls import *
import atr 
import sel 

# Maintain one fuction. 
# def set_rootfile(datatype, label, name, test=False):
#     if test:
#         name += '_test'
#     name_ = name + '.root'

#     path_ = os.path.join(atr.datpath, 'sel', datatype, label)
#     file_ = check_and_join(path_, name_)
#     return file_

def rootfile(datatype, label, test=False):
    com_name = get_name_from_label(label, ver='x')
    function = getattr(sel, com_name)
    return function(datatype, label)


# def rootfile(datatype, label, name,
#              remote=False, create=False, test=False):
#     if test:
#         name += '_test'

#     name_ = name + '.root'
#     path_ = os.path.join(atr.datpath, 'sel', datatype, label)

#     if remote:
#         path_ = os.path.join(atr.rafbpath, 'dat', 'sel', datatype, label)

#     f = os.path.join(path_, name_)

#     if create: 
#         f = check_and_join(path_, name_)

#     if not create and not remote and not os.path.isfile(f):
#         raise IOError(f)
    
#     return f
    
def Run2011v10(datatype, label):
    filepath = atr.datpath 
    f = os.path.join(atr.datpath, 'sel', datatype, label, '*.root')
    return f 


def Run2011v11(datatype, label):
    return Run2011v10(datatype, label)


def root_tree_obj(datatype, label):
    obj = 'tree'
    return obj


def BuToKstarMuMu_7TeV_2p5E3_v1_10(datatype, label):
    return Run2011v10(datatype, label)


def ntp_labels(datatype, label):
    if datatype == 'mc' and 'BuToKstarJPsi_7TeV_5E5_v1_run2011v0' in label: 
        ntp_labels = [ label ] 

    elif datatype == 'data' and 'run2011v0' in label: 
        ntp_labels = [
            'Run2011A_May10ReReco_v1_run2011v0_2',
            'Run2011A_PromptReco_v4_run2011v0_1', 
            'Run2011A_PromptReco_v5_run2011v0', 
            'Run2011A_PromptReco_v6_run2011v0', 
            'Run2011B_PromptReco_v1_run2011v0'
        ]

    else:
        raise NameError(label)

    return ntp_labels


def procdir(label):
    if 'run2011v0' in label: 
        procdir = os.path.join(atr.afbpath, 'rel', 'CMSSW_4_2_8_patch7', 
                               'src/BphAna/BToKstarMuMu_run2011v0/plugins') 
    else: 
        raise NameError(label)

    return procdir
