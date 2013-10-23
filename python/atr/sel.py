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

def rootfile(datatype, label, test=False):
    com_name = get_name_from_label(label, ver='x')
    function = getattr(sel, com_name)
    return function(datatype, label)

def Run2011v1(datatype, label):
    filepath = atr.datpath 
    f = os.path.join(atr.datpath, 'sel', datatype, label, '*.root')
    return f 


def root_tree_obj(datatype, label):
    obj = 'tree'
    return obj


def BuToKstarMuMu_7TeV_2p5E3_v1_10(datatype, label):
    return Run2011v10(datatype, label)


def ntp_labels(datatype, label):
    if datatype == 'mc' and 'BuToKstarJPsi_7TeV_5E5_v1_run2011v0' in label: 
        ntp_labels = [ label ] 

    elif datatype == 'mc' and 'BuToKstarMuMu_7TeV_2E7_v1_run2011v0' in label: 
        ntp_labels = [ label ] 

    elif datatype == 'data' and 'run2011v1' in label: 
        ntp_labels = [
            'Run2011A_May10ReReco_v1_run2011v1',
            'Run2011A_PromptReco_v4_run2011v1', 
            'Run2011A_PromptReco_v5_run2011v1', 
            'Run2011A_PromptReco_v6_run2011v1', 
            'Run2011B_PromptReco_v1_run2011v1'
        ]

    else:
        raise NameError(label)

    return ntp_labels


def procdir(label):
    if 'run2011v1' in label: 
        procdir = os.path.join(atr.afbpath, 'rel', 'CMSSW_4_2_8_patch7', 
                               'src/BphAna/BToKstarMuMu_run2011v1/plugins') 
    else: 
        raise NameError(label)

    return procdir
