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


def ntp_labels(label):
    if '/' in label: 
        return [label]

    if 'Run2011v10.' in label: 
        subversion = label.replace('Run2011v10.', '')
        ntp_labels = [
            'Run2011v10.%s/May10ReReco_v1.%s'       %(subversion, subversion), 
            'Run2011v10.%s/PromptReco_v4_10.%s'     %(subversion, subversion), 
            'Run2011v10.%s/PromptReco_v5_10.%s'     %(subversion, subversion), 
            'Run2011v10.%s/PromptReco_v6_10.%s'     %(subversion, subversion), 
            'Run2011v10.%s/B_PromptReco_v1_10_1.%s' %(subversion, subversion), 
        ]

    elif 'Run2011v11.' in label: 
        mainversion = '11'
        subversion = label.replace('Run2011v11.', '')
        ntp_labels = [
            '%s/May10ReReco_v1_%s.%s'    %(label, mainversion, subversion), 
            '%s/PromptReco_v4_%s.%s'     %(label, mainversion, subversion), 
            '%s/PromptReco_v5_%s.%s'     %(label, mainversion, subversion), 
            '%s/PromptReco_v6_%s.%s'     %(label, mainversion, subversion), 
            #'Run2011v10.%s/B_PromptReco_v1_10_1.%s' %(subversion, subversion), 
        ]

    else:
        raise NameError(label)

    return ntp_labels


def njobs(label):
    if 'Run2011v10.' in label and '/May10ReReco_v1.' in label:
        njobs = 3

    elif 'Run2011v10.' in label and '/PromptReco_v4_10.' in label: 
        njobs = 10

    elif 'Run2011v10.' in label and '/PromptReco_v5_10.' in label: 
        njobs = 4

    elif 'Run2011v10.' in label and '/PromptReco_v6_10.' in label: 
        njobs = 5

    elif 'Run2011v10.' in label and '/B_PromptReco_v1_10_1.' in label: 
        njobs = 33
    
    if 'Run2011v11.' in label and '/May10ReReco_v1_11.' in label:
        njobs = 20

    elif 'Run2011v11.' in label and '/PromptReco_v4_11.' in label: 
        njobs = 68

    elif 'Run2011v11.' in label and '/PromptReco_v5_11.' in label: 
        njobs = 25

    elif 'Run2011v11.' in label and '/PromptReco_v6_11.' in label: 
        njobs = 35

    elif 'Run2011v11.' in label and '/B_PromptReco_v1_11_1.' in label: 
        njobs = 0 
    
    else:
        raise NameError(label)

    return njobs

def procdir(label):
    if 'Run2011v11' in label: 
        procdir = os.path.join(atr.afbpath, 'rel', 'CMSSW_4_2_8_patch7', 
                               'src/BphAna/BToKstarMuMu/plugins') 
    elif 'BuToKstarMuMu/7TeV_2E7_v1' in label: 
        procdir = os.path.join(atr.afbpath, 'rel', 'CMSSW_4_2_8_patch7', 
                               'src/BphAna/BToKstarMuMu_v1/plugins') 
    else: 
        raise NameError(label)

    return procdir
