"""
Attributes for the NTuples.

"""

__author__ = "Xin Shi <Xin.Shi@cern.ch>"
__copyright__ = "Copyright (c) Xin Shi"


import os
import sys
from tls import *
import atr 
import ntp 
# from tls import * # !! This is not going to work, since the main program
# used atr.ntppath!! So, use import inside the function through the atr module.
# Now, this has changed. The idea is to use tls everywhere, make the tls generic! 


def ntppath(label, datatype):
    from atr import castor_datpath, datpath

    if label == '5ifbv1':
        path_ = os.path.join(datpath, 'ntp')

    elif label in ['B2KstarMuMu/GEN_1M_v1.1', 'B2KstarMuMu/RECO_1M_v2.2',
                   'B2KstarMuMu/RECO_1M_v2.2/GEN', 'B2KstarMuMu/RECO_1M_v2.3']:
        path_ = os.path.join(datpath, 'ntp', datatype, label)

    else:
        raise NameError(label)

    return path_


def rootfile(datatype, label, test=False):
    com_name = get_name_from_label(label)
    function = getattr(ntp, com_name)
    return function(datatype, label)


def rootpath(datatype, label, test=False):
    #eosbase = 'root://eoscms//eos/cms/store/user/xshi/'
    datpath = atr.datpath 

    if datatype == 'data': 
        ntp_label = get_label_by_version(label, 'x')
        inpath = os.path.join(datpath, 'ntp/data', ntp_label)
    else:
        inpath = os.path.join(datpath, 'ntp/mc', label)

    return inpath 


def rootname(datatype, label, batch=False):
    if datatype == 'mc': 
        name = 'BToKstarMuMu_merged_1'
    else:
         name = 'BToKstarMuMu_merged_*'
        
    if batch :
        name = 'BToKstarMuMu_merged_${LSB_JOBINDEX}'
        
    return name


def datasets(datatype, label, test=False):
    com_name = get_name_from_label(label)
    function = getattr(ntp, '%s' %com_name)
    return function(datatype, label)


def root_tree_obj(datatype, label):
    if datatype == 'mc':
        if  'B2KstarMuMu/GEN_1M_v1.1' in label or \
            'B2KstarMuMu/RECO_100M_v1.6' in label:
            obj = 'B0Cand/B0KstMuMuNTuple'

        elif 'B2KstarMuMu/RECO_1M_v2.2' in label or \
        'B2KstarMuMu/RECO_1M_v2.2' in label or \
        'B2KstarMuMu/RECO_1M_v2.3' in label or \
        'B2KstarMuMu/RECO_100M_v1.5' in label: 
             obj = 'data'
        else:
            raise NameError(label)

    else:
        raise NameError(label)
    return obj


def B2KstarMuMu_RECO_Brian_prod_2_1_v3(datatype, label, com_name):
    eosbase = 'root://eoscms//eos/cms/store/user/xshi/'
    primarydataset = 'B0dKstMuMuKPi_MuMuEtaFilter_pythia6_evtgen_prod2_1'
    psethash = 'b5872eb8acdaef414564c8748b88ada2'
    #filepath = os.path.join('/store/user/xshi', primarydataset, comname, psethash)
    rootfile = os.path.join(eosbase, primarydataset, com_name, psethash, '*.root')
    return rootfile
        
def Run2011v1(datatype, label, com_name):
    rootfile = os.path.join(atr.afbpath, 'dat/ntp/data/Run2011v1', '*.root')
    return rootfile


def Run2011v10(datatype, label):
    eosbase = 'root://eoscms//eos/cms/store/user/xshi/'
    subs = [('May10ReReco_v1', 9),
            ('PromptReco_v4_10', 10), 
            ('PromptReco_v5_10', 4),
            ('PromptReco_v6_10', 5),
            ('B_PromptReco_v1_10_1', 33)]
    p = [ ('%s/Run2011v10/%s' % (eosbase, sub), njobs) for sub, njobs in subs ]
    return p 

    

def datasets_njobs(label):
    base_label = get_label_by_version(label, 'x')
    if base_label == 'Run2011v10':
        datasets_njobs = [('May10ReReco_v1', 3),
                          ('PromptReco_v4_10', 10), 
                          ('PromptReco_v5_10', 4),
                          ('PromptReco_v6_10', 5),
                          ('B_PromptReco_v1_10_1', 33)]

    elif base_label == 'Run2011v10/May10ReReco_v1':
        datasets_njobs = [('May10ReReco_v1', 3)]

    elif base_label == 'Run2011v10/PromptReco_v4_10':
        datasets_njobs = [('PromptReco_v4_10', 10)]

    elif base_label == 'Run2011v10/PromptReco_v5_10':
        datasets_njobs = [('PromptReco_v5_10', 4)]

    elif base_label == 'Run2011v10/PromptReco_v6_10':
        datasets_njobs = [('PromptReco_v6_10', 5)]
    
    elif base_label == 'Run2011v10/B_PromptReco_v1_10_1':
        datasets_njobs = [('B_PromptReco_v1_10_1', 33)]

    else:
        raise NameError(label)

    return datasets_njobs

        
def grid_path(label):
    com_name = get_name_from_label(label)
    eosbase = '/afs/cern.ch/user/x/xshi/eos/cms/store/user/xshi/'
    if label in ['Run2011A-PromptReco-v4_run2011v1', 
                 'Run2011A-PromptReco-v5_run2011v1', 
                 'Run2011A-PromptReco-v6_run2011v1', 
                 'Run2011B-PromptReco-v1_run2011v1'
                 ] : 
        srcdir = os.path.join(
            eosbase, 'MuOnia',
            com_name, '12d093a44a0e90cf594f6f76582ad92e')
    elif label in ['Run2011A-May10ReReco-v1_run2011v1'] : 
        srcdir = os.path.join(
            eosbase, 'MuOnia',
            com_name, '832650ad13378e362824b793e74b5ef')

    elif label in ['BuToKstarJPsi-7TeV-5E5-v1_run2011v1', 
                   ] : 
        srcdir = os.path.join(
            eosbase, 'BuToKstarJPsi_EtaPtFilter_7TeV-pythia6-evtgen',
            com_name, 'a2cbd59f90c5e3694010a7612c1b81d6')

    elif 'Run2012' in label:
        srcdir = os.path.join(eosbase, com_name)

    else:
        raise NameError(label)
    return srcdir 
    
    
num_rootfiles = {
    'Run2011A_May10ReReco_v1_run2011v0_2' : 2, 
    'Run2011A_PromptReco_v4_run2011v0_1' : 3, 
    'Run2011A_PromptReco_v5_run2011v0': 1, 
    'Run2011A_PromptReco_v6_run2011v0': 1, 
    'Run2011B_PromptReco_v1_run2011v1': 1, 
    'BuToKstarMuMu_7TeV_2E7_v1_run2011v1_1': 1, 
    }

