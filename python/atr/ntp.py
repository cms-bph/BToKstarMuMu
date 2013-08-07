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

    elif label in ['5ifbv2', '5ifbv2.2']:
        path_ = ['/castor/cern.ch/user/z/zhlinl/b2KstarMuMu/TTree/10Dec2011/May10ReReco-v1',
                 '/castor/cern.ch/user/z/zhlinl/b2KstarMuMu/TTree/10Dec2011/Run2011A-PromptReco-v4',
                 '/castor/cern.ch/user/z/zhlinl/b2KstarMuMu/TTree/10Dec2011/Run2011A-PromptReco-v5'
                 ]
    
        
    elif label == '5ifbv2.1':
        path_ = os.path.join(datpath, 'skm')


    elif label in ['5ifbv2.3', '5ifbv2.3.1', '5ifbv2.3.2',
                   '5ifbv2.6.1', '5ifbv2.6.2']:
        path_label = get_label_by_version(label, 'x.x')
        ds = [ 'Run2011A-PromptReco-v5', 'Run2011A-PromptReco-v4',
               'Run2011A-May10ReReco-v1', 'Run2011A-PromptReco-v6',
               'Run2011B-PromptReco-v1']

        path_ = []
        for dataset in ds:
            path_.append('%s/srm/%s/%s/%s' %(
                castor_datpath, datatype, path_label, dataset))

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
    if datatype == 'data': 
        eosbase = 'root://eoscms//eos/cms/store/user/xshi/'
        ntp_label = get_label_by_version(label, 'x')
        inpath = os.path.join(eosbase, 'dat/ntp/data', ntp_label)
    else:
        inpath = os.path.join(atr.afbpath, 'dat/ntp/mc', label)

    return inpath 


def datasets(datatype, label, test=False):
    com_name = get_name_from_label(label)
    function = getattr(ntp, '%s' %com_name)
    return function(datatype, label)

    
    if 'B2KstarMuMu/RECO_100M_v1.1' in label: 
        sys.stdout.write('Using only one file!!! \n')
        rootfile = '/afs/cern.ch/user/x/xshi/work/cms/afb/dat/eos/B2KstarMuMu_GEN_SIM_HLT_100M_v1/B2KstarMuMu_RECO_100M_v1_1/ee707b714b6111f91c879d76d7f3a27d/tree_B2KstarMuMu_50_1_vbq.root'

    elif 'B2KstarMuMu/RECO_100M_v1.2' in label:
        rootfile = 'root://eoscms//eos/cms/store/user/xshi/B2KstarMuMu_GEN_SIM_HLT_100M_v1/B2KstarMuMu_RECO_100M_v1_2/ee707b714b6111f91c879d76d7f3a27d/*.root'
        if test:
            sys.stdout.write('Testing only one file. \n')
            rootfile = 'root://eoscms//eos/cms/store/user/xshi/B2KstarMuMu_GEN_SIM_HLT_100M_v1/B2KstarMuMu_RECO_100M_v1_2/ee707b714b6111f91c879d76d7f3a27d/tree_B2KstarMuMu_68_1_30S.root'
            #rootfile = '/afs/cern.ch/user/x/xshi/work/cms/afb/dat/eos/B2KstarMuMu_RECO_100M_v1_2/tree_B2KstarMuMu_68_1_30S.root'

    elif 'B2KstarMuMu/RECO_100M_v1.3' in label:
        rootpath = 'root://eoscms//eos/cms/store/user/xshi/B2KstarMuMu_GEN_SIM_HLT_100M_v1/B2KstarMuMu_RECO_100M_v1_3/a3995a8bd4369ba4e1dcae10384ea80e/'
        rootfile = rootpath + '*.root'
        if test:
            sys.stdout.write('Testing only one file. \n')
            rootfile = rootpath + 'tree_B2KstarMuMu_9_1_5l4.root'

    elif 'B2KstarMuMu/RECO_100M_v1.4' in label:
        rootpath = 'root://eoscms//eos/cms/store/user/xshi/B2KstarMuMu_GEN_SIM_HLT_100M_v1/B2KstarMuMu_RECO_100M_v1_4/faed562cb83316ea0bcafcfac78d74da/'
        rootfile = rootpath + '*.root'
        if test:
            sys.stdout.write('Testing only one file. \n')
            rootfile = rootpath + 'tree_B2KstarMuMu_21_1_YQK.root'

    elif 'B2KstarMuMu/RECO_100M_v1.5' in label:
        rootpath = 'root://eoscms//eos/cms/store/user/xshi/B2KstarMuMu_GEN_SIM_HLT_100M_v1/B2KstarMuMu_RECO_100M_v1_5/faed562cb83316ea0bcafcfac78d74da/'
        rootfile = rootpath + '*.root'
        if test:
            sys.stdout.write('Testing only one file. \n')
            rootfile = rootpath + 'tree_B2KstarMuMu_34_1_NLL.root'

    elif 'B2KstarMuMu/RECO_100M_v1.6' in label:
        rootpath = 'root://eoscms//eos/cms/store/user/xshi/B2KstarMuMu_GEN_SIM_HLT_100M_v1/B2KstarMuMu_RECO_100M_v1_6/a204b02393b0f85107bbc203039e4560/'
        rootfile = rootpath + '*.root'
        if test:
            sys.stdout.write('Testing only one file. \n')
            rootfile = rootpath + 'B0ToKstMuMu_97_1_UBI.root'

    else:
        raise NameError(label)

    return rootfile


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
    if label == 'Run2011A-PromptReco-v6_run2011v0': 
        srcdir = os.path.join(
            eosbase, 'MuOnia',
            com_name, '09dd54ed3307c6d768a6853667b85e6a')
    else:
        raise NameError(label)
    return srcdir 
    
    

