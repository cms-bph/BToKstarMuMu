"""
Module for Merging ROOT files 

"""

__author__ = "Xin Shi <Xin.Shi@cern.ch>"
__copyright__ = "Copyright (c) Xin Shi"

import sys 
from tls import * 
import atr
import mrg 

def main(args):
    
    if args[0] == 'SingleBuToKstarMuMu':
        import SingleBuToKstarMuMu
        return SingleBuToKstarMuMu.main(args[1:])
    
    datatype = args[0]
    label = args[1]
    
    com_name = get_name_from_label(label)
    function = getattr(mrg, com_name)
    return function(com_name)

 
def make_src_str(files, test=False):
    if test:
        sys.stdout.write('Testing for 10 files. \n')
        fi = files[:10]
        src_str = ' '.join(fi)
    else:
        src_str = ' '.join(files)
    return src_str


def BuToKstarJPsi_7TeV_5E5_v1_run2011v1(com_name):
    eosbase = '/afs/cern.ch/user/x/xshi/eos/cms/store/user/xshi/'
    srcdir = os.path.join(
        eosbase, 'BuToKstarJPsi_EtaPtFilter_7TeV-pythia6-evtgen',
        com_name, 'a2cbd59f90c5e3694010a7612c1b81d6')
    dstdir = os.path.join(eosbase, 'dat/ntp/mc', com_name)
    merge_root_files(srcdir, dstdir)


def Run2011A_May10ReReco_v1_run2011v1(com_name):
    eosbase = '/afs/cern.ch/user/x/xshi/eos/cms/store/user/xshi/'
    srcdir = os.path.join(
        eosbase, 'MuOnia',
        com_name, '7832650ad13378e362824b793e74b5ef')
    dstdir = os.path.join(eosbase, 'dat/ntp/data', com_name)
    merge_root_files(srcdir, dstdir)


def Run2011A_PromptReco_v4_run2011v1(com_name):
    eosbase = '/afs/cern.ch/user/x/xshi/eos/cms/store/user/xshi/'
    srcdir = os.path.join(
        eosbase, 'MuOnia',
        com_name, '12d093a44a0e90cf594f6f76582ad92e')
    dstdir = os.path.join(eosbase, 'dat/ntp/data', com_name)
    merge_root_files(srcdir, dstdir)

def Run2011A_PromptReco_v5_run2011v1(com_name):
    Run2011A_PromptReco_v4_run2011v1(com_name)


def Run2011A_PromptReco_v6_run2011v1(com_name):
    Run2011A_PromptReco_v5_run2011v1(com_name)


def Run2011B_PromptReco_v1_run2011v1(com_name):
    Run2011A_PromptReco_v6_run2011v1(com_name)


def BuToKstarMuMu_7TeV_2E7_v1_run2011v1(com_name):
    eosbase = '/afs/cern.ch/user/x/xshi/eos/cms/store/user/xshi/'
    srcdir = os.path.join(
        eosbase, 'BuToKstarMuMu_EtaPtFilter_7TeV-pythia6-evtgen',
        com_name, 'a8d57e0034258aee57fcad0fa4e53647')
    dstdir = os.path.join(eosbase, 'dat/ntp/mc', com_name)
    merge_root_files(srcdir, dstdir)

def Run2012A_22Jan2013_v1(com_name):
    eosbase = '/afs/cern.ch/user/x/xshi/eos/cms/store/user/xshi/'
    srcdir = os.path.join(eosbase, com_name)
    dstdir = os.path.join(atr.datpath, 'ntp/data', atr.version, com_name)
    merge_root_files(srcdir, dstdir)

def Run2012B_22Jan2013_v1(com_name):
    return Run2012A_22Jan2013_v1(com_name)

def Run2012C_22Jan2013_v1(com_name):
    return Run2012A_22Jan2013_v1(com_name)

def Run2012D_22Jan2013_v1(com_name):
    return Run2012A_22Jan2013_v1(com_name)

