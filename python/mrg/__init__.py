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
    return function(args)

    # if (label == 'Run2012B-PromptReco-v2.1' or
    #     label == 'Run2011A-PromptReco-v4.10' or
    #     label == 'Run2011A-PromptReco-v5.10' or 
    #     label == 'Run2011A-PromptReco-v6.10' or 
    #     label == 'Run2011B-PromptReco-v1.10' 
    #     ) :
    #     com_name = get_name_from_label(label)
    #     function = getattr(mrg, com_name)
    #     return function(args)

    # else:
    #     raise NameError(label)


def make_src_str(files, test=False):
    if test:
        sys.stdout.write('Testing for 10 files. \n')
        fi = files[:10]
        src_str = ' '.join(fi)
    else:
        src_str = ' '.join(files)
        
    return src_str


def Run2012B_PromptReco_v2_1(args):
    test = get_options(args, 'test')
    
    label_name = get_name_from_label(label)
    datpath = os.path.join(atr.datpath, 'ntp', datatype, label_name)
    
    com_name = 'B0ToKstMuMu'
    targetname = '%s.root' % com_name
    if test:
        targetname = '%s_%s.root' % (com_name, 'test')
        
    targetfile = os.path.join(datpath, targetname)
    sourcepath = '/afs/cern.ch/user/x/xshi/eos/cms/store/user/xshi/ParkingMonitor/Run2012B-PromptReco-v2_1/46dfe96f10eeec1a9fbac57a0581cf20/'        

    dbfile = os.path.join(datpath, 'src_files.db')
    
    if os.access(dbfile, os.F_OK):
        files = get_obj_db(dbfile)
    else:
        cmd = 'ls %s' % sourcepath            
        files = get_files_from_ls(cmd)
        set_obj_db(files, dbfile)
        
    sourcefiles = make_src_str(files, test=test)
    cmd = 'hadd %s %s' %(targetfile, sourcefiles)
    proc_cmd(cmd, procdir=sourcepath)


def Run2011A_PromptReco_v4_10(args):
    eosbase = '/afs/cern.ch/user/x/xshi/eos/cms/store/user/xshi/'
    srcdir = os.path.join(eosbase, 'MuOnia', 'Run2011A_PromptReco_v4_10',
                          '08a8a08e54bfc47841ac51e8738ca48b')
    dstdir = os.path.join(eosbase, 'dat/ntp/data/Run2011v10/PromptReco_v4_10')

    if option_exists(args, '-force_update_db'):
        merge_root_files(srcdir, dstdir, force_update_db=True)
    else: 
        merge_root_files(srcdir, dstdir)


def Run2011A_PromptReco_v5_10(args):
    eosbase = '/afs/cern.ch/user/x/xshi/eos/cms/store/user/xshi/'
    srcdir = os.path.join(eosbase, 'MuOnia', 'Run2011A_PromptReco_v5_10',
                          '08a8a08e54bfc47841ac51e8738ca48b')
    dstdir = os.path.join(eosbase, 'dat/ntp/data/Run2011v10/PromptReco_v5_10')
    merge_root_files(srcdir, dstdir)
    

def Run2011A_PromptReco_v6_10(args):
    eosbase = '/afs/cern.ch/user/x/xshi/eos/cms/store/user/xshi/'
    srcdir = os.path.join(eosbase, 'MuOnia', 'Run2011A_PromptReco_v6_10',
                          '08a8a08e54bfc47841ac51e8738ca48b')
    dstdir = os.path.join(eosbase, 'dat/ntp/data/Run2011v10/PromptReco_v6_10')
    merge_root_files(srcdir, dstdir)
    

def Run2011B_PromptReco_v1_10(args):
    eosbase = '/afs/cern.ch/user/x/xshi/eos/cms/store/user/xshi/'
    srcdir = os.path.join(eosbase, 'MuOnia', 'Run2011B_PromptReco_v1_10_1',
                          '08a8a08e54bfc47841ac51e8738ca48b')
    dstdir = os.path.join(eosbase, 'dat/ntp/data/Run2011v10/B_PromptReco_v1_10_1')

    if option_exists(args, '-job'):
        job = get_option(args, '-job')
        merge_root_files(srcdir, dstdir, job=job)
    else: 
        merge_root_files(srcdir, dstdir)
    

def Run2011A_May10ReReco_v1_11(args):
    eosbase = '/afs/cern.ch/user/x/xshi/eos/cms/store/user/xshi/'
    srcdir = os.path.join(eosbase, 'MuOnia', 'Run2011A_May10ReReco_v1_11',
                          '24878667c3a4f4b7bb31cff82dfc3539')
    dstdir = os.path.join(eosbase, 
                          'dat/ntp/data/Run2011v11/May10ReReco_v1_11')
    merge_root_files(srcdir, dstdir)
    

def Run2011A_PromptReco_v4_11(args):
    eosbase = '/afs/cern.ch/user/x/xshi/eos/cms/store/user/xshi/'
    srcdir = os.path.join(eosbase, 'MuOnia', 'Run2011A_PromptReco_v4_11',
                          '22d243fc4e0d9b61112c0c52bbb77f98')
    dstdir = os.path.join(eosbase, 'dat/ntp/data/Run2011v11/PromptReco_v4_11')
    merge_root_files(srcdir, dstdir)

    
def Run2011A_PromptReco_v5_11(args):
    eosbase = '/afs/cern.ch/user/x/xshi/eos/cms/store/user/xshi/'
    srcdir = os.path.join(eosbase, 'MuOnia', 'Run2011A_PromptReco_v5_11',
                          '22d243fc4e0d9b61112c0c52bbb77f98')
    dstdir = os.path.join(eosbase, 'dat/ntp/data/Run2011v11/PromptReco_v5_11')
    merge_root_files(srcdir, dstdir)
    

def Run2011A_PromptReco_v6_11(args):
    eosbase = '/afs/cern.ch/user/x/xshi/eos/cms/store/user/xshi/'
    srcdir = os.path.join(eosbase, 'MuOnia', 'Run2011A_PromptReco_v6_11',
                          '22d243fc4e0d9b61112c0c52bbb77f98')
    dstdir = os.path.join(eosbase, 'dat/ntp/data/Run2011v11/PromptReco_v6_11')
    merge_root_files(srcdir, dstdir)


def BuToKstarMuMu_7TeV_2E7_v1_2_1(args):
    eosbase = '/afs/cern.ch/user/x/xshi/eos/cms/store/user/xshi/'
    srcdir = os.path.join(eosbase, 'BuToKstarMuMu_EtaPtFilter_7TeV-pythia6-evtgen',
                          'BuToKstarMuMu_7TeV_2E7_v1_2_1',
                          '577561224fe5925b9a09db74f385005a')
    dstdir = os.path.join(eosbase, 'dat/ntp/mc/BuToKstarMuMu/7TeV_2E7_v1_2_1')
    merge_root_files(srcdir, dstdir)

def BuToKstarMuMu_7TeV_2E7_v1_3(args):
    eosbase = '/afs/cern.ch/user/x/xshi/eos/cms/store/user/xshi/'
    srcdir = os.path.join(eosbase, 'BuToKstarMuMu_EtaPtFilter_7TeV-pythia6-evtgen',
                          'BuToKstarMuMu_7TeV_2E7_v1_3',
                          '577561224fe5925b9a09db74f385005a')
    dstdir = os.path.join(eosbase, 'dat/ntp/mc/BuToKstarMuMu/7TeV_2E7_v1_3')
    merge_root_files(srcdir, dstdir)
