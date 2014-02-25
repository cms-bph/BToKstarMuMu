"""
Attributes for the Cfg files. 

"""

__author__ = "Xin Shi <Xin.Shi@cern.ch>"
__copyright__ = "Copyright (c) Xin Shi"

import os
import sys 
from tls.filetools import CfgFile 
from tls import * 
import atr
from dat import (get_dataset_run2011, get_lumi_mask_run2011,
                 get_dataset_run2012, get_lumi_mask_run2012,
                 get_output_name)

import cfg 

def crab_cfg(label, pset):
    f = CfgFile()
    f.add_section('CRAB')
    f.set('CRAB', 'jobtype', 'cmssw')
    f.set('CRAB', 'scheduler', 'remoteGlidein')
    f.set('CRAB', 'use_server', 0)

    f.add_section('CMSSW')
    f.set('CMSSW', 'datasetpath', 'None')
    f.set('CMSSW', 'pset', pset)
    f.set('CMSSW', 'total_number_of_events', 1000)
    f.set('CMSSW', 'number_of_jobs', 10)
    f.set('CMSSW', 'output_file', 'test.root')
    
    f.add_section('USER')
    f.set('USER', 'email', 'Xin.Shi@cern.ch')
    f.set('USER', 'return_data', 0)
    f.set('USER', 'copy_data', 1)
    f.set('USER', 'ui_working_dir', 'test_dir')
    f.set('USER', 'storage_element', 'T2_CH_CERN')
    f.set('USER', 'user_remote_dir', 'test_dir')
    f.set('USER', 'publish_data', 0)
    #f.set('USER', 'publish_data_name', 'test_data_name')
    #f.set('USER', 'dbs_url_for_publication',
    #      'https://cmsdbsprod.cern.ch:8443/cms_dbs_ph_analysis_02_writer/servlet/DBSServlet')
    
    f.add_section('GRID')
    f.set('GRID', 'virtual_organization', 'cms')

    com_name = get_name_from_label(label, ver='x.')
    com_name = filter_com_func_name(com_name)
    function = getattr(cfg, com_name)
    f, cfg_file = function(f, label)
    f.output(cfg_file)
    return cfg_file


def filter_com_func_name(com_name):
    if ('B2KstarMuMu_RECO_Brian_prod_2_' in com_name or
        'B0JPsiKstMuMuKPi_RECO_Brian_prod1' in com_name or
        'B0Psi2SKstMuMuKPi_RECO_Brian_prod1' in com_name):
        ver = com_name.split('_')[-1]
        com_name = 'B2KstarMuMu_RECO_Brian_prod_2_1_%s' %ver
        
    return com_name
    
def Run2011A_May10ReReco_v1_run2011v1(f, label):
    com_name = get_name_from_label(label)
    cfg_name = 'crab_%s.cfg' % com_name
    cmssw = 'CMSSW_4_2_8_patch7'
    src_name = 'BphAna/BToKstarMuMu_run2011v1' 
    cfg_path = os.path.join(atr.afbpath, 'rel', cmssw, 'src', src_name, 'python/')
    cfg_file = os.path.join(cfg_path, cfg_name)
    
    f.set('CMSSW', 'pset', 'btokstarmumu_Run2011A-May10ReReco-v1.py')
    f.set('CMSSW', 'total_number_of_lumis', -1)
    f.set('CMSSW', 'lumis_per_job', 1)
    f.remove_option('CMSSW', 'total_number_of_events')
    f.remove_option('CMSSW', 'number_of_jobs')
    f.set('CMSSW', 'output_file', 'BToKstarMuMu.root')
    f.set('CMSSW', 'datasetpath', get_dataset_run2011(label))
    f.set('CMSSW', 'lumi_mask', get_lumi_mask_run2011(label))
    f.set('USER', 'user_remote_dir', com_name)
    f.set('USER', 'ui_working_dir', 'crab_%s' %com_name)
    f.set('USER', 'publish_data_name', com_name)
    f.set('CRAB', 'scheduler', 'remoteGlidein')
    f.set('CRAB', 'use_server', 0)
    return f, cfg_file



def Run2011A_PromptReco_v4_run2011v1(f, label):
    f, cfg_file = Run2011A_May10ReReco_v1_run2011v1(f, label)
    return f, cfg_file

def Run2011A_PromptReco_v5_run2011v1(f, label):
    return Run2011A_PromptReco_v4_run2011v1(f, label)

def Run2011A_PromptReco_v6_run2011v1(f, label):
    return Run2011A_PromptReco_v4_run2011v1(f, label)

def Run2011B_PromptReco_v1_run2011v1(f, label):
    f, cfg_file = Run2011A_PromptReco_v4_run2011v1(f, label)
    if 'run2011v1.1' in label:
        f.remove_option('CMSSW', 'lumis_per_job')
        f.set('CMSSW', 'number_of_jobs', 5000)

    return f, cfg_file


def BuToKstarJPsi_7TeV_5E5_v1_run2011v1(f, label):
    com_name = get_name_from_label(label)
    cfg_name = 'crab_%s.cfg' % com_name
    cmssw = 'CMSSW_4_2_8_patch7'
    src_name = 'BphAna/BToKstarMuMu_run2011v1'
    cfg_path = os.path.join(atr.afbpath, 'rel', cmssw, 'src', src_name, 'python')
    cfg_file = os.path.join(cfg_path, cfg_name)

    f.set('CMSSW', 'datasetpath', '/BuToKstarJPsi_EtaPtFilter_7TeV-pythia6-evtgen/Fall11-HLTMuonia_PU_S6_START42_V14B-v1/AODSIM')
    f.set('CMSSW', 'pset', 'btokstarmumu_MC.py')
    f.set('CMSSW', 'total_number_of_events', get_number_from_label(label))
    f.set('CMSSW', 'number_of_jobs', 200)
    f.set('CMSSW', 'output_file', 'BToKstarMuMu.root')
    f.set('USER', 'user_remote_dir', com_name)
    f.set('USER', 'ui_working_dir', 'crab_%s' %com_name)
    f.set('USER', 'publish_data_name', com_name)
    f.set('CRAB', 'scheduler', 'remoteGlidein')
    f.set('CRAB', 'use_server', 0)
    return f, cfg_file



def BuToKstarMuMu_7TeV_2E7_v1_run2011v1(f, label):
    f, cfg_file = BuToKstarJPsi_7TeV_5E5_v1_run2011v1(f, label)
    f.set('CMSSW', 'datasetpath', '/BuToKstarMuMu_EtaPtFilter_7TeV-pythia6-evtgen/Fall11-HLTMuonia_PU_S6_START42_V14B-v1/AODSIM')
    f.set('CMSSW', 'number_of_jobs', 1200)
    return f, cfg_file

def Run2012A_22Jan2013_v1(f, label):
    com_name = get_name_from_label(label)
    cfg_name = 'crab_%s.cfg' % com_name
    #cmssw = 'CMSSW_4_2_8_patch7'
    #src_name = 'BphAna/BToKstarMuMu_run2011v1' 
    cfg_path = os.getcwd()
    cfg_file = os.path.join(cfg_path, cfg_name)
    
    #f.set('CMSSW', 'pset', 'btokstarmumu_Run2011A-May10ReReco-v1.py')
    f.set('CMSSW', 'total_number_of_lumis', -1)
    f.set('CMSSW', 'lumis_per_job', 20)
    f.remove_option('CMSSW', 'total_number_of_events')
    f.remove_option('CMSSW', 'number_of_jobs')
    f.set('CMSSW', 'output_file', 'BToKstarMuMu.root')
    f.set('CMSSW', 'datasetpath', get_dataset_run2012(label))
    f.set('CMSSW', 'lumi_mask', get_lumi_mask_run2012(label))
    f.set('USER', 'user_remote_dir', com_name)
    f.set('USER', 'ui_working_dir', 'crab_%s' %com_name)
    return f, cfg_file

def Run2012B_22Jan2013_v1(f, label):
    return Run2012A_22Jan2013_v1(f, label)

def Run2012C_22Jan2013_v1(f, label):
    return Run2012A_22Jan2013_v1(f, label)

def Run2012D_22Jan2013_v1(f, label):
    return Run2012A_22Jan2013_v1(f, label)
