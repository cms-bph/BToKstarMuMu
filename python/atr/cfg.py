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
                 get_output_name)

import cfg 

def crab_cfg(label):
    f = CfgFile()
    f.add_section('CRAB')
    f.set('CRAB', 'jobtype', 'cmssw')
    f.set('CRAB', 'scheduler', 'glidein')
    f.set('CRAB', 'use_server', 1)

    f.add_section('CMSSW')
    f.set('CMSSW', 'datasetpath', 'None')
    f.set('CMSSW', 'pset', 'test.py')
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
    f.set('USER', 'publish_data', 1)
    f.set('USER', 'publish_data_name', 'test_data_name')
    f.set('USER', 'dbs_url_for_publication',
          'https://cmsdbsprod.cern.ch:8443/cms_dbs_ph_analysis_02_writer/servlet/DBSServlet')
    
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
    
 

def Run2011A_May10ReReco_v1_7(f, label):
    com_name = get_name_from_label(label)
    cfg_name = 'crab_%s.cfg' % com_name
    cmssw = 'CMSSW_4_2_8_patch7'
    src_name = 'BphAna/B0KstMuMu_V00_07_01'
    cfg_path = os.path.join(atr.afbpath, 'rel', cmssw, 'src', src_name, 'python/')
    cfg_file = os.path.join(cfg_path, cfg_name)
    
    f.set('CMSSW', 'pset', 'Run2011A-May10ReReco-v1_3.py')
    f.set('CMSSW', 'total_number_of_lumis', -1)
    f.set('CMSSW', 'lumis_per_job', 20)
    f.remove_option('CMSSW', 'total_number_of_events')
    f.remove_option('CMSSW', 'number_of_jobs')
    f.set('CMSSW', 'output_file', 'B0ToKstMuMu.root')
    f.set('CMSSW', 'datasetpath', get_dataset_run2011(label))
    f.set('CMSSW', 'lumi_mask', get_lumi_mask_run2011(label))
    f.set('USER', 'user_remote_dir', com_name)
    f.set('USER', 'ui_working_dir', 'crab_%s' %com_name)
    f.set('USER', 'publish_data_name', com_name)
    f.set('USER', 'additional_input_files', os.path.join(
        atr.afbpath, 'rel', cmssw, 'src', src_name, 'python', 'ParameterFile.txt'))

    return f, cfg_file


def Run2011A_May10ReReco_v1_8(f, label):
    f, cfg_file = Run2011A_May10ReReco_v1_7(f, label)
    f.set('CRAB', 'scheduler', 'remoteGlidein')
    f.set('CRAB', 'use_server', 0)
    return f, cfg_file


def Run2011A_May10ReReco_v1_10(f, label):
    f, cfg_file = Run2011A_May10ReReco_v1_8(f, label)
    com_name = get_name_from_label(label)
    cfg_name = 'crab_%s.cfg' % com_name
    cmssw = 'CMSSW_4_2_8_patch7'
    src_name = 'BphAna/BToKstarMuMu'
    cfg_path = os.path.join(atr.afbpath, 'rel', cmssw, 'src', src_name, 'python')
    cfg_file = os.path.join(cfg_path, cfg_name)
    
    f.set('CMSSW', 'output_file', 'BToKstarMuMu.root')
    f.set('CMSSW', 'pset', 'btokstarmumu_Run2011A-May10ReReco-v1_10.py')
    f.remove_option('USER', 'additional_input_files')
    return f, cfg_file


def Run2011A_May10ReReco_v1_11(f, label):
    f, cfg_file = Run2011A_May10ReReco_v1_10(f, label)
    f.set('CMSSW', 'pset', 'btokstarmumu_Run2011A-May10ReReco-v1_11.py')

    com_name = get_name_from_label(label)
    cfg_name = 'crab_%s.cfg' % com_name
    cmssw = 'CMSSW_4_2_8_patch7'
    src_name = 'BphAna/BToKstarMuMu_V00_00_02'
    cfg_path = os.path.join(atr.afbpath, 'rel', cmssw, 'src', src_name, 'python')
    cfg_file = os.path.join(cfg_path, cfg_name)
    return f, cfg_file


def Run2011A_May10ReReco_v1_run2011v0(f, label):
    com_name = get_name_from_label(label)
    cfg_name = 'crab_%s.cfg' % com_name
    cmssw = 'CMSSW_4_2_8_patch7'
    src_name = 'BphAna/BToKstarMuMu_run2011v0' 
    cfg_path = os.path.join(atr.afbpath, 'rel', cmssw, 'src', src_name, 'python/')
    cfg_file = os.path.join(cfg_path, cfg_name)
    
    f.set('CMSSW', 'pset', 'btokstarmumu_Run2011A-May10ReReco-v1.py')
    f.set('CMSSW', 'total_number_of_lumis', -1)
    f.set('CMSSW', 'lumis_per_job', 20)
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


def Run2011A_PromptReco_v4_4(f, label):
    f, cfg_file = Run2011A_May10ReReco_v1_7(f, label)
    f.set('CMSSW', 'pset', 'Run2011A-PromptReco-v4_3.py')
    return f, cfg_file


def Run2011A_PromptReco_v4_5(f, label):
    f, cfg_file = Run2011A_PromptReco_v4_4(f, label)
    f.set('CRAB', 'scheduler', 'remoteGlidein')
    f.set('CRAB', 'use_server', 0)
    return f, cfg_file


def Run2011A_PromptReco_v4_10(f, label):
    f, cfg_file = Run2011A_May10ReReco_v1_10(f, label)
    f.set('CMSSW', 'pset', 'btokstarmumu_Run2011A-PromptReco-v4_10.py')
    return f, cfg_file

def Run2011A_PromptReco_v4_11(f, label):
    f, cfg_file = Run2011A_May10ReReco_v1_11(f, label)
    f.set('CMSSW', 'pset', 'btokstarmumu_Run2011A-PromptReco-v4_11.py')
    return f, cfg_file

def Run2011A_PromptReco_v4_run2011v0(f, label):
    f, cfg_file = Run2011A_May10ReReco_v1_run2011v0(f, label)
    f.set('CMSSW', 'pset', 'btokstarmumu_Run2011A-PromptReco.py')
    return f, cfg_file

def Run2011A_PromptReco_v5_run2011v0(f, label):
    return Run2011A_PromptReco_v4_run2011v0(f, label)

def Run2011A_PromptReco_v6_run2011v0(f, label):
    return Run2011A_PromptReco_v4_run2011v0(f, label)

def Run2011B_PromptReco_v1_run2011v0(f, label):
    return Run2011A_PromptReco_v4_run2011v0(f, label)

def Run2011A_PromptReco_v5_10(f, label):
    return Run2011A_PromptReco_v4_10(f, label)

def Run2011A_PromptReco_v6_10(f, label):
    return Run2011A_PromptReco_v4_10(f, label)

def Run2011B_PromptReco_v1_10(f, label):
    return Run2011A_PromptReco_v4_10(f, label)

def Run2011A_PromptReco_v5_11(f, label):
    return Run2011A_PromptReco_v4_11(f, label)

def Run2011A_PromptReco_v6_11(f, label):
    return Run2011A_PromptReco_v4_11(f, label)

def Run2011B_PromptReco_v1_11(f, label):
    f, cfg_file = Run2011A_PromptReco_v4_11(f, label)
    f.set('CMSSW', 'lumis_per_job', 10)
    return f, cfg_file


def Run2011B_PromptReco_v1_12(f, label):
    f, cfg_file = Run2011A_PromptReco_v4_11(f, label)
    com_name = get_name_from_label(label)
    cfg_name = 'crab_%s.cfg' % com_name
    cmssw = 'CMSSW_4_2_8_patch7'
    src_name = 'BphAna/BToKstarMuMu_Run2011v12'
    cfg_path = os.path.join(atr.afbpath, 'rel', cmssw, 'src', src_name, 'python')
    cfg_file = os.path.join(cfg_path, cfg_name)
    
    f.set('CMSSW', 'pset', 'btokstarmumu_Run2011A-PromptReco.py')
    f.set('CMSSW', 'output_file', get_output_name('btokstarmumu_Run2011A-PromptReco.py'))

    return f, cfg_file

    

def ParkingMonitor_Run2012B_PromptReco_v1_1(f, label):
    com_name = get_name_from_label(label)
    cfg_name = 'crab_%s.cfg' % com_name
    cmssw = 'CMSSW_5_3_8_patch3'
    src_name = 'BphAna/BToKstarMuMu'
    cfg_path = os.path.join(atr.afbpath, 'rel', cmssw, 'src', src_name, 'python/')
    cfg_file = os.path.join(cfg_path, cfg_name)
    
    f.set('CMSSW', 'pset', 'btokstarmumu_Run2012A-PromptReco.py')
    f.set('CMSSW', 'total_number_of_lumis', -1)
    f.set('CMSSW', 'lumis_per_job', 20)
    f.remove_option('CMSSW', 'total_number_of_events')
    f.remove_option('CMSSW', 'number_of_jobs')
    f.set('CMSSW', 'output_file', 'BToKstarMuMu.root')
    f.set('CMSSW', 'datasetpath', get_dataset_run2012(label))
    f.set('USER', 'user_remote_dir', com_name)
    f.set('USER', 'ui_working_dir', 'crab_%s' %com_name)
    f.set('USER', 'publish_data_name', com_name)
    f.set('CRAB', 'scheduler', 'remoteGlidein')
    f.set('CRAB', 'use_server', 0)
    return f, cfg_file


def DoubleMuParked_Run2012B_22Jan2013_v1_1(f, label):
    f, cfg_file = ParkingMonitor_Run2012B_PromptReco_v1_1(f, label)
    f.set('CMSSW', 'lumi_mask', get_lumi_mask_run2012(label))
    return f, cfg_file
   

def B2KstarMuMu_RECO_Brian_prod_2_1_v5(f, label):
    com_name = get_name_from_label(label)
    cfg_name = 'crab_%s.cfg' % com_name
    cmssw = 'CMSSW_4_2_8_patch7'
    src_name = 'BphAna/B0KstMuMu_V00_07_01'
    cfg_path = os.path.join(atr.afbpath, 'rel', cmssw, 'src', src_name, 'python/')
    cfg_file = os.path.join(cfg_path, cfg_name)

    f.set('CRAB', 'scheduler', 'remoteGlidein')
    f.set('CRAB', 'use_server', 0)
    f.set('CMSSW', 'pset', 'B2KstarMuMu_RECO_Brian_v3.py')    
    f.set('CMSSW', 'total_number_of_events', -1)
    f.set('CMSSW', 'number_of_jobs', 100)
    f.set('CMSSW', 'output_file', 'B0ToKstMuMu.root')
    f.set('CMSSW', 'dbs_url', 'https://cmsdbsprod.cern.ch:8443/cms_dbs_ph_analysis_02_writer/servlet/DBSServlet')
    f.set('CMSSW', 'datasetpath', get_dataset_reco_brian_prod(label))
    f.set('USER', 'user_remote_dir', com_name)
    f.set('USER', 'ui_working_dir', 'crab_%s' %com_name)
    f.set('USER', 'publish_data_name', com_name)
    f.set('USER', 'additional_input_files', os.path.join(
        atr.afbpath, 'rel', cmssw, 'src', src_name, 'python', 'ParameterFile.txt'))
    return f, cfg_file


def B2KstarMuMu_RECO_Brian_prod_2_1_v6(f, label):
    com_name = get_name_from_label(label)
    cfg_name = 'crab_%s.cfg' % com_name
    cmssw = 'CMSSW_4_2_8_patch7'
    src_name = 'BphAna/B0KstMuMu_V00_08_00'
    cfg_path = os.path.join(atr.afbpath, 'rel', cmssw, 'src', src_name, 'python/')
    cfg_file = os.path.join(cfg_path, cfg_name)

    f.set('CRAB', 'scheduler', 'remoteGlidein')
    f.set('CRAB', 'use_server', 0)
    f.set('CMSSW', 'pset', 'B2KstarMuMu_RECO_Brian_v3.py')    
    f.set('CMSSW', 'total_number_of_events', -1)
    f.set('CMSSW', 'number_of_jobs', 100)
    f.set('CMSSW', 'output_file', 'B0ToKstMuMu.root')
    f.set('CMSSW', 'dbs_url', 'https://cmsdbsprod.cern.ch:8443/cms_dbs_ph_analysis_02_writer/servlet/DBSServlet')
    f.set('CMSSW', 'datasetpath', get_dataset_reco_brian_prod(label))
    f.set('USER', 'user_remote_dir', com_name)
    f.set('USER', 'ui_working_dir', 'crab_%s' %com_name)
    f.set('USER', 'publish_data_name', com_name)
    f.set('USER', 'additional_input_files', os.path.join(
        atr.afbpath, 'rel', cmssw, 'src', src_name, 'python', 'ParameterFile.txt'))
    return f, cfg_file


def B2KstarMuMu_RECO_Brian_prod_2_1_v7(f, label):
    com_name = get_name_from_label(label)
    cfg_name = 'crab_%s.cfg' % com_name
    cmssw = 'CMSSW_4_2_8_patch7'
    src_name = 'BphAna/B0KstMuMu_V00_08_01'
    cfg_path = os.path.join(atr.afbpath, 'rel', cmssw, 'src', src_name, 'python/')
    cfg_file = os.path.join(cfg_path, cfg_name)

    f, tmp = B2KstarMuMu_RECO_Brian_prod_2_1_v6(f, label)
    f.set('USER', 'additional_input_files', os.path.join(
        atr.afbpath, 'rel', cmssw, 'src', src_name, 'python', 'ParameterFile.txt'))
    return f, cfg_file

def BuToKstarMuMu_GEN_SIM_1M_v1(f, label):
    cmssw = 'CMSSW_4_2_9_HLT1_bphpatch4'
    com_name = get_name_from_label(label)
    cfg_name = 'crab_%s.cfg' % com_name
    cfg_path = os.path.join(atr.afbpath, 'rel', cmssw, 'src')
    cfg_file = os.path.join(cfg_path, cfg_name)

    f.set('CMSSW', 'pset', '%s.py' %com_name)
    f.set('CMSSW', 'total_number_of_events', get_number_from_label(label))
    f.set('CMSSW', 'number_of_jobs', 100)
    f.set('CMSSW', 'output_file', 'BuToKstarMuMu_7TeV_PYTHIA6.root')
    f.set('USER', 'user_remote_dir', com_name)
    f.set('USER', 'ui_working_dir', 'crab_%s' %com_name)
    f.set('USER', 'publish_data_name', com_name)
    f.set('CRAB', 'scheduler', 'remoteGlidein')
    f.set('CRAB', 'use_server', 0)

    return f, cfg_file

def BuToKstarMuMu_HLT_1M_v1(f, label):
    f, cfg_file = BuToKstarMuMu_GEN_SIM_1M_v1(f, label)
    f.set('CMSSW', 'number_of_jobs', 1)
    f.set('CMSSW', 'dbs_url', 'https://cmsdbsprod.cern.ch:8443/\
cms_dbs_ph_analysis_02_writer/servlet/DBSServlet')
    f.set('CMSSW', 'datasetpath', '/BuToKstarMuMu_GEN_SIM_1M_v1/\
xshi-BuToKstarMuMu_GEN_SIM_1M_v1-16e627c3ed3fc8fde633727fdddd8391/USER')

    return f, cfg_file

def BuToKstarMuMu_RECO_1M_v1(f, label):
    f, cfg_file = BuToKstarMuMu_HLT_1M_v1(f, label)
    f.set('CMSSW', 'datasetpath', '/BuToKstarMuMu_GEN_SIM_1M_v1/\
xshi-BuToKstarMuMu_HLT_1M_v1-af3f88d919223c97e5a5e67b7294db4e/USER')

    return f, cfg_file
    
    
def BuToKstarMuMu_GEN_SIM_10M_v1_1(f, label):
    return BuToKstarMuMu_GEN_SIM_1M_v1(f, label)

def BuToKstarMuMu_HLT_10M_v1(f, label):
    f, cfg_file = BuToKstarMuMu_HLT_1M_v1(f, label)
    f.set('CMSSW', 'datasetpath', '/BuToKstarMuMu_GEN_SIM_10M_v1_1/\
xshi-BuToKstarMuMu_GEN_SIM_10M_v1_1-16e627c3ed3fc8fde633727fdddd8391/USER')
    f.set('CMSSW', 'number_of_jobs', 4)
    
    return f, cfg_file
    
def BuToKstarMuMu_RECO_10M_v1(f, label):
    f, cfg_file = BuToKstarMuMu_HLT_10M_v1(f, label)
    f.set('CMSSW', 'datasetpath', '/BuToKstarMuMu_GEN_SIM_10M_v1_1/\
xshi-BuToKstarMuMu_HLT_10M_v1-af3f88d919223c97e5a5e67b7294db4e/USER')
    return f, cfg_file
    
def BuToKstarMuMu_7TeV_2E7_v1_1(f, label):
    com_name = get_name_from_label(label)
    cfg_name = 'crab_%s.cfg' % com_name
    cmssw = 'CMSSW_4_2_8_patch7'
    src_name = 'BphAna/BToKstarMuMu_v1'
    cfg_path = os.path.join(atr.afbpath, 'rel', cmssw, 'src', src_name, 'python')
    cfg_file = os.path.join(cfg_path, cfg_name)
 
    f.set('CMSSW', 'datasetpath', '/BuToKstarMuMu_EtaPtFilter_7TeV-pythia6-evtgen/Fall11-HLTMuonia_PU_S6_START42_V14B-v1/AODSIM')
    #f.set('CMSSW', 'datasetpath', '/BuToKstarMuMu_EtaPtFilter_7TeV-pythia6-evtgen/Fall11-HLTMuonia_PU_S6_START44_V9B-v1/AODSIM')

    f.set('CMSSW', 'pset', 'btokstarmumu_MC.py')
    f.set('CMSSW', 'total_number_of_events', get_number_from_label(label))
    f.set('CMSSW', 'number_of_jobs', 2000)
    f.set('CMSSW', 'output_file', 'BToKstarMuMu.root')
    f.set('USER', 'user_remote_dir', com_name)
    f.set('USER', 'ui_working_dir', 'crab_%s' %com_name)
    f.set('USER', 'publish_data_name', com_name)
    f.set('CRAB', 'scheduler', 'remoteGlidein')
    f.set('CRAB', 'use_server', 0)

    return f, cfg_file


def BuToKstarMuMu_7TeV_2E7_v1_2(f, label):
    com_name = get_name_from_label(label)
    cfg_name = 'crab_%s.cfg' % com_name
    cmssw = 'CMSSW_4_2_8_patch7'
    src_name = 'BphAna/BToKstarMuMu_v2'
    cfg_path = os.path.join(atr.afbpath, 'rel', cmssw, 'src', src_name, 'python')
    cfg_file = os.path.join(cfg_path, cfg_name)

    f.set('CMSSW', 'datasetpath', '/BuToKstarMuMu_EtaPtFilter_7TeV-pythia6-evtgen/Fall11-HLTMuonia_PU_S6_START42_V14B-v1/AODSIM')

    f.set('CMSSW', 'pset', 'btokstarmumu_MC.py')
    f.set('CMSSW', 'total_number_of_events', get_number_from_label(label))
    f.set('CMSSW', 'number_of_jobs', 400)
    f.set('CMSSW', 'output_file', 'BToKstarMuMu.root')
    f.set('USER', 'user_remote_dir', com_name)
    f.set('USER', 'ui_working_dir', 'crab_%s' %com_name)
    f.set('USER', 'publish_data_name', com_name)
    f.set('CRAB', 'scheduler', 'remoteGlidein')
    f.set('CRAB', 'use_server', 0)
    return f, cfg_file

def BuToKstarMuMu_7TeV_2E7_v1_3(f, label):
    f, cfg_file = BuToKstarMuMu_7TeV_2E7_v1_2(f, label)
    com_name = get_name_from_label(label)
    cfg_name = 'crab_%s.cfg' % com_name
    cmssw = 'CMSSW_4_2_8_patch7'
    src_name = 'BphAna/BToKstarMuMu'
    cfg_path = os.path.join(atr.afbpath, 'rel', cmssw, 'src', src_name, 'python')
    cfg_file = os.path.join(cfg_path, cfg_name)
    return f, cfg_file

    
def Run2012B_22Jan2013_v1_1(f, label):
    com_name = get_name_from_label(label)
    cfg_name = 'crab_%s.cfg' % com_name
    cmssw = 'CMSSW_5_3_9_patch3'
    src_name = 'BphAna/BToKstarMuMu_Run2012v1'
    cfg_path = os.path.join(atr.afbpath, 'rel', cmssw, 'src', src_name, 'python/')
    cfg_file = os.path.join(cfg_path, cfg_name)
    
    f.set('CMSSW', 'pset', 'btokstarmumu_Run2012.py')
    f.set('CMSSW', 'total_number_of_lumis', -1)
    f.set('CMSSW', 'lumis_per_job', 20)
    f.remove_option('CMSSW', 'total_number_of_events')
    f.remove_option('CMSSW', 'number_of_jobs')
    #f.set('CMSSW', 'output_file', 'BToKstarMuMu.root')
    f.set('CMSSW', 'output_file', get_output_name('btokstarmumu_Run2012.py'))
    f.set('CMSSW', 'datasetpath', get_dataset_run2012(label))
    f.set('CMSSW', 'lumi_mask', get_lumi_mask_run2012(label))
    f.set('USER', 'user_remote_dir', com_name)
    f.set('USER', 'ui_working_dir', 'crab_%s' %com_name)
    f.set('USER', 'publish_data_name', com_name)
    f.set('CRAB', 'scheduler', 'remoteGlidein')
    f.set('CRAB', 'use_server', 0)

    return f, cfg_file

def BuToKstarJPsi_7TeV_5E5_v1_run2011v0(f, label):
    com_name = get_name_from_label(label)
    cfg_name = 'crab_%s.cfg' % com_name
    cmssw = 'CMSSW_4_2_8_patch7'
    src_name = 'BphAna/BToKstarMuMu_run2011v0'
    cfg_path = os.path.join(atr.afbpath, 'rel', cmssw, 'src', src_name, 'python')
    cfg_file = os.path.join(cfg_path, cfg_name)

    f.set('CMSSW', 'datasetpath', '/BuToKstarJPsi_EtaPtFilter_7TeV-pythia6-evtgen/Fall11-HLTMuonia_PU_S6_START42_V14B-v1/AODSIM')
    f.set('CMSSW', 'pset', 'btokstarmumu_MC.py')
    f.set('CMSSW', 'total_number_of_events', get_number_from_label(label))
    f.set('CMSSW', 'number_of_jobs', 20)
    f.set('CMSSW', 'output_file', 'BToKstarMuMu.root')
    f.set('USER', 'user_remote_dir', com_name)
    f.set('USER', 'ui_working_dir', 'crab_%s' %com_name)
    f.set('USER', 'publish_data_name', com_name)
    f.set('CRAB', 'scheduler', 'remoteGlidein')
    f.set('CRAB', 'use_server', 0)
    return f, cfg_file

def BuToKstarMuMu_7TeV_2E7_v1_run2011v0(f, label):
    f, cfg_file = BuToKstarJPsi_7TeV_5E5_v1_run2011v0(f, label)
    f.set('CMSSW', 'datasetpath', '/BuToKstarMuMu_EtaPtFilter_7TeV-pythia6-evtgen/Fall11-HLTMuonia_PU_S6_START42_V14B-v1/AODSIM')
    f.set('CMSSW', 'number_of_jobs', 1200)
    return f, cfg_file

