"""
Attributes for the Datasets.

"""

__author__ = "Xin Shi <Xin.Shi@cern.ch>"
__copyright__ = "Copyright (c) Xin Shi"


import os 
from tls import * 

def datasets(datatype, label):
    raise NameError(label)
    return datasets


def get_dataset_reco_brian_prod(label):
    if 'B2KstarMuMu/RECO_Brian_prod_2_1' in label:
        dataset = '/B0dKstMuMuKPi_MuMuEtaFilter_pythia6_evtgen_prod2_1/drell-B0dKstMuMuKPi_MuMuEtaFilter_pythia6_evtgen_RECO_prod2_1-d80aca2d30d1865a7fb9254b7a4518c6/USER'

    elif  'B2KstarMuMu/RECO_Brian_prod_2_2' in label:
        dataset =  '/B0dKstMuMuKPi_MuMuEtaFilter_pythia6_evtgen_prod2_2/drell-B0dKstMuMuKPi_MuMuEtaFilter_pythia6_evtgen_RECO_prod2_2-d80aca2d30d1865a7fb9254b7a4518c6/USER'

    elif 'B2KstarMuMu/RECO_Brian_prod_2_3' in label:
        dataset = '/B0dKstMuMuKPi_MuMuEtaFilter_pythia6_evtgen_prod2_3/drell-B0dKstMuMuKPi_MuMuEtaFilter_pythia6_evtgen_RECO_prod2_3-d80aca2d30d1865a7fb9254b7a4518c6/USER'
        
    elif 'B2KstarMuMu/RECO_Brian_prod_2_4' in label:
        dataset = '/B0dKstMuMuKPi_MuMuEtaFilter_pythia6_evtgen_prod2_4/drell-B0dKstMuMuKPi_MuMuEtaFilter_pythia6_evtgen_RECO_prod2_4-d80aca2d30d1865a7fb9254b7a4518c6/USER'
        
    elif 'B2KstarMuMu/RECO_Brian_prod_2_6' in label:
        dataset = '/B0dKstMuMuKPi_MuMuEtaFilter_pythia6_evtgen_prod2_6/drell-B0dKstMuMuKPi_MuMuEtaFilter_pythia6_evtgen_RECO_prod2_6-d80aca2d30d1865a7fb9254b7a4518c6/USER'
        
    elif 'B2KstarMuMu/RECO_Brian_prod_2_7' in label:
        dataset = '/B0dKstMuMuKPi_MuMuEtaFilter_pythia6_evtgen_prod2_7/drell-B0dKstMuMuKPi_MuMuEtaFilter_pythia6_evtgen_RECO_prod2_7-d80aca2d30d1865a7fb9254b7a4518c6/USER'
        
    elif 'B2KstarMuMu/RECO_Brian_prod_2_8' in label:
        dataset = '/B0dKstMuMuKPi_MuMuEtaFilter_pythia6_evtgen_prod2_8/drell-B0dKstMuMuKPi_MuMuEtaFilter_pythia6_evtgen_RECO_prod2_8-d80aca2d30d1865a7fb9254b7a4518c6/USER'
        
    elif 'B0JPsiKstMuMuKPi/RECO_Brian_prod1' in label:
        dataset = '/B0JPsiKstMuMuKPi_MuMuEtaFilter_pythia6_evtgen_prod1/drell-B0dJPsiKstMuMuKPi_MuMuEtaFilter_pythia6_evtgen_RAW2DIGI_L1Reco_RECO_reproc2-7a521434649527a1364a93839a9fbd45/USER'
        
    elif 'B0Psi2SKstMuMuKPi/RECO_Brian_prod1' in label:
        dataset = '/B0Psi2SKstMuMuKPi_MuMuEtaFilter_pythia6_evtgen_prod1/drell-B0dPsi2SKstMuMuKPi_MuMuEtaFilter_pythia6_evtgen_RAW2DIGI_L1Reco_RECO_reproc2-7a521434649527a1364a93839a9fbd45/USER'
        
    else:
        raise NameError(label)

    return dataset


def get_dataset_run2011(label):
    if 'Run2011A-May10ReReco-v1' in label:
        dataset = '/MuOnia/Run2011A-May10ReReco-v1/AOD'
    elif 'Run2011A-PromptReco-v4' in label:
        dataset = '/MuOnia/Run2011A-PromptReco-v4/AOD'
    elif 'Run2011A-PromptReco-v5' in label:
        dataset = '/MuOnia/Run2011A-PromptReco-v5/AOD'
    elif 'Run2011A-PromptReco-v6' in label:
        dataset = '/MuOnia/Run2011A-PromptReco-v6/AOD'
    elif 'Run2011B-PromptReco-v1' in label:
        dataset = '/MuOnia/Run2011B-PromptReco-v1/AOD'
    else:
        raise NameError(label)

    return dataset

def get_dataset_run2012(label):
    if 'ParkingMonitor/Run2012B-PromptReco-v1' in label:
	dataset = '/ParkingMonitor/Run2012B-PromptReco-v1/AOD'
    elif 'Run2012B-22Jan2013-v1' in label:
        dataset = '/DoubleMuParked/Run2012B-22Jan2013-v1/AOD'
    else:
        raise NameError(label)

    return dataset


def get_lumi_mask_run2011(label):
    if 'Run2011A-May10ReReco-v1' in label:
        lumi_mask = 'Cert_160404-163869_7TeV_May10ReReco_Collisions11_JSON_MuonPhys_v3.txt'
    elif 'Run2011A-PromptReco-v' in label or \
        'Run2011B-PromptReco-v' in label:
        lumi_mask = 'Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON_MuonPhys.txt'
    else:
        raise NameError(label)

    return lumi_mask


def get_lumi_mask_run2012(label):
    if 'Run2012B-22Jan2013-v1' in label:
        #lumi_mask = 'Cert_190456-190688_8TeV_PromptReco_Collisions12_JSON_MuonPhys.txt'
        #lumi_mask = 'Cert_190456-203742_8TeV_22Jan2013ReReco_Collisions12_JSON_MuonPhys.txt'
        lumi_mask =  'Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON_MuonPhys.txt'
    else:
        raise NameError(label)
    return lumi_mask

    
def jobs_created(datatype, label):
    if datatype == 'data' and label in [
            'Run2011A-May10ReReco-v1_run2011v0_2'
            ]:
        total_jobs = 1028
    elif datatype == 'data' and label in [
            'Run2011A-PromptReco-v4.6',
            'Run2011A-PromptReco-v4.10', 
            'Run2011A-PromptReco-v4.11'
            ]:
        total_jobs = 1279
    elif datatype == 'data' and label in [
            'Run2011A-PromptReco-v5_run2011v0'
            ]:
        total_jobs = 326
    elif datatype == 'data' and label in [
            'Run2011A-PromptReco-v6.5'
            ]:
        total_jobs = 474 
    elif datatype == 'data' and label in [
            'Run2011A-PromptReco-v6_run2011v0', 
            ]:
        total_jobs = 473 
    elif datatype == 'data' and label in [
            'Run2011B-PromptReco-v1.6'
            ]:
        total_jobs = 1504 

    elif datatype == 'data' and label in [
            'Run2011B-PromptReco-v1_run2011v0', 
            ]:
        total_jobs = 1578

    elif datatype == 'mc' and label in [
            'B2KstarMuMu/RECO_Brian_prod_2_1_v5', 
            'B2KstarMuMu/RECO_Brian_prod_2_2_v5', 
            'B2KstarMuMu/RECO_Brian_prod_2_3_v5', 
            'B2KstarMuMu/RECO_Brian_prod_2_4_v5', 
            'B2KstarMuMu/RECO_Brian_prod_2_8_v5',
            'B0JPsiKstMuMuKPi/RECO_Brian_prod1_v5',
            'B0Psi2SKstMuMuKPi/RECO_Brian_prod1_v5',
            'B2KstarMuMu/RECO_Brian_prod_2_1_v7', 
            'B2KstarMuMu/RECO_Brian_prod_2_2_v7', 
            'B2KstarMuMu/RECO_Brian_prod_2_3_v7', 
            'B2KstarMuMu/RECO_Brian_prod_2_4_v7', 
            'B2KstarMuMu/RECO_Brian_prod_2_8_v7',
            'B0JPsiKstMuMuKPi/RECO_Brian_prod1_v7',
            'B0Psi2SKstMuMuKPi/RECO_Brian_prod1_v7',
            ]:
        total_jobs = 101 

    elif datatype == 'mc' and label in [
            'B2KstarMuMu/RECO_Brian_prod_2_6_v5', 
            'B2KstarMuMu/RECO_Brian_prod_2_7_v5', 
            'B2KstarMuMu/RECO_Brian_prod_2_6_v7', 
            'B2KstarMuMu/RECO_Brian_prod_2_7_v7', 
            ]:
        total_jobs = 100 

    elif datatype == 'data' and label in [
            'Run2011A-May10ReReco-v1.10.6', 
            'Run2011A-May10ReReco-v1.11'
            ]:
        total_jobs = 1028
    elif label in [
            'BuToKstarJPsi-7TeV-5E5-v1_run2011v0_1', 
            ]:
        total_jobs = 20 
    else:
        raise NameError(label)
    
    jobs_created = set(range(1,total_jobs+1))
    return jobs_created

    
def get_filepath(datatype, label):
    filebase = '/store/user/xshi'
    primarydataset = ''
    psethash = ''    
    comname = get_name_from_label(label)
    if 'Run2011A-May10ReReco-v1_run2011v0_2' in label:
        primarydataset = 'MuOnia'
        psethash = 'f379a209b0c58c9fe256f1d4da070fb3'

    elif 'Run2011A-PromptReco-v5_run2011v0' in label:
        primarydataset = 'MuOnia'
        psethash = '09dd54ed3307c6d768a6853667b85e6a'
    
    elif 'Run2011A-PromptReco-v6_run2011v0' in label:
        primarydataset = 'MuOnia'
        psethash = '09dd54ed3307c6d768a6853667b85e6a'

    elif 'Run2011B-PromptReco-v1_run2011v0' in label:
        primarydataset = 'MuOnia'
        psethash = '09dd54ed3307c6d768a6853667b85e6a'

    elif ('BuToKstarJPsi-7TeV-5E5-v1_run2011v0_1' in label): 
        primarydataset = 'BuToKstarJPsi_EtaPtFilter_7TeV-pythia6-evtgen'
        psethash = 'a8d57e0034258aee57fcad0fa4e53647'

    else:
        raise NameError(label)
    
    filepath = os.path.join(filebase, primarydataset, comname, psethash)
    return filepath

def get_output_name(pyname):
    if pyname in ['btokstarmumu_Run2012.py', 
                  'btokstarmumu_Run2011A-PromptReco.py'
                  ] :
        output_name = 'BToKstarMuMu.root'
    else: 
        raise NameError(pyname)
    
    return output_name

