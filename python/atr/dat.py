"""
Attributes for the Datasets.

"""

__author__ = "Xin Shi <Xin.Shi@cern.ch>"
__copyright__ = "Copyright (c) Xin Shi"


import os 
from tls import * 

def datasets(datatype, label):
    if datatype == 'dt' and label in ['5ifbv2.3', '5ifbv2.4']:
        datasets = {
            'Run2011A-May10ReReco-v1': {
                'datasetpath': '/MuOnia/Run2011A-May10ReReco-v1/AOD',
                'lumi_mask': 'Cert_160404-163869_7TeV_May10ReReco_Collisions11_JSON_MuonPhys_v3.txt',
                'pset': 'b2MuMuX_data_May10ReReco_Run2011A.py',
                'output_file': 'tree_b2MuMuX.root', 
                'number_of_jobs': '200'                
            },
            
            'Run2011A-PromptReco-v4': {
                'datasetpath': '/MuOnia/Run2011A-PromptReco-v4/AOD',
                'lumi_mask': 'Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON_MuonPhys.txt',
                'pset': 'b2MuMuX_data_promptReco_Run2011A.py',
                'output_file': 'tree_b2MuMuX.root', 
                'number_of_jobs': '300'     
            },

            'Run2011A-PromptReco-v5': {
                'datasetpath': '/MuOnia/Run2011A-PromptReco-v5/AOD',
                'lumi_mask': 'Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON_MuonPhys.txt',
                'pset': 'b2MuMuX_data_promptReco_Run2011A.py',
                'output_file': 'tree_b2MuMuX.root', 
                'number_of_jobs': '200'     
            },

            'Run2011A-PromptReco-v6': {
                'datasetpath': '/MuOnia/Run2011A-PromptReco-v6/AOD',
                'lumi_mask': 'Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON_MuonPhys.txt',
                'pset': 'b2MuMuX_data_promptReco_Run2011A.py',
                'output_file': 'tree_b2MuMuX.root', 
                'number_of_jobs': '180'     
            },

            'Run2011B-PromptReco-v1': {
                'datasetpath': '/MuOnia/Run2011B-PromptReco-v1/AOD',
                'lumi_mask': 'Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON_MuonPhys.txt',
                'pset': 'b2MuMuX_data_promptReco_Run2011A.py',
                'output_file': 'tree_b2MuMuX.root', 
                'number_of_jobs': '500'     
            },
            }

    elif datatype == 'dt' and label in ['5ifbv2.5', '5ifbv2.6']:
        datasets = {
            'Run2011A-May10ReReco-v1': {
                'datasetpath': '/MuOnia/Run2011A-May10ReReco-v1/AOD',
                'lumi_mask': 'Cert_160404-163869_7TeV_May10ReReco_Collisions11_JSON_MuonPhys_v3.txt',
                'pset': 'b2MuMuX_data_May10ReReco_Run2011A.py',
                'output_file': 'tree_b2MuMuX.root', 
                'lumis_per_job': '100'
            },
            
            'Run2011A-PromptReco-v4': {
                'datasetpath': '/MuOnia/Run2011A-PromptReco-v4/AOD',
                'lumi_mask': 'Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON_MuonPhys.txt',
                'pset': 'b2MuMuX_data_promptReco_Run2011A.py',
                'output_file': 'tree_b2MuMuX.root', 
                'lumis_per_job': '100'
            },

            'Run2011A-PromptReco-v5': {
                'datasetpath': '/MuOnia/Run2011A-PromptReco-v5/AOD',
                'lumi_mask': 'Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON_MuonPhys.txt',
                'pset': 'b2MuMuX_data_promptReco_Run2011A.py',
                'output_file': 'tree_b2MuMuX.root', 
                'lumis_per_job': '100'
            },

            'Run2011A-PromptReco-v6': {
                'datasetpath': '/MuOnia/Run2011A-PromptReco-v6/AOD',
                'lumi_mask': 'Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON_MuonPhys.txt',
                'pset': 'b2MuMuX_data_promptReco_Run2011A.py',
                'output_file': 'tree_b2MuMuX.root', 
                'lumis_per_job': '100'
            },

            'Run2011B-PromptReco-v1': {
                'datasetpath': '/MuOnia/Run2011B-PromptReco-v1/AOD',
                'lumi_mask': 'Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON_MuonPhys.txt',
                'pset': 'b2MuMuX_data_promptReco_Run2011A.py',
                'output_file': 'tree_b2MuMuX.root', 
                'lumis_per_job': '100'
            },
            }
                
    else:
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
            'Run2011A-May10ReReco-v1.9'
            ]:
        total_jobs = 1027
    elif datatype == 'data' and label in [
            'Run2011A-PromptReco-v4.6',
            'Run2011A-PromptReco-v4.10', 
            'Run2011A-PromptReco-v4.11'
            ]:
        total_jobs = 1279
    elif datatype == 'data' and label in [
            'Run2011A-PromptReco-v5.6',
            'Run2011A-PromptReco-v5.10', 
            'Run2011A-PromptReco-v5.11'
            ]:
        total_jobs = 326
    elif datatype == 'data' and label in [
            'Run2011A-PromptReco-v6.5'
            ]:
        total_jobs = 474 
    elif datatype == 'data' and label in [
            'Run2011A-PromptReco-v6.10', 
            'Run2011A-PromptReco-v6.11'
            ]:
        total_jobs = 473 
    elif datatype == 'data' and label in [
            'Run2011B-PromptReco-v1.6'
            ]:
        total_jobs = 1504 

    elif datatype == 'data' and label in [
            'Run2011B-PromptReco-v1.10.1', 
            'Run2011B-PromptReco-v1.11'
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
    else:
        raise NameError(label)
    
    jobs_created = set(range(1,total_jobs+1))
    return jobs_created

    
def get_filepath(datatype, label):
    filebase = '/store/user/xshi'
    primarydataset = ''
    psethash = ''    
    comname = get_name_from_label(label)
    if 'B2KstarMuMu/RECO_Brian_prod_2' in label:
        suffix = label.replace('B2KstarMuMu/RECO_Brian_prod_2', '')
        sub = suffix.split('_')[1]
        ver = suffix.split('_')[2]

        primarydataset = 'B0dKstMuMuKPi_MuMuEtaFilter_pythia6_evtgen_prod2_%s' % sub
        if (ver == 'v5' or ver == 'v7'):
            psethash = 'b5872eb8acdaef414564c8748b88ada2'
        else:
            raise ValueError(ver)

    elif ('B0JPsiKstMuMuKPi/RECO_Brian_prod1_v5' in label or
          'B0JPsiKstMuMuKPi/RECO_Brian_prod1_v7' in label):
          primarydataset = 'B0JPsiKstMuMuKPi_MuMuEtaFilter_pythia6_evtgen_prod1'
          psethash = 'b5872eb8acdaef414564c8748b88ada2'

    elif ('B0Psi2SKstMuMuKPi/RECO_Brian_prod1_v5' in label or
          'B0Psi2SKstMuMuKPi/RECO_Brian_prod1_v7' in label) :
          primarydataset = 'B0Psi2SKstMuMuKPi_MuMuEtaFilter_pythia6_evtgen_prod1'
          psethash = 'b5872eb8acdaef414564c8748b88ada2'

    elif 'Run2011A-May10ReReco-v1.9' in label:
        primarydataset = 'MuOnia'
        psethash = '5d47dca21a399f8cf6d7328cd7a784e8'

    elif ('Run2011A-PromptReco-v4.6' in label or
           'Run2011A-PromptReco-v5.6' in label or
           'Run2011A-PromptReco-v6.5' in label or
           'Run2011B-PromptReco-v1.6' in label):
        primarydataset = 'MuOnia'
        psethash = '623e14783eef8f7e818e025b3a872da1'

    elif 'Run2011A-May10ReReco-v1.10.6' in label:
        primarydataset = 'MuOnia'
        psethash = '18f93bc048bbabe22a8425cd3f8ea708'

    elif ('Run2011A-PromptReco-v4.10' in label or
          'Run2011A-PromptReco-v5.10' in label or 
          'Run2011A-PromptReco-v6.10' in label or
          'Run2011B-PromptReco-v1.10.1' in label): 
        primarydataset = 'MuOnia'
        psethash = '08a8a08e54bfc47841ac51e8738ca48b'
    
    elif 'Run2011A-May10ReReco-v1.11' in label:
        primarydataset = 'MuOnia'
        psethash = '24878667c3a4f4b7bb31cff82dfc3539'
    elif ('Run2011A-PromptReco-v4.11' in label or
          'Run2011A-PromptReco-v5.11' in label or 
          'Run2011A-PromptReco-v6.11' in label or
          'Run2011B-PromptReco-v1.11' in label): 
        primarydataset = 'MuOnia'
        psethash = '22d243fc4e0d9b61112c0c52bbb77f98'
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

