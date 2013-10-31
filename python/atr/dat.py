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
    elif 'Run2012A-22Jan2013-v1' in label:
        dataset = '/DoubleMuParked/Run2012A-22Jan2013-v1/AOD'
    elif 'Run2012B-22Jan2013-v1' in label:
        dataset = '/MuOniaParked/Run2012B-22Jan2013-v1/AOD'
    elif 'Run2012C-22Jan2013-v1' in label:
        dataset = '/MuOniaParked/Run2012C-22Jan2013-v1/AOD'
    elif 'Run2012D-22Jan2013-v1' in label:
        dataset = '/MuOniaParked/Run2012D-22Jan2013-v1/AOD'
    else:
        raise NameError(label)

    return dataset


def get_lumi_mask_run2011(label):
    if 'Run2011A-May10ReReco-v1' in label:
        lumi_mask = '../data/Cert_160404-163869_7TeV_May10ReReco_Collisions11_JSON_MuonPhys_v3.txt'
    elif 'Run2011A-PromptReco-v' in label or \
        'Run2011B-PromptReco-v' in label:
        lumi_mask = '../data/Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON_MuonPhys.txt'
    else:
        raise NameError(label)

    return lumi_mask


def get_lumi_mask_run2012(label):
    if 'Run2012' in label:
        lumi_mask =  '../data/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON_MuonPhys.txt'
    else:
        raise NameError(label)
    return lumi_mask

    
def jobs_created(datatype, label):
    if datatype == 'data' and label in ['Run2011A-May10ReReco-v1_run2011v1']:
        total_jobs = 2725

    elif datatype == 'data' and label in ['Run2011A-PromptReco-v4_run2011v1']:
        total_jobs = 1291

    elif datatype == 'data' and label in ['Run2011A-PromptReco-v5_run2011v1']:
        total_jobs = 327

    elif datatype == 'data' and label in ['Run2011A-PromptReco-v6_run2011v1']:
        total_jobs = 480 

    elif datatype == 'data' and label in ['Run2011B-PromptReco-v1_run2011v1']:
        total_jobs = 1798

    elif label in ['BuToKstarJPsi-7TeV-5E5-v1_run2011v1', 
            ]:
        total_jobs = 200 
    elif label in [
            'BuToKstarMuMu-7TeV-2E7-v1_run2011v0_2', 
            ]:
        total_jobs = 1203 

    elif label in ['Run2012B-22Jan2013-v1']: 
        total_jobs = 1835

    elif label in ['Run2012C-22Jan2013-v1']: 
        total_jobs = 3305 

    elif label in ['Run2012D-22Jan2013-v1']: 
        total_jobs = 3287 

    else:
        raise NameError(label)
    
    jobs_created = set(range(1,total_jobs+1))
    return jobs_created

    
def get_filepath(datatype, label):
    filebase = '/store/user/xshi'
    primarydataset = ''
    psethash = ''    
    comname = get_name_from_label(label)
    if 'Run2011A-May10ReReco-v1_run2011v1' in label:
        primarydataset = 'MuOnia'
        psethash = '7832650ad13378e362824b793e74b5ef'

    elif 'Run2011A-PromptReco-v4_run2011v1' in label:
        primarydataset = 'MuOnia'
        psethash = '12d093a44a0e90cf594f6f76582ad92e'
    
    elif 'Run2011A-PromptReco-v5_run2011v1' in label:
        primarydataset = 'MuOnia'
        psethash = '12d093a44a0e90cf594f6f76582ad92e'
    
    elif 'Run2011A-PromptReco-v6_run2011v1' in label:
        primarydataset = 'MuOnia'
        psethash = '12d093a44a0e90cf594f6f76582ad92e'

    elif 'Run2011B-PromptReco-v1_run2011v1' in label:
        primarydataset = 'MuOnia'
        psethash = '12d093a44a0e90cf594f6f76582ad92e'

    elif ('BuToKstarJPsi-7TeV-5E5-v1_run2011v1' in label): 
        primarydataset = 'BuToKstarJPsi_EtaPtFilter_7TeV-pythia6-evtgen'
        psethash = 'a2cbd59f90c5e3694010a7612c1b81d6'
    
    elif 'Run2012' in label:
        primarydataset = ''
        psethash = ''

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

