"""
Module for producing NTuple (Before Selection) 

"""

__author__ = "Xin Shi <Xin.Shi@cern.ch>"
__copyright__ = "Copyright (c) Xin Shi"

import sys
import atr 
from tls import *
from time import time 


def main(args):
    datatype = args[0]
    label = args[1]
    test = option_exists(args, '-t')
    pbar = option_exists(args, '-p')
    batch = option_exists(args, '-b')
    job_opt = option_exists(args, '-job')
    grid = option_exists(args, '-grid')

    stage = 'srm'
    
    # if label in ['B2KstarMuMu/GEN_1M_v1',
    #              'B2KstarMuMu/RECO_1M_v2',
    #              'B2KstarMuMu/RECO_1M_v2.4',
    #              'B2KstarMuMu/RECO_Brian_v1',
    #              'B2KstarMuMu/RECO_Brian_prod_2_1_v1',
    #              'B2KstarMuMu/RECO_Brian_prod_2_2_v1',
    #              'B2KstarMuMu/RECO_Brian_prod_2_6_v1',
    #              'B2KstarMuMu/RECO_Brian_prod_2_7_v1',
    #              'B2KstarMuMu/RECO_Brian_prod_2_8_v1',
    #              'B2KstarMuMu/RECO_Brian_prod_2_3_v1',
    #              'B2KstarMuMu/RECO_Brian_prod_2_4_v1']: 
    #     stage = 'srm'
    #     #size_max = 100 # do not group jobs
        
    # else:
    #     raise NameError(label)


    if job_opt == 'create':
        dbfile = create_jobs_db(stage, datatype, label, test=test)
        return

    elif job_opt == 'submit':
        jobs = get_jobs(stage, datatype, label)        
        range_ = option_exists(args, 'submit')
        if range_ == 'all':
            range_ = range(1, len(jobs)+1)

        for job in range_:
            j = jobs[str(job)]
            status = j['status']
            rootfile = j['rootfile']
            rootfile = [r[0].encode('ascii') for r in rootfile]

            proc_one_job(rootfile, datatype, label, job, test, batch)
    
    elif job_opt == 'status':
        jobs = get_jobs(stage, datatype, label)
        range_ = option_existsxs(args, 'status')
        if range_ == 'all':
            range_ = range(1, len(jobs)+1)

        sys.stdout.write('| Job | Duration |\n')
        sys.stdout.write('|-----+----------|\n')

        bad_jobs = []

        for job in range_:
            job, processed, selected, duration =  check_job_log(
                'ntp', datatype, label, job)

            if processed == 'N/A':
                bad_jobs.append(job)
                continue

        sys.stdout.write('|-----+----------|\n')
        if bad_jobs:
            bad_jobs = numbers_to_string(bad_jobs)
            sys.stdout.write('Bad jobs: %s\n' % bad_jobs)

    else:
        rootfile = ''
        proc_one_job(rootfile, datatype, label, test=test, grid=grid)

        
def proc_one_job(rootfile, datatype, label, job=None, test=False, 
                 batch=False, grid=False):
    if batch:
        cmd = create_batch_cmd(job)
        bashname = 'ntp_%s.sh' %job
        bashfile = create_bashfile_cmd(cmd, bashname, label, test=test)
        logfile = set_logfile('ntp', datatype, label, 'evt', job)
        jobname = 'ntp%s' %job 
        bsub_jobs(logfile, jobname, bashfile, test)
        return


    if grid:
        crab_cfg = atr.cfg.crab_cfg(label)
        cfg_path, cfg_name = os.path.split(crab_cfg)

        sys.stdout.write('Please create grid jobs with command: \n\n')
        sys.stdout.write('cd %s \n\n' %cfg_path)
        sys.stdout.write('crab -cfg %s -create\n\n' %cfg_name)
        return 


    if job != None: 
        sys.stdout.write('Start processing job %s ...' %job)
    
    time_start = time()
    if label in ['B2KstarMuMu/GEN_1M_v1', 'B2KstarMuMu/RECO_1M_v2']:
        procdir = os.path.join(os.getenv('rel'),
                               'BphAna/B0KstMuMu/python/')
        rfname = 'readFiles_%s.py' %label
        rfname = rfname.replace('/', '_')
        readfile = os.path.join(procdir, rfname)
        list_to_file(rootfile, readfile, listname='readFiles', prefix='rfio:')
        cmd = 'cmsRun B2KstarMuMu_GEN_1M_v1.py'
            
    elif label in ['B2KstarMuMu/RECO_1M_v2.4',
                   'B2KstarMuMu/RECO_10M_v2.1',
                   'B2KstarMuMu/RECO_100M_v1.1']:
        cmssw = 'CMSSW_4_2_8_patch7'
        procdir = os.path.join(os.getenv('afb'), 'rel', cmssw,
                               'src/BphAna/B2KstarMuMu_V00_01_02/test')
        cmd = 'cmsRun %s.py' %label_to_name(label)
        rfname = 'readFiles_%s.py' %label_to_name(label)
        readfile = os.path.join(procdir, rfname)
        list_to_file(rootfile, readfile, listname='readFiles',
                          prefix='root://eoscms//eos/cms')

    elif label in ['B2KstarMuMu/RECO_Brian_v1']:
        procdir = os.path.join(os.getenv('afb'),
                               'rel/CMSSW_4_2_8_patch7/src/BphAna/B0KstMuMu_V00_05_00/python')
        cmd = 'cmsRun %s.py' %label_to_name(label)

    elif label in ['B2KstarMuMu/RECO_Brian_prod_2_1_v1',
                   'B2KstarMuMu/RECO_Brian_prod_2_4_v1',
                   'B2KstarMuMu/RECO_Brian_prod_2_6_v1',
                   'B2KstarMuMu/RECO_Brian_prod_2_7_v1',
                   'B2KstarMuMu/RECO_Brian_prod_2_8_v1',                   
                   'B2KstarMuMu/RECO_Brian_prod_2_2_v1',
                   'B2KstarMuMu/RECO_Brian_prod_2_3_v1',
                   'B0JPsiKstMuMuKPi/RECO_Brian_prod1_v1', 
                   'B0Psi2SKstMuMuKPi/RECO_Brian_prod1_v1'
                   ]:
        procdir = os.path.join(os.getenv('afb'),
        'rel/CMSSW_4_2_8_patch7/src/BphAna/B0KstMuMu_V00_05_02/python')
        cmd = 'cmsRun %s.py' % 'B2KstarMuMu_RECO_Brian_v1'

    elif label in ['B2KstarMuMu/RECO_Brian_prod_2_1_v2']:
        procdir = os.path.join(atr.afbpath,
        'rel/CMSSW_4_2_8_patch7/src/BphAna/B2KstarMuMu_V00_01_02/test')
        cmd = 'cmsRun %s.py' % 'B2KstarMuMu_RECO_Brian_v2'

    elif label in ['B2KstarMuMu/RECO_100M_v1.2',
                   'B2KstarMuMu/RECO_100M_v1.3',
                   'B2KstarMuMu/RECO_100M_v1.4', 
                   'B2KstarMuMu/RECO_100M_v1.5', 
                   ]:
        cmssw = 'CMSSW_4_2_8_patch7'
        procdir = os.path.join(os.getenv('afb'), 'rel', cmssw,
                               'src/BphAna/B2KstarMuMu/test')
        cmd = 'cmsRun %s.py' %label_to_name(label) 
        #rfname = 'readFiles_%s.py' %label_to_name(label)
        #readfile = os.path.join(procdir, rfname)
        #list_to_file(rootfile, readfile, listname='readFiles',
        #                  prefix='root://eoscms//eos/cms')

    elif label in ['B2KstarMuMu/RECO_100M_v1.6',
                   ]:
        cmssw = 'CMSSW_4_2_8_patch7'
        procdir = os.path.join(os.getenv('afb'), 'rel', cmssw,
                               'src/BphAna/B0KstMuMu/python')
        cmd = 'cmsRun %s.py' %label_to_name(label) 

    elif label in ['B2KstarMuMu/RECO_Brian_prod_2_1_v3',
                   'B2KstarMuMu/RECO_Brian_prod_2_2_v3',
                   'B2KstarMuMu/RECO_Brian_prod_2_3_v3',
                   'B2KstarMuMu/RECO_Brian_prod_2_4_v3',
                   'B2KstarMuMu/RECO_Brian_prod_2_6_v3',
                   'B2KstarMuMu/RECO_Brian_prod_2_7_v3',
                   'B2KstarMuMu/RECO_Brian_prod_2_8_v3',
                   'B0JPsiKstMuMuKPi/RECO_Brian_prod1_v3', 
                   'B0Psi2SKstMuMuKPi/RECO_Brian_prod1_v3'
                   ]:
        procdir = os.path.join(os.getenv('afb'),
        'rel/CMSSW_4_2_8_patch7/src/BphAna/B0KstMuMu_V00_06_01/python')
        cmd = 'cmsRun %s.py' % 'B2KstarMuMu_RECO_Brian_v3'

    elif label in ['Run2012B-PromptReco-v1.1',
                   'Run2012B-PromptReco-v1.2',
                   'Run2012B-PromptReco-v1.3',
                   'Run2012B-PromptReco-v1.3.1',
                   ]:
        procdir = os.path.join(atr.afbpath,
        'rel/CMSSW_5_3_2_patch4/src/BphAna/B0KstMuMu/python')
        cmd = 'cmsRun %s.py' % label_to_name(label) 

    elif label in ['Run2012B-PromptReco-v1.4']:
        procdir = os.path.join(atr.afbpath,
        'rel/CMSSW_5_2_5_patch3/src/BphAna/B0KstMuMu/python')
        cmd = 'cmsRun %s.py' % label_to_name(label) 
        
    elif label in ['Run2012B-PromptReco-v2.1']:
        # this is a name error, should be v1.5
        raise NameError(label)
        procdir = os.path.join(atr.afbpath,
        'rel/CMSSW_5_2_5_patch3/src/BphAna/B0KstMuMu_V00_06_07/python')
        cmd = 'cmsRun %s.py' % label_to_name(label) 

    elif label in ['Run2011A-May10ReReco-v1.6']:
        procdir = os.path.join(atr.afbpath,
        'rel/CMSSW_4_2_8_patch7/src/BphAna/B0KstMuMu_V00_06_01/python')
        cmd = 'cmsRun Run2011A-May10ReReco-v1_3.py'

    elif label in ['Run2011A-PromptReco-v4.3',
                   'Run2011A-PromptReco-v5.4',
                   'Run2011A-PromptReco-v6.3',
                   'Run2011B-PromptReco-v1.4',
                   ]:
        procdir = os.path.join(atr.afbpath,
        'rel/CMSSW_4_2_8_patch7/src/BphAna/B0KstMuMu_V00_06_01/python')
        cmd = 'cmsRun Run2011A-PromptReco-v4_3.py'
        
    elif label in ['Run2011A-May10ReReco-v1.7',
                   'Run2011A-May10ReReco-v1.8',
                   'Run2011A-May10ReReco-v1.9',
                   ]:
        procdir = os.path.join(atr.afbpath,
        'rel/CMSSW_4_2_8_patch7/src/BphAna/B0KstMuMu_V00_07_01/python')
        cmd = 'cmsRun Run2011A-May10ReReco-v1_3.py'

    elif label in ['Run2011A-PromptReco-v4.4',
                   'Run2011A-PromptReco-v4.5',
                   'Run2011A-PromptReco-v4.6',
                   'Run2011A-PromptReco-v5.5',
                   'Run2011A-PromptReco-v5.6',
                   'Run2011A-PromptReco-v6.4',
                   'Run2011A-PromptReco-v6.5',
                   'Run2011B-PromptReco-v1.5',
                   'Run2011B-PromptReco-v1.6',
                   ]:        
        procdir = os.path.join(atr.afbpath,
        'rel/CMSSW_4_2_8_patch7/src/BphAna/B0KstMuMu_V00_07_01/python')
        cmd = 'cmsRun Run2011A-PromptReco-v4_3.py'

    elif label in ['B2KstarMuMu/RECO_Brian_prod_2_1_v5',
                   'B2KstarMuMu/RECO_Brian_prod_2_2_v5',
                   'B2KstarMuMu/RECO_Brian_prod_2_3_v5',
                   'B2KstarMuMu/RECO_Brian_prod_2_4_v5',
                   'B2KstarMuMu/RECO_Brian_prod_2_6_v5',
                   'B2KstarMuMu/RECO_Brian_prod_2_7_v5',
                   'B2KstarMuMu/RECO_Brian_prod_2_8_v5',
                   'B0JPsiKstMuMuKPi/RECO_Brian_prod1_v5', 
                   'B0Psi2SKstMuMuKPi/RECO_Brian_prod1_v5'
                   ]:
        procdir = os.path.join(os.getenv('afb'),
        'rel/CMSSW_4_2_8_patch7/src/BphAna/B0KstMuMu_V00_07_01/python')
        cmd = 'cmsRun %s.py' % 'B2KstarMuMu_RECO_Brian_v3'
        
    elif label in ['B2KstarMuMu/RECO_Brian_prod_2_1_v6',
                   'B2KstarMuMu/RECO_Brian_prod_2_2_v6',
                   'B2KstarMuMu/RECO_Brian_prod_2_3_v6',
                   'B2KstarMuMu/RECO_Brian_prod_2_4_v6',
                   'B2KstarMuMu/RECO_Brian_prod_2_6_v6',
                   'B2KstarMuMu/RECO_Brian_prod_2_7_v6',
                   'B2KstarMuMu/RECO_Brian_prod_2_8_v6',
                   'B0JPsiKstMuMuKPi/RECO_Brian_prod1_v6', 
                   'B0Psi2SKstMuMuKPi/RECO_Brian_prod1_v6'
                   ]:
        procdir = os.path.join(os.getenv('afb'),
        'rel/CMSSW_4_2_8_patch7/src/BphAna/B0KstMuMu_V00_08_00/python')
        cmd = 'cmsRun %s.py' % 'B2KstarMuMu_RECO_Brian_v3'
        
    elif label in ['B2KstarMuMu/RECO_Brian_prod_2_1_v7',
                   'B2KstarMuMu/RECO_Brian_prod_2_2_v7',
                   'B2KstarMuMu/RECO_Brian_prod_2_3_v7',
                   'B2KstarMuMu/RECO_Brian_prod_2_4_v7',
                   'B2KstarMuMu/RECO_Brian_prod_2_6_v7',
                   'B2KstarMuMu/RECO_Brian_prod_2_7_v7',
                   'B2KstarMuMu/RECO_Brian_prod_2_8_v7',
                   'B0JPsiKstMuMuKPi/RECO_Brian_prod1_v7.1', 
                   'B0JPsiKstMuMuKPi/RECO_Brian_prod1_v7.2', 
                   'B0Psi2SKstMuMuKPi/RECO_Brian_prod1_v7'
                   ]:
        procdir = os.path.join(os.getenv('afb'),
        'rel/CMSSW_4_2_8_patch7/src/BphAna/B0KstMuMu_V00_08_01/python')
        cmd = 'cmsRun %s.py' % 'B2KstarMuMu_RECO_Brian_v3'

    elif label in ['Run2011A-May10ReReco-v1.10']:
        procdir = os.path.join(atr.afbpath,
        'rel/CMSSW_4_2_8_patch7/src/BphAna/BToKstarMuMu/python')
        cmd = 'cmsRun %s.py' % 'btokstarmumu_Run2011A-May10ReReco-v1_10' 

    elif label in ['Run2011A-PromptReco-v4.10',
                   'Run2011A-PromptReco-v5.10',
                   'Run2011A-PromptReco-v6.10',
                   'Run2011B-PromptReco-v1.10',
                   ]:
        procdir = os.path.join(atr.afbpath,
        'rel/CMSSW_4_2_8_patch7/src/BphAna/BToKstarMuMu/python')
        cmd = 'cmsRun %s.py' % 'btokstarmumu_Run2011A-PromptReco-v4_10' 

    elif label in ['BuToKstarMuMu/7TeV_2p5E3']:
        procdir = os.path.join(atr.afbpath,
        'rel/CMSSW_4_2_8_patch7/src/BphAna/BToKstarMuMu/python')
        cmd = 'cmsRun %s.py' % 'btokstarmumu_Run2011A-PromptReco-v4_10' 

    elif label in ['DoubleMuParked/Run2012B-22Jan2013-v1.1',
                   'ParkingMonitor/Run2012B-PromptReco-v1.1'
                   ]:
        procdir = os.path.join(atr.afbpath,
        'rel/CMSSW_5_3_8_patch3/src/BphAna/BToKstarMuMu/python')
        cmd = 'cmsRun %s.py' % 'btokstarmumu_Run2012A-PromptReco' 
        

    elif label in ['Run2011A-May10ReReco-v1.11']:
        procdir = os.path.join(atr.afbpath,
        'rel/CMSSW_4_2_8_patch7/src/BphAna/BToKstarMuMu_V00_00_02/python')
        cmd = 'cmsRun %s.py' % 'btokstarmumu_Run2011A-May10ReReco-v1_11' 

    elif label in ['Run2011A-PromptReco-v4.11',
                   'Run2011A-PromptReco-v5.11',
                   'Run2011A-PromptReco-v6.11',
                   'Run2011B-PromptReco-v1.11',
                   'Run2011B-PromptReco-v1.11.1',
                   ]:
        procdir = os.path.join(atr.afbpath,
        'rel/CMSSW_4_2_8_patch7/src/BphAna/BToKstarMuMu_V00_00_02/python')
        cmd = 'cmsRun %s.py' % 'btokstarmumu_Run2011A-PromptReco-v4_11' 
 
    elif label in ['BuToKstarMuMu-7TeV-2E7-v1']:
        procdir = os.path.join(
            atr.afbpath, 'rel/CMSSW_4_2_8_patch7/src/BphAna/BToKstarMuMu_v1/python')
        cmd = 'cmsRun %s.py' % 'btokstarmumu_MC' 

    elif label in ['Run2011A-PromptReco-v4.11',
                   'Run2011A-PromptReco-v5.11',
                   'Run2011A-PromptReco-v6.11',
                   #'Run2011B-PromptReco-v1.12',
                   ]:
        procdir = os.path.join(atr.afbpath,
        'rel/CMSSW_4_2_8_patch7/src/BphAna/BToKstarMuMu_V00_00_04/python')
        cmd = 'cmsRun %s.py' % 'btokstarmumu_Run2011A-PromptReco' 

    elif label in ['Run2011B-PromptReco-v1.12',
                   ]:
        procdir = os.path.join(atr.afbpath,
        'rel/CMSSW_4_2_8_patch7/src/BphAna/BToKstarMuMu_Run2011v12/python')
        cmd = 'cmsRun %s.py' % 'btokstarmumu_Run2011A-PromptReco' 

    elif label in ['Run2011B-PromptReco-v1.12',
                   ]:
        procdir = os.path.join(atr.afbpath,
        'rel/CMSSW_4_2_8_patch7/src/BphAna/BToKstarMuMu_Run2011v12/python')
        cmd = 'cmsRun %s.py' % 'btokstarmumu_Run2011A-PromptReco' 

    elif label in ['BuToKstarMuMu-7TeV-2E7-v1.12']:
        procdir = os.path.join(
            atr.afbpath, 'rel/CMSSW_5_3_9_patch3/src/BphAna/BToKstarMuMu_Run2012v1/python')
        cmd = 'cmsRun %s.py' % 'btokstarmumu_MC'  

    elif label in ['BuToKstarMuMu-7TeV-2E7-v1.2']:
        procdir = os.path.join(
            atr.afbpath, 'rel/CMSSW_4_2_8_patch7/src/BphAna/BToKstarMuMu_v2/python')
        cmd = 'cmsRun %s.py' % 'btokstarmumu_MC'  

    elif label in ['BuToKstarMuMu-7TeV-2E7-v1.3']:
        procdir = os.path.join(
            atr.afbpath, 'rel/CMSSW_4_2_8_patch7/src/BphAna/BToKstarMuMu/python')
        cmd = 'cmsRun %s.py' % 'btokstarmumu_MC'  

    elif label in ['Run2011A-May10ReReco-v1_run2011v0']:
        procdir = os.path.join(atr.afbpath,
        'rel/CMSSW_4_2_8_patch7/src/BphAna/BToKstarMuMu_run2011v0/python')
        cmd = 'cmsRun btokstarmumu_Run2011A-May10ReReco-v1.py'

    elif label in ['Run2011A-PromptReco-v4_run2011v0', 
                   'Run2011A-PromptReco-v5_run2011v0', 
                   'Run2011A-PromptReco-v6_run2011v0', 
                   'Run2011B-PromptReco-v1_run2011v0', 
                   ]:
        procdir = os.path.join(atr.afbpath,
        'rel/CMSSW_4_2_8_patch7/src/BphAna/BToKstarMuMu_run2011v0/python')
        cmd = 'cmsRun btokstarmumu_Run2011A-PromptReco.py'
 
    elif label in ['BuToKstarJPsi-7TeV-5E5-v1_run2011v0', 
                   'BuToKstarMuMu-7TeV-2E7-v1_run2011v0', 
                   ]:
        procdir = os.path.join(atr.afbpath,
        'rel/CMSSW_4_2_8_patch7/src/BphAna/BToKstarMuMu_run2011v0/python')
        cmd = 'cmsRun btokstarmumu_MC.py'
 
    elif label in [ 'Run2011B-PromptReco-v1_run2011v0.1' 
                   ]:
        procdir = os.path.join(atr.afbpath,
        'rel/CMSSW_4_2_8_patch7/src/BphAna/BToKstarMuMu_run2011v0/python')
        cmd = 'cmsRun btokstarmumu_Run2011A-PromptReco_v1.py'

    else:
        raise NameError(label)

    proc_cmd(cmd, test, procdir=procdir)
    if test:
        sys.stdout.write('Please test with the above command.\n')
        return 
    dur = duration(time()-time_start)
    sys.stdout.write(' done. \n duration %s \n' %dur)
    sys.stdout.flush()

