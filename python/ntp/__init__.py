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




    if job != None: 
        sys.stdout.write('Start processing job %s ...' %job)
    
    time_start = time()

    procdir = os.getcwd()

    if 'BuToKstarJPsi-7TeV-5E5-v1' in label or 'BuToKstarMuMu' in label: 
        pset = 'btokstarmumu_MC.py'
        
    elif 'Run2011A-May10ReReco' in label: 
        pset = 'btokstarmumu_Run2011A-May10ReReco.py'

    elif 'Run2011B' in label and 'PromptReco' in label: 
        pset = 'btokstarmumu_Run2011A-PromptReco.py'

    elif 'Run2012' in label: 
        pset = 'btokstarmumu_Run2012.py'

    else:
        raise NameError(label)

    cmd = 'cmsRun %s' % pset 


    if test:
        proc_cmd(cmd, test, procdir=procdir)
        sys.stdout.write('Please test with the above command.\n')
        return 
 
    if grid:
        crab_cfg = atr.cfg.crab_cfg(label, pset)
        cfg_path, cfg_name = os.path.split(crab_cfg)

        sys.stdout.write('Please create grid jobs with command: \n\n')
        #sys.stdout.write('cd %s \n\n' %cfg_path)
        sys.stdout.write('crab -cfg %s -create\n\n' %cfg_name)
        return 

    dur = duration(time()-time_start)
    sys.stdout.write(' done. \n duration %s \n' %dur)
    sys.stdout.flush()

