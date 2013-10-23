"""
Module for Generating Monte Carlo 

"""

__author__ = "Xin Shi <Xin.Shi@cern.ch>"
__copyright__ = "Copyright (c) Xin Shi"

import sys 
import atr
from tls import *
from time import time, sleep
import gen 

def main(args):
    datatype = args[0]
    label = args[1]
    test = option_exists(args, '-t')
    pbar = option_exists(args, '-p')
    batch = option_exists(args, '-b')
    job_opt = option_exists(args, 'job')
    grid = option_exists(args, '-grid')

    if label in ['B2KstarMuMu/HLT_1M_v1', 'B2KstarMuMu/RECO_1M_v1',
                 'B2KstarMuMu/HLT_1M_v2', 'B2KstarMuMu/RECO_1M_v2',
                 'B2KstarMuMu/GEN_SIM_10M_v2']:
        stage = 'srm'


    if job_opt == 'create':
        size_max = 100 # do not group jobs
        if label in ['B2KstarMuMu/HLT_1M_v1', 'B2KstarMuMu/HLT_1M_v2']: 
            input_label = 'B2KstarMuMu/GEN_1M_v1'

        elif label in ['B2KstarMuMu/RECO_1M_v1']: 
            input_label = 'B2KstarMuMu/HLT_1M_v1'

        elif label in ['B2KstarMuMu/RECO_1M_v2']: 
            input_label = 'B2KstarMuMu/HLT_1M_v2'

        elif label in ['B2KstarMuMu/GEN_SIM_10M_v2']: 
            input_label = label 

        else:
            raise NameError(label)
        
        dbfile = create_jobs_db(stage, datatype, label, test=test,
                                size_max=size_max, input_label=input_label)
        return

    elif job_opt == 'stager_get' or job_opt == 'stager_qry':
        cmdname = job_opt
        jobs = get_jobs(stage, datatype, label)
        range_ = option_exists(args, job_opt)
        verbose = option_exists(args, 'verbose')

        if range_ == 'all':
            range_ = range(1, len(jobs)+1)
        
        for job in range_:
            nfiles = 0 
            staged = 0
            stagein = 0

            j = jobs[str(job)]
            status = j['status']
            rootfile = j['rootfile']
            rootfiles = [r[0].encode('ascii') for r in rootfile]

            for rootfile in rootfiles:
                nfiles += 1 
                stdout = stager_cmd(cmdname, rootfile, test, verbose=verbose)
                if 'STAGED' in stdout:
                    staged += 1
                if 'STAGIN' in stdout:
                    stagein += 1

            if nfiles == staged:
                status = 'OK.'
            else:
                status =  'Stagein %s' % stagein
                
            sys.stdout.write('Job %s, files %s ... %s \n' %(
                job, nfiles, status))
        
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
        range_ = option_exists(args, 'status')
        if range_ == 'all':
            range_ = range(1, len(jobs)+1)

        sys.stdout.write('| Job | Duration |\n')
        sys.stdout.write('|-----+----------|\n')

        bad_jobs = []

        for job in range_:
            job, processed, selected, duration =  check_job_log(
                'gen', datatype, label, job)

            if processed == 'N/A':
                bad_jobs.append(job)
                continue

        sys.stdout.write('|-----+----------|\n')
        if bad_jobs:
            bad_jobs = numbers_to_string(bad_jobs)
            sys.stdout.write('Bad jobs: %s\n' % bad_jobs)

    else:
        rootfile = ''
        proc_one_job(rootfile, datatype, label, test=test, batch=batch, grid=grid)
        #raise NameError(job_opt)

        
def proc_one_job(rootfile, datatype, label, job=None, test=False, batch=False,
                 grid=False):
    if batch:
        cmd = create_batch_cmd(job)
        cmssw = ''
        if label in ['B2KstarMuMu/HLT_1M_v2', 'B2KstarMuMu/RECO_1M_v2']:
            cmssw = 'CMSSW_4_2_9_HLT1_bphpatch4'
        bashname = 'gen_%s.sh' %job
        bashfile = create_bashfile_cmd(cmd, bashname, label, test=test,
                                       cmssw=cmssw)
        logfile = set_logfile('gen', datatype, label, 'evt', job)
        jobname = 'gen%s' %job 
        bsub_jobs(logfile, jobname, bashfile, test)
        return

    if grid:
        #cfg_name = 'crab_' + label.replace('/', '_') + '.cfg'
        cfg_name = 'crab_%s.cfg' %get_name_from_label(label)
        cfg_path = get_cfg_path(label)
        cfg_file = os.path.join(cfg_path, cfg_name)
        #crab_cfg = atr.crab_cfg(label)
        crab_cfg = atr.cfg.crab_cfg(label)
        sys.stdout.write('Please create grid jobs with command: \n\n')
        sys.stdout.write('cd %s \n\n' %cfg_path)
        sys.stdout.write('crab -cfg %s -create\n\n' %cfg_name)
        return 

        
    if job != None: 
        sys.stdout.write('Start processing job %s ...' %job)
    
    time_start = time()

    if label == 'B2KstarMuMu/HLT_1M_v1':
        rf = rootfile[0]
        suffix = rf.split('GEN')[-1]
        cmd = 'cmsRun B2KstarMuMu_HLT_1M_v1.py %s' % suffix
        procdir = os.getenv('test')
    elif label == 'B2KstarMuMu/RECO_1M_v1':
        rf = rootfile[0]
        suffix = rf.split('7TeV_PYTHIA6')[-1]
        cmd = 'cmsRun B2KstarMuMu_RECO_1M_v1.py %s' % suffix
        procdir = os.getenv('test')
    elif label == 'B2KstarMuMu/HLT_1M_v2':
        rf = rootfile[0]
        suffix = rf.split('GEN')[-1]
        cmd = 'cmsRun B2KstarMuMu_HLT_1M_v2.py %s' % suffix
        procdir = os.path.join(os.getenv('rel'), 'prod')

    elif label == 'B2KstarMuMu/RECO_1M_v2':
        rf = rootfile[0]
        suffix = rf.split('7TeV_PYTHIA6')[-1]
        cmd = 'cmsRun B2KstarMuMu_RECO_1M_v2.py %s' % suffix
        procdir = os.path.join(os.getenv('rel'), 'prod')

    else: 
        config_file_name, procdir = create_config_file(label, test=test) 
        cmd = 'cmsRun %s' % config_file_name


    output = proc_cmd(cmd, procdir=procdir, test=test)
    if test:
        sys.stdout.write('Please test with command: \n\n')
        sys.stdout.write('cd %s && cmsenv && %s\n' %(procdir, cmd))
        return 

    if job != None: 
        dur = duration(time()-time_start)
        sys.stdout.write(' done. \n duration %s \n' %dur)
        sys.stdout.flush()

   

def create_config_file(label, test=False):
    com_name = get_name_from_label(label)
    function = getattr(gen, com_name)
    return function(label)


    nevent = get_number_from_label(label)
    if test:
        nevent = 10

    if label == 'B2KstarMuMu/GEN_SIM_1M_v4':
        procdir = os.path.join(os.getenv('afb'),
                               'rel/CMSSW_4_1_8_patch9/src/prod')
        config_file_name = label.replace('/', '_') + '.py'
        input_file = 'Configuration/GenProduction/python/PYTHIA6_B0dToKstMuMuKPi_7TeV_nofilter_cff.py'
        cmd = 'cmsDriver.py %s  --step GEN,SIM  --beamspot Realistic7TeV2011Collision --conditions START311_V2::All --pileup NoPileUp --datamix NODATAMIXER --eventcontent RAWSIM --datatier GEN-SIM  --fileout=7TeV_PYTHIA6.root --python_filename=%s --number %s --no_exec' %(input_file, config_file_name, nevent)

    elif label in ['B2KstarMuMu/GEN_SIM_1M_v5',
                   'B2KstarMuMu/GEN_SIM_10M_v1',
                   'B2KstarMuMu/GEN_SIM_10M_v2']:

        cmssw = 'CMSSW_4_1_8_patch9'
        procdir = os.path.join(os.getenv('afb'), 'rel', cmssw, 'src/prod')
        config_file_name = label.replace('/', '_') + '.py'
        input_file = config_genprod_file(label, cmssw)
        cmd = 'cmsDriver.py %s  --step GEN,SIM  --beamspot Realistic7TeV2011Collision --conditions START311_V2::All --pileup NoPileUp --datamix NODATAMIXER --eventcontent RAWSIM --datatier GEN-SIM  --fileout=7TeV_PYTHIA6.root --python_filename=%s --number %s --no_exec' %(input_file, config_file_name, nevent)

    elif label in ['B2KstarMuMu/HLT_1M_v5']:
        cmssw = 'CMSSW_4_2_9_HLT1_bphpatch4'
        procdir = os.path.join(os.getenv('afb'), 'rel', cmssw, 'src/prod')
        if not os.access(procdir, os.F_OK):
            sys.stdout.write('Please create the cmssw area first: \n')
            sys.stdout.write('cd $afb/rel && cmsrel %s\n' %cmssw)
            sys.exit()
        config_file_name = label.replace('/', '_') + '.py'
        cmd = 'cmsDriver.py HLT  --step DIGI,L1,DIGI2RAW,HLT:BPh2011  --conditions START42_V14A::All --pileup E7TeV_FlatDist10_2011EarlyData_50ns_PoissonOOT --eventcontent RAWSIM --datatier GEN-SIM-RAW  --fileout=7TeV_PYTHIA6.root --python_filename=%s --number %s --no_exec' %(config_file_name, nevent)

    elif label in ['B2KstarMuMu/HLT_10M_v2']:
        cmssw = 'CMSSW_4_2_9_HLT1_bphpatch4'
        procdir = os.path.join(os.getenv('afb'), 'rel', cmssw, 'src/prod')
        if not os.access(procdir, os.F_OK):
            sys.stdout.write('Please create the cmssw area first: \n')
            sys.stdout.write('cd $afb/rel && cmsrel %s\n' %cmssw)
            sys.exit()
        config_file_name = label.replace('/', '_') + '.py'
        cmd = 'cmsDriver.py HLT  --step DIGI,L1,DIGI2RAW,HLT:BPh2011  --conditions START42_V14A::All --eventcontent RAWSIM --datatier GEN-SIM-RAW  --fileout=7TeV_PYTHIA6.root --python_filename=%s --number %s --no_exec' %(config_file_name, nevent)

    elif label in ['B2KstarMuMu/RECO_10M_v2', 'B2KstarMuMu/RECO_100M_v1']:
        cmssw = 'CMSSW_4_2_9_HLT1_bphpatch4' #'CMSSW_4_4_2_patch8'
        procdir = os.path.join(os.getenv('afb'), 'rel', cmssw, 'src/prod')
        config_file_name = label.replace('/', '_') + '.py'
        cmd = 'cmsDriver.py RECO  --step RAW2DIGI,L1Reco,RECO  --conditions START42_V14A::All --eventcontent RECOSIM --datatier GEN-SIM-RECO  --fileout=7TeV_PYTHIA6.root --python_filename=%s --number %s --no_exec' %(config_file_name, nevent)

    elif label in ['B2KstarMuMu/GEN_SIM_HLT_100M_v1']:
        cmssw = 'CMSSW_4_2_9_HLT1_bphpatch4'
        procdir = os.path.join(os.getenv('afb'), 'rel', cmssw, 'src/prod')
        config_file_name = label.replace('/', '_') + '.py'
        input_file = config_genprod_file(label, cmssw)
        cmd = 'cmsDriver.py %s  --step GEN,SIM,DIGI,L1,DIGI2RAW,HLT:BPh2011  --beamspot Realistic7TeV2011Collision --conditions START42_V14A::All --pileup NoPileUp --datamix NODATAMIXER --eventcontent RAWSIM --datatier GEN-SIM-RAW  --fileout=7TeV_PYTHIA6.root --python_filename=%s --number %s --no_exec' %(input_file, config_file_name, nevent)

    elif label in [ 'B2KstarMuMu/GEN_SIM_10M_v5_3_2_patch4']:
        cmssw = 'CMSSW_5_3_2_patch4'
        procdir = os.path.join(os.getenv('afb'), 'rel', cmssw, 'src/prod')
        config_file_name = label.replace('/', '_') + '.py'
        input_file = config_genprod_file(label, cmssw)
        cmd = 'cmsDriver.py %s  --step GEN,SIM --beamspot Realistic8TeVCollision --conditions START53_V7A::All --pileup 2012_Summer_50ns_PoissonOOTPU --datamix NODATAMIXER --eventcontent RAWSIM --datatier GEN-SIM-RAW  --fileout=8TeV_PYTHIA6.root --python_filename=%s --number %s --no_exec' %(input_file, config_file_name, nevent)

    elif label in ['B2KstarMuMu/GEN_SIM_HLT_10M_v5_3_2_patch4',
                   'B2KstarMuMu/RECO_10M_v5_3_2_patch4']:
        cmssw = 'CMSSW_5_3_2_patch4'
        procdir = os.path.join(os.getenv('afb'), 'rel', cmssw, 'src/prod')
        conditions = 'START53_V7A::All'
        fileout = '8TeV_PYTHIA6.root' 
        config_file_name = label.replace('/', '_') + '.py'

        if label == 'B2KstarMuMu/GEN_SIM_HLT_10M_v5_3_2_patch4' :
            input_file = config_genprod_file(label, cmssw)
            cmd = 'cmsDriver.py %s  --step GEN,SIM,DIGI,L1,DIGI2RAW,HLT:7E33v2 --beamspot Realistic8TeVCollision --conditions %s --datamix NODATAMIXER --eventcontent RAWSIM --datatier GEN-SIM-RAW  --fileout=%s --python_filename=%s --number %s --no_exec' %(
                input_file, conditions, fileout, config_file_name, nevent)
        else:
            cmd = 'cmsDriver.py RECO  --step RAW2DIGI,L1Reco,RECO  --conditions %s --eventcontent RECOSIM --datatier GEN-SIM-RECO  --fileout=%s --python_filename=%s --number %s --no_exec' %(
                conditions, fileout, config_file_name, nevent)

    else:
        raise NameError(label)

    config_file = check_and_join(procdir, config_file_name)
    stdout = proc_cmd(cmd, procdir=procdir, test=False)
        
    sys.stdout.write(stdout)    
    check_update_status(config_file)
    return config_file_name, procdir 


def config_genprod_file(label, cmssw):
    if label in ['B2KstarMuMu/GEN_SIM_HLT_10M_v1',
                 'B2KstarMuMu/GEN_SIM_10M_v1',
                 'B2KstarMuMu/GEN_SIM_10M_v2',
                 'B2KstarMuMu/GEN_SIM_HLT_100M_v1']:
        cfg_name = 'PYTHIA6_B0dToKstMuMuKPi_7TeV_b0filter_cff.py'
        src = os.path.join(os.getenv('afb'), 'src/py/atr/decs', cfg_name)
        dstpath = os.path.join(os.getenv('afb'), 'rel', cmssw, 'src',
                               'Configuration/GenProduction/python')
        dst = check_and_join(dstpath, cfg_name)
        check_and_copy(src, dst)
        config_file = 'Configuration/GenProduction/python/%s' % cfg_name

    elif label in ['B2KstarMuMu/GEN_SIM_HLT_10M_v5_3_2_patch4']:
        cfg_name = 'PYTHIA6_B0dToKstMuMuKPi_8TeV_b0filter_cff.py'
        src = os.path.join(os.getenv('afb'), 'src/py/atr/decs', cfg_name)
        dstpath = os.path.join(os.getenv('afb'), 'rel', cmssw, 'src',
                               'Configuration/GenProduction/python')
        dst = check_and_join(dstpath, cfg_name)
        check_and_copy(src, dst)
        scram_build(cmssw)
        config_file = 'Configuration/GenProduction/python/%s' % cfg_name

    elif label in ['BuToKstarMuMu/GEN_SIM_1M_v1',
                   'BuToKstarMuMu/GEN_SIM_10M_v1.1' ]:
         config_file = 'Configuration/GenProduction/python/SevenTeV/\
PYTHIA6_BuToKstarMuMu_EtaPtFilter_TuneZ2_7TeV_cff.py'

    else:
        raise NameError(label)

    return config_file 


def get_cfg_path(label):
    if label in ['B2KstarMuMu/GEN_SIM_1M_v4',
                 'B2KstarMuMu/GEN_SIM_1M_v5',
                 'B2KstarMuMu/GEN_SIM_10M_v1', 
                 'B2KstarMuMu/GEN_SIM_10M_v2']:
        cfg_path = os.path.join(os.getenv('afb'),
                                'rel/CMSSW_4_1_8_patch9/src/prod/')
    elif label in ['B2KstarMuMu/HLT_1M_v5',
                   'B2KstarMuMu/HLT_10M_v2',
                   'B2KstarMuMu/RECO_10M_v2',
                   'B2KstarMuMu/GEN_SIM_HLT_100M_v1', 
                   'B2KstarMuMu/RECO_100M_v1']:
        cfg_path = os.path.join(os.getenv('afb'),
                                'rel/CMSSW_4_2_9_HLT1_bphpatch4/src/prod/')
    elif label in ['B2KstarMuMu/GEN_SIM_HLT_10M_v5_3_2_patch4',
                   'B2KstarMuMu/RECO_10M_v5_3_2_patch4']:
        cfg_path = os.path.join(atr.afbpath, 'rel/CMSSW_5_3_2_patch4/src/prod/')

    elif label in ['BuToKstarMuMu/GEN_SIM_1M_v1',
                   'BuToKstarMuMu/HLT_1M_v1',
                   'BuToKstarMuMu/RECO_1M_v1',
                   'BuToKstarMuMu/GEN_SIM_10M_v1.1',
                   'BuToKstarMuMu/HLT_10M_v1',
                   'BuToKstarMuMu/RECO_10M_v1',
                   ]:
        cfg_path = os.path.join(atr.afbpath, 'rel/CMSSW_4_2_9_HLT1_bphpatch4/src/')

    else:
        raise NameError(label)

    return cfg_path


def BuToKstarMuMu_GEN_SIM_1M_v1(label):
    nevent = get_number_from_label(label)
    cmssw = 'CMSSW_4_2_9_HLT1_bphpatch4'
    procdir = os.path.join(atr.afbpath, 'rel', cmssw, 'src')
    input_file = config_genprod_file(label, cmssw)
    #config_file_name = 'BuToKstarMuMu_GEN_SIM_1M_v1.py'
    config_file_name = '%s.py' %get_name_from_label(label)

    cmd = 'cmsDriver.py %s  --step GEN,SIM  --beamspot Realistic7TeV2011Collision --conditions START311_V2::All --pileup NoPileUp --datamix NODATAMIXER --eventcontent RAWSIM --datatier GEN-SIM  --fileout=BuToKstarMuMu_7TeV_PYTHIA6.root --python_filename=%s --number %s --no_exec' %(input_file, config_file_name, nevent)

    config_file = check_and_join(procdir, config_file_name)
    stdout = proc_cmd(cmd, procdir=procdir, test=False)
        
    sys.stdout.write(stdout)    
    check_update_status(config_file)
    return config_file_name, procdir 
    

def BuToKstarMuMu_HLT_1M_v1(label):
    cmssw = 'CMSSW_4_2_9_HLT1_bphpatch4'
    procdir = os.path.join(atr.afbpath, 'rel', cmssw, 'src')
    nevent = get_number_from_label(label)
    config_file_name = '%s.py' %get_name_from_label(label)
    config_file = check_and_join(procdir, config_file_name)
 
    cmd = 'cmsDriver.py HLT  --step DIGI,L1,DIGI2RAW,HLT:BPh2011  --conditions START42_V14A::All --eventcontent RAWSIM --datatier GEN-SIM-RAW  --fileout=BuToKstarMuMu_7TeV_PYTHIA6.root --python_filename=%s --number %s --no_exec' %(config_file_name, nevent)

    stdout = proc_cmd(cmd, procdir=procdir, test=False)
        
    sys.stdout.write(stdout)    
    check_update_status(config_file)
    return config_file_name, procdir 


def BuToKstarMuMu_RECO_1M_v1(label):
    cmssw = 'CMSSW_4_2_9_HLT1_bphpatch4'
    procdir = os.path.join(atr.afbpath, 'rel', cmssw, 'src')
    nevent = get_number_from_label(label)
    config_file_name = '%s.py' %get_name_from_label(label)
    config_file = check_and_join(procdir, config_file_name)
    
    cmd = 'cmsDriver.py RECO  --step RAW2DIGI,L1Reco,RECO  --conditions START42_V14A::All --eventcontent RECOSIM --datatier GEN-SIM-RECO  --fileout=BuToKstarMuMu_7TeV_PYTHIA6.root --python_filename=%s --number %s --no_exec' %(config_file_name, nevent)
    
    stdout = proc_cmd(cmd, procdir=procdir, test=False)
        
    sys.stdout.write(stdout)    
    check_update_status(config_file)
    return config_file_name, procdir 
    

def BuToKstarMuMu_GEN_SIM_10M_v1_1(label):
    return BuToKstarMuMu_GEN_SIM_1M_v1(label)
    
def BuToKstarMuMu_HLT_10M_v1(label):
    return BuToKstarMuMu_HLT_1M_v1(label)
    
def BuToKstarMuMu_RECO_10M_v1(label):
    return BuToKstarMuMu_RECO_1M_v1(label)
