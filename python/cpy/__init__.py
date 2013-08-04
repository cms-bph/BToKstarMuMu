"""
Module for Copy Data files 

"""

__author__ = "Xin Shi <Xin.Shi@cern.ch>"
__copyright__ = "Copyright (c) Xin Shi"

import sys 
from tls import * 
import atr

def main(args):
    stage = args[0]
    datatype = args[1]
    label = args[2]
    test = get_options(args, 'test')
    pbar = get_options(args, 'pbar')
    #verbose = get_options(args, 'verbose')
    verbose = False
    if '-v' in args:
        verbose = True 

    if stage == 'sel':
        return cpy_sel(args[1:])
        
    job_opt = get_options(args, 'job')
    batch = get_options(args, 'batch')


    if job_opt == 'create':
        dbfile = create_jobs_db(stage, datatype, label, test=test)
        return


    elif job_opt == 'stager_get' or job_opt == 'stager_qry':
        cmdname = job_opt
        jobs = get_jobs(stage, datatype, label)
        range_ = get_options(args, job_opt)
        verbose = get_options(args, 'verbose')

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
                if 'STAGEIN' in stdout:
                    stagein += 1

            if nfiles == staged:
                status = 'OK.'
            else:
                status =  'STAGED: %s, STAGEIN: %s' % (staged, stagein)
                
            sys.stdout.write('Job %s, files %s ... %s \n' %(
                job, nfiles, status))
        
    
    elif job_opt == 'submit':
        jobs = get_jobs(stage, datatype, label)
        range_ = get_options(args, 'submit')#, start=1, stop=len(jobs)+1)
        
        for job in range_:
            j = jobs[str(job)]
            rootfile = j['rootfile']
            #rootfile = [r.encode('ascii') for r in rootfile]
            sys.stdout.write('Process job %s ... \n' % job)
            sys.stdout.flush()
            proc_one_job(stage, rootfile, datatype, label, job, test, batch, pbar)

    else:
        raise NameError(job_opt)


def cpy_sel(args): 
    name = args[0]
    datatype = args[1]
    label = args[2]
    test = get_options(args, 'test')

    srcfile = atr.sel.rootfile(datatype, label, name, remote=True)
    dstfile = atr.sel.rootfile(datatype, label, name, create=True)

    cmd = 'scp %s %s' %(srcfile, dstfile)
    output = proc_cmd(cmd, test=test)
    
    #print output
    
    #sys.stdout.write(output)

    return

    #srcfile = get_selfile(datatype, label)
    # if label == '5ifbv2.3.1':
    #     dst = os.path.join(atr.datpath, 'sel', datatype, label)
    # else:
    #     raise NameError(label)
    if pbar:
        pb = get_progressbar(maxval=len(srcfile))

    i = 0 
    for src in srcfile:
        i += 1
        cmd = 'rfcp %s %s' %(src, dst) 
        output = proc_cmd(cmd, test=test)
        if verbose: 
            sys.stdout.write(output)
            
        if pbar:
            pb.update(i)

    if pbar:
        pb.finish()
        

def proc_one_job(stage, rootfile, datatype, label, job, test, batch, pbar):
    if batch:
        cmd = create_batch_cmd(job)
        bashname = 'cpy_%s_%s.sh' %(stage, job )
        bashfile = create_bashfile_cmd(cmd, bashname, test=test)
        logfile = set_logfile('cpy', datatype, label, stage, job)
        jobname = 'cpy%s' % job
        bsub_jobs(logfile, jobname, bashfile, test)
        return

    if pbar:
        pb = get_progressbar(maxval=len(rootfile))

    i = 0

    for rf in rootfile:
        i += 1
        #print rf
        f = rf[0].encode('ascii')
        print i, ':  ', f

        if i < 20:
            continue 
        
        if label == 'B2KstarMuMu/RECO_1M_v2' :
            copy_to_cern_t2(f, test)
        else:
            raise NameError(label)

        # cf = create_castor_file(stage, f, datatype, label)
        # dst_size = get_file_size(cf)

        # if dst_size != s:
        #     src = 'srm://srm.ihep.ac.cn:8443/srm/managerv2?SFN=%s' %f
        #     dst = 'srm://srm-cms.cern.ch:8443/srm/managerv2?SFN=%s' %cf 
        #     cmd = 'lcg-cp %s %s' %(src, dst) 
        #     output = proc_cmd(cmd, test=test)
        #     new_dst_size = get_file_size(cf)
        #     if new_dst_size != s:
        #         sys.stdout.write('problem file: %s.\n' %f)
        #         print output
        #         print cmd 
        #         sys.exit()
        if pbar:
            pb.update(i)

    if pbar:
        pb.finish()


def create_castor_file(stage, f, datatype, label):
    if label in ['5ifbv2.3']:
        subdir, name = shorten_subdir_name(f)
        p = os.path.join(atr.castor_datpath, stage, datatype, label, subdir)
    else:
        raise NameError(label)

    cf = check_and_join(p, name)
    return cf 


def shorten_subdir_name(f):
    path, name = os.path.split(f)
    subdir = path.split('/')[-2]
    return subdir, name


def copy_to_cern_t2(f, test):
    src = f.replace('/castor/cern.ch', '')
    dst = src.replace('/user/x', '/store/user')
    src_size = get_file_size(f)
    #print 'src size: ', src_size
    dst_size = get_file_size(dst, 'cmsLs')

    #print 'dst size: ', dst_size
    
    if src_size == dst_size:
        return

    #print 'here'
    cmd = 'cmsStage %s %s' %(src, dst)
    sys.stdout.write('%s ...' %cmd)
    sys.stdout.flush()
    
    stdout = proc_cmd(cmd, test=test)

    if stdout and 'does not exist' in stdout:
        dst_path = os.path.split(dst)[0]
        sys.stdout.write('creating dir %s ... ' % dst_path)
        proc_cmd('cmsMkdir %s' % dst_path)
        sys.stdout.write(' OK.\n')
        
        stdout = proc_cmd(cmd, test=test)
        if 'does not exist' in stdout:
            raise ValueError(stdout)

    sys.stdout.write(' OK.\n')
