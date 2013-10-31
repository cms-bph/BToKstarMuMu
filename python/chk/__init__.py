"""
Module for Checking Data files 

"""

__author__ = "Xin Shi <Xin.Shi@cern.ch>"
__copyright__ = "Copyright (c) Xin Shi"

import sys 
from tls import * 
import atr 

def main(args):
    if args[0] == 'diff':
        return diff(args[1:])

    if args[0] == 'SingleBuToKstarMuMu':
        return SingleBuToKstarMuMu(args[1:])

    datatype = args[0]
    label = args[1]
    jobs_created = atr.dat.jobs_created(datatype, label)
    filepath = atr.dat.get_filepath(datatype, label)

    cmd = 'cmsLs %s' %filepath
    files = get_files_from_ls(cmd)

    goodfiles = get_goodfiles(label)

    jobs_finished = [] 
    for f in goodfiles:
        try: 
            job = f.split('_')[1]
        except IndexError:
            sys.stderr.write('Error in the goodfiles: \n')
            sys.stderr.write(f + '\n') 
            sys.exit()
        jobs_finished.append(int(job))

    jobs_missing = jobs_created.difference(jobs_finished)

    # if '-rmdupl' in args: 
    #     remove_duplicates(filepath, files, goodfiles)
    #     return 
        
    if '-stat' in args: 
        save_stat(filepath, files)
        return 

    sys.stdout.write('Number of missing files : %s \n' % len(jobs_missing))
    remove_duplicates(filepath, files, goodfiles, jobs_missing)
    
    if len(jobs_missing) == 0:
        sys.stdout.write('All jobs are done correctly. \n')
        #remove_duplicates(filepath, files, goodfiles)

    else: 
        sys.stdout.write('Please check and resubmit these jobs: \n')
        crabcmd = 'crab -c crab_%s -resubmit %s' %(
        get_name_from_label(label), ','.join(str(s) for s in jobs_missing))
        sys.stdout.write(crabcmd + '\n')
        #sys.stdout.write('You can use -rmdupl to remove duplicates. \n')


def get_goodfiles(label):
    comname = get_name_from_label(label) 
    cmd = 'find_goodfiles.py -c crab_%s -q' %comname

    if 'run2011v1' in label: 
        procdir = os.path.join(
            atr.afbpath,
            'rel/CMSSW_4_2_8_patch7/src/BphAna/BToKstarMuMu_run2011v1/python')

    elif 'Run2012' in label:
        procdir = os.getcwd()

    else: 
        raise NameError(label)
        
    output = proc_cmd(cmd, procdir=procdir)
    goodfiles = []
    for line in output.split():
        goodfiles.append(line.split('/')[-1])

    return goodfiles


def remove_duplicates(filepath, files, goodfiles, jobs_missing):
    duplicatefiles = set(files).difference(set(goodfiles))
    duplicated_jobs = get_job_number_from_files(duplicatefiles)

    #print duplicated_jobs 
    
    sys.stdout.write('Number of files by cmsLs : %s \n' % len(files))
    sys.stdout.write('Number of goodfiles : %s \n' % len(goodfiles))
    sys.stdout.write('Number of duplicates : %s \n' % len(duplicatefiles))

    if len(duplicatefiles) == 0:
        return

    s = raw_input('Do you want to list the duplicates names? (yes or no) ')
    if s == 'yes': 
        for f in duplicatefiles:
            sys.stdout.write(' %s \n' %f )

    s = raw_input('Do you want to remove all of them? (yes or no) ')
    if s != 'yes':
        return

    sys.stdout.write('Removing duplicates: \n')
    for f in duplicatefiles:
        sys.stdout.write(' %s ...' %f )
        sys.stdout.flush()
        cmd = 'cmsRm %s/%s' %(filepath, f)
        proc_cmd(cmd)
        sys.stdout.write(' OK.\n')

    sys.stdout.write('Total of %s duplicated files removed. \n' %len(duplicatefiles))


def save_stat(filepath, files):
    print filepath, files
    sys.exit()

def get_job_number_from_files(fs):
    jobs = []
    for f in fs:
        job = f.split('_')[1]
        jobs.append(int(job))
    return jobs


def diff(args):
    if args[0] == 'task1':
        f1 = '/Users/xshi/Downloads/B0JPsiKstMuMuKPi_MC_cff.py'
        f2 = '/Users/xshi/Downloads/aaa.py'

        flist1 = get_root_file_names_from_py(f1)
        flist2 = get_root_file_names_from_py(f2)

        diff = set(flist1)-set(flist2)
        print '\n'.join(list(diff))
        
    elif args[0] == 'task2':
        f0 = '/Users/xshi/work/cms/afb/dat/sel/data/Run2011v1.0.1/SingleB2KstarMuMu.root'
        f1 = '/Users/xshi/work/cms/afb/dat/sel/data/Run2011v1.3/SingleB2KstarMuMu.root'
        fout0 = 'output0.txt'
        fout1 = 'output1.txt'

        t0 = root_chain(f0, 'B0SingleCand/B0KstMuMuSingleCandNTuple')
        var0 = 'runN:eventN'
        tree_scan_var_to_file(t0, var0, fout0)

        t1 = root_chain(f1, 'sel')
        var1 = 'run:event'
        tree_scan_var_to_file(t1, var1, fout1)

        list0 = get_run_event_list(fout0)
        list1 = get_run_event_list(fout1)

        s0 = set(list0)
        s1 = set(list1)
        
        u0 = s0.difference(s1)
        u1 = s1.difference(s0)
        print 'Common:  ', len(s0.intersection(s1))
        print 'Unique in 1.0.1: ', len(u0)
        print 'Unique in 1.3: ', len(u1)

        output_list_to_file(u0, 'unique0.txt')
        output_list_to_file(u1, 'unique1.txt')
        
    else:
        raise NameError(args)
    
def get_root_file_names_from_py(filename):
    f = open(filename, 'r')
    fnames = []
    for line in f:
        if '.root' in line:
            fname = line.split('.root')[0]
            fname = fname.split('/')[-1] + '.root'
            fnames.append(fname)

    return fnames


def tree_scan_var_to_file(t, var, f):
    t.GetPlayer().SetScanRedirect(1)
    # You may have to cast t.GetPlayer() to a TTreePlayer*
    t.GetPlayer().SetScanFileName(f)
    t.Scan(var)


def get_run_event_list(filename):
    f = open(filename, 'r')
    run_event_list = []
    line_num = 0
    
    for line in f:
        line_num += 1
        if line_num < 4:
            continue
        
        items=line.split('*')
        run = items[2].strip()
        event = items[3].strip()
        run_event = '%s,%s' %(run, event)
        run_event_list.append(run_event)

    return run_event_list
    
def output_list_to_file(lst, filename):
    f = open(filename, 'w')
    for li in lst:
        f.write(li + '\n')
    f.close()

def SingleBuToKstarMuMu(args):
    datatype = args[0]
    label = args[1]
    comname = 'SingleBuToKstarMuMu'
    #datasets_njobs = atr.ntp.datasets_njobs(label)
    ntp_labels = atr.sel.ntp_labels(label)

    #for dataset, njobs in datasets_njobs:
    for ntp_label in ntp_labels:
        #label = '%s/%s' %(args[1], dataset)
        njobs = atr.sel.njobs(ntp_label)
        for n in range(1, njobs+1):
            logfile = set_file(atr.logpath, ntp_label, comname, '.log.%s' %n)
            parsed = LogFile(logfile)
            sys.stdout.write('%s:%s\t%s\t%s\t%s \n' % (
                ntp_label, n, parsed.processed,   
                parsed.selected, parsed.duration))
            
            
    
