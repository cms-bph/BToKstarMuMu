"""
Module for SingleBuToKstarMuMu Selection

"""

__author__ = "Xin Shi <Xin.Shi@cern.ch>"
__copyright__ = "Copyright (c) Xin Shi"

import sys
import atr
from time import time
from tls import *

def main(args):
    datatype = args[0]
    label = args[1]

    #ntp_labels = [label]
    ntp_labels = atr.sel.ntp_labels(label)

    for ntp_label in ntp_labels:
        proc_dataset(args, ntp_label)
	    

def proc_dataset(args, label): 
    datatype = args[0]

    test = option_exists(args, '-t')
    batch = option_exists(args, '-b')

    #inname = 'BToKstarMuMu*'
    inname = 'BToKstarMuMu_merged_1'
    comname = 'SingleBuToKstarMuMu'

    outname = comname    

    # if datatype == 'data' and test:
    #     inname = 'BToKstarMuMu_1'

    #if batch and option_exists(args, '-J'):
    if batch :
        inname = 'BToKstarMuMu_merged_${LSB_JOBINDEX}'
        outname = comname+'_${LSB_JOBINDEX}'
        
    # eosbase = 'root://eoscms//eos/cms/store/user/xshi/'
    # inpath = os.path.join(eosbase, 'dat/ntp/data', label)

    inpath = atr.ntp.rootpath(datatype, label)
    
    infile = os.path.join(inpath, inname+'.root')
    
    outpath = os.path.join(atr.datpath, 'sel', datatype)
    outfile = set_file(outpath, label, outname, '.root', test=test)

    binname = './SingleBuToKstarMuMuSelector'
    #procdir = os.path.join(os.environ['HOME'], 'work/cms/afb/src/cc')
    procdir = atr.sel.procdir(label)

    selector_label = label.split('/')[0]

    cmd = '%s %s %s %s' %(binname, selector_label, infile, outfile)

    if option_exists(args, '-n'): 
        nentries = get_option(args, '-n') 
        cmd = '%s -n %s' %(cmd, nentries)

    if option_exists(args, '-j'): 
        nworkers = get_option(args, '-j')
        cmd = '%s -j %s' %(cmd, nworkers)

    time_start = time()

    if batch:
        afb = atr.get_afb_from_label(label)
        pre = 'setafb %s\n\ncd %s' % (afb, procdir) 
        bashfile = set_file(atr.bashpath, label, comname, '.sh', test=test)
	update_bashfile_cmd(bashfile, cmd, pre=pre, test=test)
        logfile = set_file(atr.logpath, label, comname, '.log', test=test)

        #jobname = 'selbu'
        njobs = atr.sel.njobs(label)
        jobname = 'selb[1-%s]' %njobs  
        if option_exists(args, '-J'):
            jobname = get_option(args, '-J')

        bsub_jobs(logfile, jobname, bashfile,
                  test=option_exists(args, '-btest'), queue='8nh')
        return

    output = proc_cmd(cmd, procdir=procdir, test=option_exists(args, '-ctest'))
    if output is not None: 
        sys.stdout.write(output)

    if option_exists(args, '-ctest'):
        return 

    dur = duration(time()-time_start)
    sys.stdout.write(' \nDone in %s. \n' % dur)
    sys.stdout.flush()
