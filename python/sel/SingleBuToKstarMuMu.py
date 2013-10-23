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
    cut = args[2]

    ntp_labels = atr.sel.ntp_labels(datatype, label)

    for ntp_label in ntp_labels:
        proc_ntuple(args, ntp_label, cut)
	    

def proc_ntuple(args, label, cut): 
    datatype = args[0]

    test = option_exists(args, '-t')
    batch = option_exists(args, '-b')

    inname = atr.ntp.rootname(datatype, label, batch) 
    comname = 'SingleBuToKstarMuMu'

    outname = comname    
    if batch :
        outname = comname+'_${LSB_JOBINDEX}'

    inpath = atr.ntp.rootpath(datatype, label)
    infile = os.path.join(inpath, inname+'.root')
    
    outpath = os.path.join(atr.datpath, 'sel', datatype, label)
    outfile = set_file(outpath, cut, outname, '.root', test=test)

    sel_datatype = datatype 
    if datatype == 'mc' and 'BuToKstarJPsi' in label:
        sel_datatype =  'BuToKstarJPsi' 
    if datatype == 'mc' and 'BuToKstarMuMu' in label:
        sel_datatype =  'BuToKstarMuMu' 

    procdir = atr.sel.procdir(label)
    cmd = './sel %s %s %s %s' %(sel_datatype, cut, infile, outfile)

    if option_exists(args, '-n') : 
        nentries = get_option(args, '-n') 
        cmd = '%s -n %s' %(cmd, nentries)
    else:
        if test:
            cmd = '%s -n 1000' % cmd

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
        njobs = atr.ntp.num_rootfiles[label]
        jobname = 'selb[1-%s]' %njobs  
        if option_exists(args, '-J'):
            jobname = get_option(args, '-J')

        bsub_jobs(logfile, jobname, bashfile,
                  test=option_exists(args, '-btest'), queue='1nh')
        return

    output = proc_cmd(cmd, procdir=procdir, test=option_exists(args, '-ctest'))
    if output is not None: 
        sys.stdout.write(output)

    if option_exists(args, '-ctest'):
        return 

    dur = duration(time()-time_start)
    sys.stdout.write(' \nDone in %s. \n' % dur)
    sys.stdout.flush()
