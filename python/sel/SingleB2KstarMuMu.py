"""
Module for SingleB2KstarMuMu Selection

"""

__author__ = "Xin Shi <Xin.Shi@cern.ch>"
__copyright__ = "Copyright (c) Xin Shi"

import sys
import atr
from time import time
from shutil import move
from tls import *
from ROOT import gROOT, kTRUE


def main(args):
    comname = 'SingleB2KstarMuMu'

    gROOT.SetBatch(kTRUE)
    datatype = args[0]
    label = args[1]
    test = get_options(args, 'test')
    batch = get_options(args, 'batch')

    if batch:
        cmd = create_batch_cmd()
        bashname = 'sel_%s.sh' %comname
        bashpath = os.path.join(atr.afbpath, 'src', 'sh', label)
        afb = atr.get_afb_from_label(label)
        bashfile = create_bashfile_cmd(cmd, bashpath, bashname, afb, test=test)

        logpath = os.path.join(atr.afbpath, 'log', 'sel', datatype, label)
        logfile = set_logfile(logpath, comname, test=test)
        jobname = 'selb2kmm' 
        bsub_jobs(logfile, jobname, bashfile, test)
        return


    ntp_label = get_label_by_version(label, 'x')
    infile = atr.ntp.rootfile(datatype, ntp_label)

    outfile = atr.sel.rootfile(datatype, label, comname, create=True, test=test)
    
    ccname = 'SingleB2KstarMuMuSelector.cc+'
    if '-fcompile' in args:
        ccname += '+'

    ccpath = os.path.join(os.environ['HOME'], 'bak/src/cms/afb/cc')
    ccfile = os.path.join(ccpath, ccname)

    t = root_chain(infile, 'B0Cand/B0KstMuMuNTuple')
    
    if test:
        nentries = 10000
    else:
        nentries = t.GetEntries()


    option = 'outfile=%s' % outfile
    if datatype == 'data':
        option += 'checkTrigger1'
        
    if '-d' in args:
        option += 'debug'
        nentries = 2 

    firstentry = 0
    if '-firstentry' in args:
        firstentry = int(args[args.index('-firstentry') + 1])

    time_start = time()

    t.Process(ccfile, option, nentries, firstentry)

    dur = duration(time()-time_start)
    sys.stdout.write(' \nDone in %s. \n' % dur)

    if test:
        estimate = duration(t.GetEntries()/float(nentries)*(time()-time_start))
        sys.stdout.write('Estimate for the total %s : %s\n' % (t.GetEntries(), estimate))
    
    sys.stdout.flush()

    
