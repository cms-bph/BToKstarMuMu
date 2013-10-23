"""
Module for Get variables from the selected events

"""

__author__ = "Xin Shi <Xin.Shi@cern.ch>"
__copyright__ = "Copyright (c) Xin Shi"

import sys 
from tls import * 
from ROOT import (TFile)#, TTree
                  
def main(args):
    if args[0] == 'run_event':
        run_event(args[1:])
    elif args[0] == 'diff_run_event':
        diff_run_event(args[1:])
    else:
        raise NameError(args)

def run_event(args):
    datatype = args[0]
    label = args[1]
    test = get_options(args, 'test')
    pbar = get_options(args, 'pbar')
    job_opt = get_options(args, 'job')

    selfile = get_selfile(datatype, label, job_opt)
    if job_opt == 'stager_get' or job_opt == 'stager_qry':
        cmdname = job_opt
        for f in selfile:
            stager_cmd(cmdname, f, test)

    chain = get_chain(selfile, 'sel')
    entries = chain.GetEntries()
    if test:
        entries = 300

    if pbar:
        pb = get_progressbar(maxval=entries)

    ntot = 0

    logfile = set_logfile('get', datatype, label, 'run_event', test=test)
    l = UserFile()

    for i in xrange(entries):
        ntot += 1 
        if pbar:
            pb.update(i+1)
            
        chain.LoadTree(i)
        chain.GetEntry(i)

        l.append('%s %s\n' %(chain.run, chain.event))
        if test and ntot > 10:
            break

    if pbar:
        pb.finish()

    l.output(logfile)


def diff_run_event(args):
    datatype = args[0]
    label = args[1]
    logfile = set_logfile('get', datatype, label, 'run_event')

    if label == '5ifbv2.3.1':
        l2 = set_logfile('get', datatype, label, 'RunEvt_Linlin')
        l3 = set_logfile('get', datatype, label, 'run_event_xin_unique')
        
    events_xin = get_events_set(logfile)
    events_lin = get_events_set(l2)

    print 'Xin: ', len(events_xin)
    print 'Lin: ', len(events_lin)

    events_union = events_xin | events_lin
    events_inter = events_xin & events_lin
    events_xin_unique = events_xin - events_inter

    print 'Xin Unique: ', len(events_xin_unique)
    
    output_set(events_xin_unique, l3)
    
    
    

    
