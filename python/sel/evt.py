"""
Module for A_FB EVT Selection

"""

__author__ = "Xin Shi <Xin.Shi@cern.ch>"
__copyright__ = "Copyright (c) Xin Shi"


import sys 
from time import time
from array import array
from tls import *
from atr.cuts import select_b0s

from ROOT import (gROOT, TObject, TClonesArray, AddressOf, TLorentzVector,
                  TFile, TTree, kTRUE)

def main(args):
    #gROOT.SetBatch(kTRUE)
    datatype = args[0]
    label = args[1]
    test = get_options(args, 'test')
    pbar = get_options(args, 'pbar')
    batch = get_options(args, 'batch')
    job_opt = get_options(args, 'job')

    main_label = label.split('/')[0]
    if main_label in ['5ifbv2.3']:
        stage = 'skm'        
    elif main_label in ['5ifbv2.3.1', '5ifbv2.3.2', '5ifbv2.8.1']:
        stage = 'ntp'
    else:
        raise NameError(label)


    if job_opt == 'create':
        dbfile = create_jobs_db(stage, datatype, label)
        return

    if job_opt == 'update':
        dbfile = update_jobs_db(stage, datatype, label)
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
            if label in ['5ifbv2.2']:
                rootfiles = [r.encode('ascii') for r in rootfile]
            else:
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
        range_ = get_options(args, 'submit')
        if range_ == 'all':
            range_ = range(1, len(jobs)+1)

        for job in range_:
            j = jobs[str(job)]
            status = j['status']
            rootfile = j['rootfile']
            if label in ['5ifbv2.2']:
                rootfile = [r.encode('ascii') for r in rootfile]
            else:
                rootfile = [r[0].encode('ascii') for r in rootfile]
            proc_one_job(rootfile, datatype, label, job, test, batch, pbar)

    elif job_opt == 'status':
        jobs = get_jobs(stage, datatype, label)
        range_ = get_options(args, 'status')
        if range_ == 'all':
            range_ = range(1, len(jobs)+1)

        sys.stdout.write('| Job | Processed | Selected | Duration |\n')
        sys.stdout.write('|-----+-----------+----------+----------|\n')

        total_processed = 0
        total_selected = 0 

        bad_jobs = []

        for job in range_:
            job, processed, selected, duration =  check_job_log(
                'sel', datatype, label, job)

            if processed == 'N/A':
                bad_jobs.append(job)
                continue
            total_processed += int(processed)
            total_selected += int(selected)

        sys.stdout.write('|-----+-----------+---------+----------|\n')
        sys.stdout.write('Total: processed %s , selected %s, (%s) \n' %(
        total_processed, total_selected, total_selected/float(total_processed)))
        if bad_jobs:
            bad_jobs = numbers_to_string(bad_jobs)
            sys.stdout.write('Bad jobs: %s\n' % bad_jobs)

    else:
        raise NameError(job_opt)

        
def proc_one_job(rootfile, datatype, label, job=None, test=False, batch=False,
                 pbar=False):
    if batch:
        cmd = create_batch_cmd(job)
        bashname = 'sel_evt_%s.sh' %job
        bashfile = create_bashfile_cmd(cmd, bashname, label, test=test)
        logfile = set_logfile('sel', datatype, label, 'evt', job)
        jobname = 'sel%s' % job
        bsub_jobs(logfile, jobname, bashfile, test)
        return

    time_start = time()

    if label in ['5ifbv2.3']:
        chain = get_chain(rootfile, 'skm')
        ntot, nsel, selfile = output_evt_mass_b0_with_skm(
            chain, datatype, label, test=test, pbar=pbar)

    elif label in ['5ifbv2.8.1']:
        chain = get_chain(rootfile, 'B0Cand/B0KstMuMuNTuple')
        ntot, nsel, selfile = output_single_cand_b0_kstmumu(
            chain, datatype, label, job, test=test, pbar=pbar)

    else:
        chain = get_chain(rootfile, 'data')
        ntot, nsel, selfile = output_evt_mass_b0(
            chain, datatype, label, job, test=test, pbar=pbar)


    dur = duration(time()-time_start)
    sys.stdout.write(' processed %s \n selected %s \n duration %s \n saved as %s\n'
                     % (ntot, nsel, dur, selfile))
    sys.stdout.flush()


def output_evt_mass_b0(chain, datatype, label, job, test=False, pbar=False):
    selfile = set_selfile(datatype, label, 'tree', job=job, test=test)

    f = TFile.Open(selfile, 'RECREATE')
    t = TTree('sel', 'sel')

    s = SelTree(t)

    cuts = select_b0s(label)

    pass_trigger = cuts.pass_trigger  
    pass_b0s = cuts.pass_b0s
    pass_jpsimass = cuts.pass_jpsimass
    pass_psi2smass = cuts.pass_psi2smass

    nsel = 0 
    ntot = 0

    xp4_ = TClonesArray('TLorentzVector')
    chain.SetBranchAddress('xP4', AddressOf(xp4_))

    kstarp4_ = TClonesArray('TLorentzVector')
    chain.SetBranchAddress('KstarP4', AddressOf(kstarp4_))

    entries = chain.GetEntries()
    if test:
        entries = 3000

    if pbar:
        pb = get_progressbar(maxval=entries)

    for i in xrange(entries):
        ntot += 1
        if pbar:
            pb.update(i+1)
            
        ientry = chain.LoadTree(i)
        if ientry < 0:
            break

        nb = chain.GetEntry(i)
        if nb <= 0:
            continue

        if not pass_trigger(chain):
            continue

        chain.LoadTree(i)
        chain.GetEntry(i)

        b0s = pass_b0s(chain)
        if b0s == []:
            continue

        s.run[0] = chain.runNb
        s.event[0] = chain.eventNb

        #nb0 = 0
        #nj = 0
        #np = 0 
        #nout = 0 
        for i in b0s:
            #nb0 += 1 
            b0p4 = chain.xP4[i]
            kstarp4 = chain.KstarP4[i]
            oniap4 = b0p4 - kstarp4

            s.b0mass[0] = b0p4.M()
            s.oniamass[0] = oniap4.M()

            s.b0massinjpsi[0] = 0 
            s.b0massinpsi2s[0] = 0
            s.b0massout[0] = 0
            
            if pass_jpsimass(oniap4.M()):
                s.b0massinjpsi[0] = b0p4.M()
                #nj += 1 
            elif pass_psi2smass(oniap4.M()):
                s.b0massinpsi2s[0] = b0p4.M()
                #np += 1
            else:
                s.b0massout[0] = b0p4.M()
                #nout += 1

        #print '%s + %s + %s = %s' %(nj, np, nout, nb0)

        t.Fill()
        nsel += 1
        if nsel > 100 and test:
            break

    f.Write()
    f.Close()

    if pbar:
        pb.finish()

    return ntot, nsel, selfile 



def output_evt_mass_b0_with_skm(chain, datatype, label, test=False, pbar=False):
    selfile = set_selfile(datatype, label, 'tree', test=test)
    f = TFile(selfile, 'RECREATE')
    t = TTree('sel', 'sel')

    s = SelTree(t)

    nsel = 0 
    ntot = 0

    entries = chain.GetEntries()
    if test:
        entries = 3000

    if pbar:
        pb = get_progressbar(maxval=entries)

    for i in xrange(entries):
        ntot += 1
        if pbar:
            pb.update(i+1)
            
        chain.LoadTree(i)
        chain.GetEntry(i)

        b0s = select_b0s_with_skm(chain, label)
        if b0s == []:
            continue

        s.run[0] = chain.run
        s.event[0] = chain.event

        for i in b0s:
            s.b0mass[0] = chain.b0mass[i]
            if pass_jpsimass(label, chain.dimumass[i]):
                s.b0massinjpsi[0] = chain.b0mass[i]
            elif pass_psi2smass(label, chain.dimumass[i]):
                s.b0massinpsi2s[0] = chain.b0mass[i]
            else:
                s.b0massout[0] = chain.b0mass[i]
                
        t.Fill()
        nsel += 1
        if nsel > 100 and test:
            break

    f.Write()
    f.Close()

    if pbar:
        pb.finish()

    return ntot, nsel, selfile 


def output_single_cand_b0_kstmumu(chain, datatype, label, job, test=False,
                                  pbar=False):
    selfile = set_selfile(datatype, label, 'tree', job, test=test)
    f = TFile.Open(selfile, 'RECREATE')
    t = TTree('sel', 'sel')
    s = SelTree(t)

    nsel = 0 
    ntot = 0
    entries = chain.GetEntries()
    if test:
        entries = 300#00

    if pbar:
        pb = get_progressbar(maxval=entries)

    for i in xrange(entries):
        ntot += 1
        if pbar:
            pb.update(i+1)
            
        ientry = chain.LoadTree(i)
        if ientry < 0:
            break

        nb = chain.GetEntry(i)
        if nb <= 0:
            continue

        bMass = chain.bMass

        # will use ConfigParser later... see crb/__init__.py .
        max_mumuVtxCL = 0.05
        
        for i in range(len(bMass)):
            mumuVtxCL = chain.mumuVtxCL[i]
            #print i, mumuVtxCL
            
            if mumuVtxCL >  max_mumuVtxCL:
                s.run[0] = chain.runN

        t.Fill()
        nsel += 1
        if nsel > 100 and test:
            break

    f.Write()
    f.Close()

    if pbar:
        pb.finish()

    return ntot, nsel, selfile 


class SelTree(object):
    '''Class to make Selection Tree  '''
    def __init__(self, t):
        self.run = array('i', [0])
        t.Branch('run', self.run, 'run/I')

        self.event = array('i', [0])        
        t.Branch('event', self.event, 'event/I')

        self.b0mass = array('d', [0.])
        t.Branch('b0mass', self.b0mass, 'b0mass/D')

        self.b0massinjpsi = array('d', [0.])
        t.Branch('b0massinjpsi', self.b0massinjpsi, 'b0massinjpsi/D')

        self.b0massinpsi2s = array('d', [0.])
        t.Branch('b0massinpsi2s', self.b0massinpsi2s, 'b0massinpsi2s/D')

        self.b0massout = array('d', [0.])
        t.Branch('b0massout', self.b0massout, 'b0massout/D')

        self.oniamass = array('d', [0.])
        t.Branch('oniamass', self.oniamass, 'oniamass/D')

        


   
    
