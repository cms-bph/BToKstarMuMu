"""
Module for A_FB EVT Skim 

"""

__author__ = "Xin Shi <Xin.Shi@cern.ch>"
__copyright__ = "Copyright (c) Xin Shi"


import sys 
from tls import *
from tls.cuts import *
from ROOT import (gROOT, TObject, TClonesArray, AddressOf, TLorentzVector,
                  TFile, TTree, kTRUE)
from time import time
from array import array


def main(args):
    datatype = args[0]
    label = args[1]
    job_opt = get_options(args, 'job')
    batch = get_options(args, 'batch')
    test = get_options(args, 'test')
    pbar = get_options(args, 'pbar')
    stage = 'ntp'

    if job_opt == 'create':
        dbfile = create_jobs_db(stage, datatype, label)
        return

    if job_opt == 'update':
        dbfile = update_jobs_db(stage, datatype, label)
        return
    
    elif job_opt == 'stager_get' or job_opt == 'stager_qry':
        cmdname = job_opt
        jobs = get_jobs(stage, datatype, label)
        range_ = range(1, len(jobs)+1)
        for job in range_:
            j = jobs[str(job)]
            status = j['status']
            rootfile = j['rootfile']
            if label in ['5ifbv2.2']:
                rootfiles = [r.encode('ascii') for r in rootfile]
            else:
                rootfiles = [r[0].encode('ascii') for r in rootfile]

            for rootfile in rootfiles:
                stager_cmd(cmdname, rootfile, test)


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

        sys.stdout.write('| Job | Processed | Skimmed | Duration |\n')
        sys.stdout.write('|-----+-----------+---------+----------|\n')

        total_processed = 0
        total_skimmed = 0 
        for job in range_:
          job, processed, skimmed, duration = check_job_log(
            'skm', datatype, label, job)
          total_processed += int(processed)
          total_skimmed += int(skimmed)
          
        sys.stdout.write('|-----+-----------+---------+----------|\n')
        sys.stdout.write('Total: processed %s , skimmed %s, (%s) \n' %(
        total_processed, total_skimmed, total_skimmed/float(total_processed)))
        
    else:
        raise NameError(job_opt)
        

def proc_one_job(rootfile, datatype, label, job=None, test=False, batch=False,
                 pbar=False):
    if batch:
        cmd = create_batch_cmd(job)
        bashname = 'skm_evt_%s.sh' %job 
        bashfile = create_bashfile_cmd(cmd, bashname, test=test)
        logfile = set_logfile('skm', datatype, label, 'evt', job)
        jobname = 'skm%s' % job
        bsub_jobs(logfile, jobname, bashfile, test)
        return

    chain = get_chain(rootfile, 'data')
    time_start = time()
    ntot, nskm, skmfile = output_evt_mass_b0(chain, datatype, label, job, test=test,
                                             pbar=pbar)
    dur = duration(time() - time_start)
    
    sys.stdout.write(' processed %s \n skimmed %s \n duration %s \n saved as %s\n' % (
        ntot, nskm, dur, skmfile))
    sys.stdout.flush()


def output_evt_mass_b0(chain, datatype, label, job, test=False, pbar=False):
    skmfile = set_skmfile(datatype, label, 'skim', job)
    if test:
        skmfile = skmfile + '.test'

    f = TFile.Open(skmfile, 'RECREATE')
    t = TTree('skm', 'skm')

    s = SkimTree(t)

    nskm = 0 
    ntot = 0

    xp4_ = TClonesArray('TLorentzVector')
    chain.SetBranchAddress('xP4', AddressOf(xp4_))

    kstarp4_ = TClonesArray('TLorentzVector')
    chain.SetBranchAddress('KstarP4', AddressOf(kstarp4_))

    entries = chain.GetEntries()

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

        if not pass_trigger(chain, label):
            continue

        if chain.nXcand < 1:
            continue

        for nx in range(chain.nXcand):
            b0p4 = chain.xP4[nx]
            b0mass = b0p4.M()
            lxysig = chain.xlxySig[nx]
            input_var =  (b0mass, lxysig)
            if not pass_b0_skim(input_var, label):
                continue

            fill_skim_tree(s, chain, nx)

        t.Fill()
        nskm += 1
        if nskm > 100 and test:
            break

    f.Write()
    f.Close()

    if pbar:
        pb.finish()

    return ntot, nskm, skmfile 


class SkimTree(object):
    '''Class to make Skim Tree  '''
    def __init__(self, t):
        max_nxcand = 10000

        self.run = array('i', [0])
        t.Branch('run', self.run, 'run/I')

        self.event = array('i', [0])        
        t.Branch('event', self.event, 'event/I')

        self.nxcand = array('i', [0])
        t.Branch('nxcand', self.nxcand, 'nxcand/I')
        
        self.b0mass = array('d', max_nxcand*[0.])
        t.Branch('b0mass', self.b0mass, 'b0mass[nxcand]/D')

        self.b0pt = array('d', max_nxcand*[0.])
        t.Branch('b0pt', self.b0pt, 'b0pt[nxcand]/D')
        
        self.kstarmass = array('d', max_nxcand*[0.])
        t.Branch('kstarmass', self.kstarmass, 'kstarmass[nxcand]/D')

        self.kstarpt = array('d', max_nxcand*[0.])
        t.Branch('kstarpt', self.kstarpt, 'kstarpt[nxcand]/D')
        
        self.dimumass = array('d', max_nxcand*[0.])
        t.Branch('dimumass', self.dimumass, 'dimumass[nxcand]/D')

        self.dimupt = array('d', max_nxcand*[0.])
        t.Branch('dimupt', self.dimupt, 'dimupt[nxcand]/D')

        self.vtxcl = array('d', max_nxcand*[0.])
        t.Branch('vtxcl', self.vtxcl, 'vtxcl[nxcand]/D')        

        self.cosalpha = array('d',  max_nxcand*[0.])
        t.Branch('cosalpha', self.cosalpha, 'cosalpha[nxcand]/D')        
        
        self.lxysig = array('d', max_nxcand*[0.])
        t.Branch('lxysig', self.lxysig, 'lxysig[nxcand]/D')        

        self.ctaupv = array('d', max_nxcand*[0.])
        t.Branch('ctaupv', self.ctaupv, 'ctaupv[nxcand]/D')        

        self.oniadca = array('d', max_nxcand*[0.])
        t.Branch('oniadca', self.oniadca, 'oniadca[nxcand]/D')        

        self.kstardca = array('d', max_nxcand*[0.])
        t.Branch('kstardca', self.kstardca, 'kstardca[nxcand]/D')        



def fill_skim_tree(s, ch, nx):
    b0p4 = ch.xP4[nx]
    kstarp4 = ch.KstarP4[nx]
    oniap4 = b0p4 - kstarp4

    s.run[0] = ch.runNb
    s.event[0] = ch.eventNb        
    s.nxcand[0] = ch.nXcand

    s.b0mass[nx] = b0p4.M()
    s.b0pt[nx] = b0p4.Pt()
    
    s.kstarmass[nx] = kstarp4.M()
    s.kstarpt[nx] = kstarp4.Pt()    
    
    s.dimumass[nx] = oniap4.M()
    s.dimupt[nx] = oniap4.Pt()

    s.vtxcl[nx] = ch.xVtxCL[nx]
    s.cosalpha[nx] = ch.xcosAlpha[nx]
    s.lxysig[nx] = ch.xlxySig[nx]
    s.ctaupv[nx] = ch.xctauPV[nx] 
    s.oniadca[nx] = ch.OniaDCA[nx]
    s.kstardca[nx] = ch.kstarDCA[nx]   

