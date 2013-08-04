"""
Module for pT distritubions Figures

"""

__author__ = "Xin Shi <Xin.Shi@cern.ch>"
__copyright__ = "Copyright (c) Xin Shi"

import sys 
from tls import *
from array import array
from ROOT import TH1F, TCanvas, TClonesArray, AddressOf, TLorentzVector, TVector3
from atr.cuts import select_b0s


def main(args, figname):
    if args[0] == 'mup':
        mup(args[1:], figname)
    elif args[0] == 'kaon':
        kaon(args[1:], figname)
    elif args[0] == 'mue':
        mue(args[1:], figname)
    elif args[0] == 'pion':
        pion(args[1:], figname)
    else:
        raise NameError(args)

def mup(args, figname):
    datatype = args[0]
    label = args[1]
    test = get_options(args, 'test')
    batch = get_options(args, 'batch')

    if batch:
        cmd = create_batch_cmd()
        bashname = '%s.sh' %figname
        bashfile = create_bashfile_cmd(cmd, bashname, label, test=test)
        logfile = set_logfile('fig', datatype, label, figname)
        jobname = 'figptmu'
        bsub_jobs(logfile, jobname, bashfile, test)
        return

    figfile = set_figfile(figname, label, '.pdf', test=test)
    rootfile = atr.rootfile(datatype, label, test=test)

    obj = atr.root_tree_obj(datatype, label)
    chain = root_chain(rootfile, obj)

    canvas = TCanvas("aCanvas", "Canvas", 600, 600)
    hist = TH1F('mupt', '#mu^{+} p_{T}', 100, 0, 20)

    if  'B2KstarMuMu/RECO_100M_v1.1' in label or \
        'B2KstarMuMu/RECO_100M_v1.2' in label:
        Gen_muonPos_P4_ = TClonesArray('TLorentzVector')
        chain.SetBranchAddress('Gen_muonPos_P4', AddressOf(Gen_muonPos_P4_))
        
        MuPP4_ = TClonesArray('TLorentzVector')
        chain.SetBranchAddress('MuPP4', AddressOf(MuPP4_))

    if 'B2KstarMuMu/RECO_100M_v1.3' in label or \
        'B2KstarMuMu/RECO_100M_v1.4' in label or \
        'B2KstarMuMu/RECO_100M_v1.5' in label:
        Gen_muonPos_P4_ = TClonesArray('TLorentzVector')
        chain.SetBranchAddress('Gen_muonPos_P4', AddressOf(Gen_muonPos_P4_))
        Gen_muonNeg_P4_ = TClonesArray('TLorentzVector')
        chain.SetBranchAddress('Gen_muonNeg_P4', AddressOf(Gen_muonNeg_P4_))

        reco_mup_p4_ = TLorentzVector()
        chain.SetBranchAddress('reco_mup_p4', AddressOf(reco_mup_p4_))

        reco_mum_p4_ = TLorentzVector()
        chain.SetBranchAddress('reco_mum_p4', AddressOf(reco_mum_p4_))
    
    ntot = chain.GetEntries()
    if test:
        ntot = 1000

    if label in ['B2KstarMuMu/RECO_100M_v1.1/HLT',
                 'B2KstarMuMu/RECO_100M_v1.2/HLT',
                 'B2KstarMuMu/RECO_100M_v1.3/HLT',
                 'B2KstarMuMu/RECO_100M_v1.3/MCmatched/HLT',
                 'B2KstarMuMu/RECO_100M_v1.4/MCmatched/HLT',
                 'B2KstarMuMu/RECO_100M_v1.5/MCmatched/OfflineHLT',
                 ]:
        cuts_label = '5ifbv2.6.2'
        cuts = select_b0s(cuts_label)
        
    sys.stdout.write('Processing %s events ...\n' %ntot)
    sys.stdout.flush()
    nfill = 0 
    for i in xrange(ntot):
        chain.LoadTree(i)
        chain.GetEntry(i)

        if label in ['B2KstarMuMu/RECO_100M_v1.1/HLT',
                     'B2KstarMuMu/RECO_100M_v1.2/HLT',
                     'B2KstarMuMu/RECO_100M_v1.3/HLT',
                     'B2KstarMuMu/RECO_100M_v1.3/MCmatched/HLT', 
                     'B2KstarMuMu/RECO_100M_v1.4/MCmatched/HLT',
                     ] and \
                     not cuts.pass_trigger(chain):
           continue

        if label in ['B2KstarMuMu/RECO_100M_v1.5/MCmatched/OfflineHLT',
                     ] and \
                     not chain.offline_hlt_passed:
           continue

        if 'GEN' in label:
            if 'B2KstarMuMu/RECO_100M_v1.5' in label:
                # need others if necessary for backward compatibility.
                mup4 = chain.Gen_muonPos_P4[0]
            elif 'B2KstarMuMu/RECO_100M_v1.6' in label:
                mup4 = TVector3(chain.genMupPx, chain.genMupPy, chain.genMupPz)

            else:
                raise NameError(label)
                
            
        if 'MCmatched' in label and not chain.mc_matched:
            continue

        if 'B2KstarMuMu/RECO_100M_v1.1' in label:
            if chain.nXcand <= 0:
                continue
            mup4 = chain.MuPP4[0]
            
        if 'B2KstarMuMu/RECO_100M_v1.2' in label or \
            'B2KstarMuMu/RECO_100M_v1.3' in label or \
            'B2KstarMuMu/RECO_100M_v1.4' in label or \
            'B2KstarMuMu/RECO_100M_v1.5' in label:
            mup4 = chain.reco_mup_p4

        hist.Fill(mup4.Pt())
        nfill += 1
        
    sys.stdout.write('Filled %s events. \n' %nfill)
    hist.GetXaxis().SetTitle('p_{T} (GeV/c)')
    hist.Draw()
    canvas.SaveAs(figfile)
    hist.Delete()
    
        
def mue(args, figname):
    datatype = args[0]
    label = args[1]
    test = get_options(args, 'test')
    batch = get_options(args, 'batch')

    if batch:
        cmd = create_batch_cmd()
        bashname = '%s.sh' %figname
        bashfile = create_bashfile_cmd(cmd, bashname, label, test=test)
        logfile = set_logfile('fig', datatype, label, figname)
        jobname = 'figptmu'
        bsub_jobs(logfile, jobname, bashfile, test)
        return

    figfile = set_figfile(figname, label, '.pdf', test=test)
    rootfile = atr.rootfile(datatype, label, test=test)
    obj = atr.root_tree_obj(datatype, label)
    chain = root_chain(rootfile, obj)
    canvas = TCanvas("aCanvas", "Canvas", 600, 600)
    hist = TH1F('mupt', '#mu^{-} p_{T}', 100, 0, 20)
 
    Gen_muonNeg_P4_ = TClonesArray('TLorentzVector')
    chain.SetBranchAddress('Gen_muonNeg_P4', AddressOf(Gen_muonNeg_P4_))

    MuMP4_ = TClonesArray('TLorentzVector')
    chain.SetBranchAddress('MuMP4', AddressOf(MuMP4_))

    ntot = chain.GetEntries()
    if test:
        ntot = 1000

    sys.stdout.write('Processing %s events ...\n' %ntot)
    sys.stdout.flush()
    nfill = 0 
    for i in xrange(ntot):
        chain.LoadTree(i)
        chain.GetEntry(i)
        if 'GEN' in label:
            mup4 = chain.Gen_muonNeg_P4[0]
        else:
            if chain.nXcand <= 0:
                continue
            mup4 = chain.MuMP4[0]

        hist.Fill(mup4.Pt())
        nfill += 1
        
    sys.stdout.write('Filled %s events. \n' %nfill)
    hist.GetXaxis().SetTitle('p_{T} (GeV/c)')
    hist.Draw()
    canvas.SaveAs(figfile)
    hist.Delete()


def kaon(args, figname):
    datatype = args[0]
    label = args[1]
    test = get_options(args, 'test')
    batch = get_options(args, 'batch')

    if batch:
        cmd = create_batch_cmd()
        bashname = '%s.sh' %figname
        bashfile = create_bashfile_cmd(cmd, bashname, label, test=test)
        logfile = set_logfile('fig', datatype, label, figname)
        jobname = 'figptkaon'
        bsub_jobs(logfile, jobname, bashfile, test)
        return

    figfile = set_figfile(figname, label, '.pdf', test=test)
    rootfile = atr.rootfile(datatype, label, test=test)
    obj = atr.root_tree_obj(datatype, label)
    chain = root_chain(rootfile, obj)

    canvas = TCanvas("aCanvas", "Canvas", 600, 600)
    hist = TH1F('kaonpt', 'K p_{T}', 100, 0, 20)
 
    Gen_Kaon_P4_ = TClonesArray('TLorentzVector')
    chain.SetBranchAddress('Gen_Kaon_P4', AddressOf(Gen_Kaon_P4_))

    KaonPP4_ = TClonesArray('TLorentzVector')
    chain.SetBranchAddress('KaonPP4', AddressOf(KaonPP4_))

    ntot = chain.GetEntries()
    if test:
        ntot = 1000

    sys.stdout.write('Processing %s events ...\n' %ntot)
    sys.stdout.flush()
    nfill = 0 
    for i in xrange(ntot):
        chain.LoadTree(i)
        chain.GetEntry(i)

        if 'GEN' in label:
            kaonp4 = chain.Gen_Kaon_P4[0]
        else:
            if chain.nXcand <= 0:
                continue
            kaonp4 = kaonp4 = chain.KaonPP4[0]
        hist.Fill(kaonp4.Pt())
        nfill += 1
        
    sys.stdout.write('Filled %s events. \n' %nfill)
    hist.GetXaxis().SetTitle('p_{T} (GeV/c)')
    hist.Draw()
    canvas.SaveAs(figfile)
    hist.Delete()
        
def pion(args, figname):
    datatype = args[0]
    label = args[1]
    test = get_options(args, 'test')
    batch = get_options(args, 'batch')

    if batch:
        cmd = create_batch_cmd()
        bashname = '%s.sh' %figname
        bashfile = create_bashfile_cmd(cmd, bashname, label, test=test)
        logfile = set_logfile('fig', datatype, label, figname)
        jobname = 'figptpion'
        bsub_jobs(logfile, jobname, bashfile, test)
        return

    figfile = set_figfile(figname, label, '.pdf', test=test)
    rootfile = atr.rootfile(datatype, label, test=test)
    obj = atr.root_tree_obj(datatype, label)
    chain = root_chain(rootfile, obj)

    canvas = TCanvas("aCanvas", "Canvas", 600, 600)
    hist = TH1F('pionpt', '#pi p_{T}', 100, 0, 20)
 
    Gen_Pion_P4_ = TClonesArray('TLorentzVector')
    chain.SetBranchAddress('Gen_Pion_P4', AddressOf(Gen_Pion_P4_))

    PionPP4_ = TClonesArray('TLorentzVector')
    chain.SetBranchAddress('PionPP4', AddressOf(PionPP4_))

    ntot = chain.GetEntries()
    if test:
        ntot = 1000

    sys.stdout.write('Processing %s events ...\n' %ntot)
    sys.stdout.flush()
    nfill = 0 
    for i in xrange(ntot):
        chain.LoadTree(i)
        chain.GetEntry(i)

        if 'GEN' in label:
            pionp4 = chain.Gen_Pion_P4[0]
        else:
            if chain.nXcand <= 0:
                continue
            pionp4 = pionp4 = chain.PionPP4[0]
        hist.Fill(pionp4.Pt())
        nfill += 1
        
    sys.stdout.write('Filled %s events. \n' %nfill)
    hist.GetXaxis().SetTitle('p_{T} (GeV/c)')
    hist.Draw()
    canvas.SaveAs(figfile)
    hist.Delete()
        
