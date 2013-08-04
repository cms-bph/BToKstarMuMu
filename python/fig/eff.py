"""
Module for Efficiency Figures

"""

__author__ = "Xin Shi <Xin.Shi@cern.ch>"
__copyright__ = "Copyright (c) Xin Shi"

import sys 
from tls import *
from array import array
from ROOT import TH1F, TCanvas, TClonesArray, AddressOf, TLorentzVector
from atr.cuts import select_b0s

def main(args, figname):
    if args[0] == 'q2mumu':
        q2mumu(args[1:], figname)
    else:
        raise NameError(args)


def q2mumu(args, figname):
    datatype = args[0]
    label = args[1]
    test = get_options(args, 'test')
    batch = get_options(args, 'batch')

    if batch:
        cmd = create_batch_cmd()
        bashname = '%s.sh' %figname
        bashfile = create_bashfile_cmd(cmd, bashname, label, test=test)
        logfile = set_logfile('fig', datatype, label, figname)
        jobname = 'effq2mm'
        bsub_jobs(logfile, jobname, bashfile, test)
        return


    figfile = set_figfile(figname, label, '.pdf', test=test)
    rootfile = atr.rootfile(datatype, label, test=test)
    obj = atr.root_tree_obj(datatype, label)
    chain = root_chain(rootfile, obj)

    canvas = TCanvas("aCanvas", "Canvas", 600, 600)
    #h_mm_gen = TH1F('mumumass_gen', '#mu^{+} #mu^{-} mass', 100, 0, 25)
    #h_mm_reco = TH1F('mumumass_reco', '#mu^{+} #mu^{-} mass', 100, 0, 25)
    lower = array('f', [0, 2, 4.3, 8.68, 10.09, 12.86, 14.18, 16, 19, 25])

    h_mm_gen = TH1F('mumumass_gen', '#mu^{+} #mu^{-} mass', 9, lower)
    h_mm_reco = TH1F('mumumass_reco', '#mu^{+} #mu^{-} mass', 9, lower)

    if 'B2KstarMuMu/RECO_100M_v1.1' in label:
        Gen_muonPos_P4_ = TClonesArray('TLorentzVector')
        chain.SetBranchAddress('Gen_muonPos_P4', AddressOf(Gen_muonPos_P4_))
        Gen_muonNeg_P4_ = TClonesArray('TLorentzVector')
        chain.SetBranchAddress('Gen_muonNeg_P4', AddressOf(Gen_muonNeg_P4_))

        MuPP4_ = TClonesArray('TLorentzVector')
        chain.SetBranchAddress('MuPP4', AddressOf(MuPP4_))
        MuMP4_ = TClonesArray('TLorentzVector')
        chain.SetBranchAddress('MuMP4', AddressOf(MuMP4_))

        KstarP4_ = TClonesArray('TLorentzVector')
        chain.SetBranchAddress('KstarP4', AddressOf(KstarP4_))

    elif 'B2KstarMuMu/RECO_100M_v1.2' in label or \
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

        
    else:
        raise NameError(label)

    ntot = chain.GetEntries()
    if test:
        ntot = 1000
    
    if  'B2KstarMuMu/RECO_100M_v1.2' in label or \
        'B2KstarMuMu/RECO_100M_v1.4' in label or \
        'B2KstarMuMu/RECO_100M_v1.5' in label: 
        cuts_label = '5ifbv2.6.2'
        cuts = select_b0s(cuts_label)

    sys.stdout.write('Processing %s events ...\n' %ntot)
    sys.stdout.flush()

    nfill_gen = 0
    nfill_reco = 0
    
    for i in xrange(ntot):
        chain.LoadTree(i)
        chain.GetEntry(i)

        if len(chain.Gen_muonPos_P4) > 0:
            mup4_gen = chain.Gen_muonPos_P4[0]
            mum4_gen = chain.Gen_muonNeg_P4[0]
            try: 
                mumu_gen = mup4_gen + mum4_gen
            except TypeError:

                continue 
            h_mm_gen.Fill(mumu_gen.M2())
            nfill_gen += 1 
        
        if '/HLT' in label and not cuts.pass_trigger(chain):
            continue

        if 'OfflineHLT' in label and not chain.offline_hlt_passed:
           continue
        
        if 'MCmatched' in label and not chain.mc_matched:
            continue

        if 'B2KstarMuMu/RECO_100M_v1.1' in label and chain.nXcand > 0:
            if label in ['B2KstarMuMu/RECO_100M_v1.1/Kstar'] and \
                not cuts.pass_kstarmass(chain, 0):
                continue

            if label in ['B2KstarMuMu/RECO_100M_v1.1/lxysig'] and \
                not cuts.pass_lxysig(chain, 0):
                continue

            mup4 = chain.MuPP4[0]
            mum4 = chain.MuMP4[0]
            mumu = mup4 + mum4
            h_mm_reco.Fill(mumu.M2())
            nfill_reco += 1 
            
        if 'B2KstarMuMu/RECO_100M_v1.2' in label:
            mup4 = chain.reco_mup_p4
            mum4 = chain.reco_mum_p4
            mumu = mup4 + mum4
            h_mm_reco.Fill(mumu.M2())
            nfill_reco += 1 

        if 'B2KstarMuMu/RECO_100M_v1.4' in label or \
            'B2KstarMuMu/RECO_100M_v1.5' in label:
            h_mm_reco.Fill(mumu_gen.M2())
            nfill_reco += 1 
        
    sys.stdout.write('Filled events: GEN: %s, RECO: %s. \n' %(nfill_gen, nfill_reco))

    hist = h_mm_reco
    hist.Divide(h_mm_gen)
                     
    hist.SetTitle('RECO Efficiency')
    hist.GetXaxis().SetTitle('q^{2} (GeV^{2}/c^{2})')
    hist.Draw()
    canvas.SaveAs(figfile)
    hist.Delete()
