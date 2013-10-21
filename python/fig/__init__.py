"""
Module for A_FB Figures 

"""

__author__ = "Xin Shi <Xin.Shi@cern.ch>"
__copyright__ = "Copyright (c) Xin Shi"

import sys 
from tls import * 
import atr
#from atr.cuts import select_b0s
from ROOT import TH1F, TCanvas, TClonesArray, AddressOf, TLorentzVector

def main(args):
    set_root_style()
    figname = set_figname(args)

    if args[0] == 'q2mumu':
        q2mumu(args[1:], figname)
    elif args[0] == 'eff':
        import eff
        eff.main(args[1:], figname)
    elif args[0] == 'pt':
        import pt
        pt.main(args[1:], figname)
    elif args[0] == 'SingleBuToKstarMuMu':
        import SingleBuToKstarMuMu
        SingleBuToKstarMuMu.main(args[1:], figname)
    elif args[0] == 'bmass':
        bmass(args[1:], figname)
    elif args[0] == 'bmass_snb':
        bmass_snb(args[1:], figname)
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
        jobname = 'figq2mm'
        bsub_jobs(logfile, jobname, bashfile, test)
        return

 
    figfile = set_figfile(figname, label, '.pdf', test=test)

    #rootfile = get_rootfile(datatype, label)
    rootfile = atr.rootfile(datatype, label)
    
    obj = atr.root_tree_obj(datatype, label)
    chain = root_chain(rootfile, obj)

    canvas = TCanvas("aCanvas", "Canvas", 600, 600)
    hist = TH1F('mumumass', '#mu^{+} #mu^{-} mass', 100, 0, 25)

    if label in ['B2KstarMuMu/RECO_1M_v2.2',
                 'B2KstarMuMu/RECO_1M_v2.3',
                 'B2KstarMuMu/RECO_100M_v1.1']:
        MuPP4_ = TClonesArray('TLorentzVector')
        chain.SetBranchAddress('MuPP4', AddressOf(MuPP4_))
        MuMP4_ = TClonesArray('TLorentzVector')
        chain.SetBranchAddress('MuMP4', AddressOf(MuMP4_))
        
        Gen_muonPos_P4_ = TClonesArray('TLorentzVector')
        chain.SetBranchAddress('Gen_muonPos_P4', AddressOf(Gen_muonPos_P4_))
        Gen_muonNeg_P4_ = TClonesArray('TLorentzVector')
        chain.SetBranchAddress('Gen_muonNeg_P4', AddressOf(Gen_muonNeg_P4_))

    elif 'B2KstarMuMu/RECO_100M_v1.2' in label or 'B2KstarMuMu/RECO_100M_v1.4' in label: 
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

    if  'B2KstarMuMu/RECO_100M_v1.1' in label or \
        'B2KstarMuMu/RECO_100M_v1.4' in label: 
        cuts_label = '5ifbv2.6.2'
        cuts = select_b0s(cuts_label)

    sys.stdout.write('Processing %s events ...\n' %ntot)
    sys.stdout.flush()

    nfill = 0 
    for i in xrange(ntot):
        chain.LoadTree(i)
        chain.GetEntry(i)

        if label == 'B2KstarMuMu/GEN_1M_v1.1': 
            mup4 = fourvecs_xyzm(chain.genMupPx, chain.genMupPy,
                                 chain.genMupPz, atr.muon_mass)
            mum4 = fourvecs_xyzm(chain.genMumPx, chain.genMumPy,
                                 chain.genMumPz, atr.muon_mass)

            mumu = mup4 + mum4
            hist.Fill(mumu.M2())
            nfill += 1 
            
        elif label in ['B2KstarMuMu/RECO_1M_v2.2',
                       'B2KstarMuMu/RECO_1M_v2.3', 
                       'B2KstarMuMu/RECO_100M_v1.1']:
            if chain.nXcand > 0:
                mup4 = chain.MuPP4[0]
                mum4 = chain.MuMP4[0]
                mumu = mup4 + mum4
                hist.Fill(mumu.M2())
                nfill += 1 

        elif label in ['B2KstarMuMu/RECO_1M_v2.2/GEN',
                       'B2KstarMuMu/RECO_1M_v2.3/GEN',
                       'B2KstarMuMu/RECO_100M_v1.1/GEN',
                       'B2KstarMuMu/RECO_100M_v1.1/GEN_HLT']:

            if label in ['B2KstarMuMu/RECO_100M_v1.1/GEN_HLT']:
                if not cuts.pass_trigger(chain):
                    continue

            if len(chain.Gen_muonPos_P4) > 0:
                mup4 = chain.Gen_muonPos_P4[0]
                mum4 = chain.Gen_muonNeg_P4[0]

                try: 
                    mumu = mup4 + mum4
                except TypeError:
                    continue 
                    # mup4 = chain.Gen_muonPos_P4[1]
                    # mum4 = chain.Gen_muonNeg_P4[1]
                    # mumu = mup4 + mum4
                hist.Fill(mumu.M2())
                nfill += 1
                
        elif 'B2KstarMuMu/RECO_100M_v1.2' in label or \
            'B2KstarMuMu/RECO_100M_v1.4' in label:

            if 'MCmatched' in label and not chain.mc_matched:
                continue

            mup4 = chain.reco_mup_p4
            mum4 = chain.reco_mum_p4
            
            if '/GEN' in label and len(chain.Gen_muonPos_P4) > 0:
                mup4 = chain.Gen_muonPos_P4[0]
                mum4 = chain.Gen_muonNeg_P4[0]

            mumu = mup4 + mum4
            hist.Fill(mumu.M2())
            nfill += 1 

        else:
            raise NameError(label)
                
    sys.stdout.write('Filled %s events. \n' %nfill)
    hist.GetXaxis().SetTitle('q^{2} (GeV^{2}/c^{2})')
    hist.Draw()
    canvas.SaveAs(figfile)
    hist.Delete()


def bmass(args, figname):
    test=option_exists(args, '-t')
    datatype = args[0] 
    label    = args[1] 
    cut      = args[2] 
    procdir = '../plugins'
    outfile = os.path.join(atr.figpath, figname+'.pdf')
    cmd = './fig bmass %s %s %s %s' %(datatype, label, cut, outfile)
    output = proc_cmd(cmd, procdir=procdir, test=test)
    if test: 
        return 

    print output 
    sys.stdout.write('Figure saved as:\n [[%s][%s]]\n' %(outfile, figname))

    
def bmass_snb(args, figname):
    test=option_exists(args, '-t')
    datatype = args[0] 
    label    = args[1] 
    cut      = args[2] 
    procdir = '../plugins'
    outfile = os.path.join(atr.figpath, figname+'.pdf')
    cmd = './fig bmass_snb %s %s %s %s' %(datatype, label, cut, outfile)
    output = proc_cmd(cmd, procdir=procdir, test=test)
    if test: 
        return 

    print output 
    sys.stdout.write('Figure saved as:\n [[%s][%s]]\n' %(outfile, figname))

    
