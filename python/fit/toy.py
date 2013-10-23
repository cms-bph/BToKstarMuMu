"""
Module for Toy MC Study 

"""

__author__ = "Xin Shi <Xin.Shi@cern.ch>"
__copyright__ = "Copyright (c) Xin Shi"

import sys 
from tls import *
#from atr.toy import *
import atr 
#from ROOT import RooRealVar, RooGaussian

#(gROOT, TFile, TTree, TCanvas, gDirectory, RooRealVar,
#    RooCategory, RooArgSet, RooDataSet, RooArgList, RooChebychev,
#    RooGaussian, RooAddPdf, RooFit, kTRUE, kFALSE, kDashed, kDot, kFullCircle,
#    kFullSquare, RooAbsData)

def main(args):
    binnums = get_range_from_str(args[4])
    
    for binnum in binnums:
        proc_bin(binnum, args)


def proc_bin(binnum, args):
    datatype = args[0]
    label = args[1]
    fittype = args[2]
    effcorr = args[3]

    test = option_exists(args, '-t')
    pbar = option_exists(args, '-pbar')
    batch = option_exists(args, '-b')
    #queue = get_option(args, 'queue')
    jobname = '' #get_option(args, 'jobname')

    parafile  = '../python/ParameterFile.txt'

    if not jobname:
        jobname = '%s%s' %(fittype, binnum)

    figname = 'ToyMC_%s_%s' % (fittype, binnum)
    jobindex = os.getenv('LSB_JOBINDEX')
    if jobindex:
        figname += '_%s' %jobindex

    plotfile = set_file(atr.figpath, label, figname, '.root', test=test)

    if label in ['ntoy4000v1', 'ntoy1000v1', 'ntoy1000v3']:
        figname = 'MassAngleToyMC_%s_%s' % (fittype, binnum)
        plotfile = set_file(atr.figpath, label, figname, '.root', test=test)


    if label in ['ntoy4000v1', 'ntoy1000v1', 'ntoy1000v3', 'ntoy1000v4',
                 'ntoy1000v5', 'ntoy1000_mult100v1']:
        procdir = '/afs/cern.ch/user/x/xshi/work/cms/afb/rel/\
CMSSW_4_2_8_patch7/src/BphAna/B0KstMuMu/plugins' 
        maincode = 'ExtractYield'
        if label == 'ntoy1000_mult100v1':
            maincode = 'ExtractYield_Mult100'
        if label == 'ntoy4000v1':
            ntoy = 4000
        if label in ['ntoy1000v1', 'ntoy1000v3', 'ntoy1000v4', 'ntoy1000v5',
                     'ntoy1000_mult100v1']:
            ntoy = 1000

    elif label in ['ntoy1000v2']:
        run_toy_mc(fittype, effcorr, binnum)
        return
    
    elif label in ['ntoy1000v6', 'ntoy1000_mult100v2']:
        procdir = '/afs/cern.ch/user/x/xshi/work/cms/afb/rel/\
CMSSW_4_2_8_patch7/src/BphAna/B0KstMuMu_V00_04_05/plugins' 
        ntoy = 1000
        maincode = 'ExtractYield'

    elif label in ['ntoy100_mult100v1']:
        procdir = '/afs/cern.ch/user/x/xshi/work/cms/afb/rel/\
CMSSW_4_2_8_patch7/src/BphAna/B0KstMuMu_V00_04_05/plugins' 
        ntoy = 100
        #ntoy = 1
        #sys.stdout.write('Attention! using ntoy = 1 !!! \n')
        maincode = 'ExtractYield'

    elif label in ['ntoy100_mult100v2']:
        procdir = '/afs/cern.ch/user/x/xshi/work/cms/afb/rel/\
CMSSW_4_2_8_patch7/src/BphAna/B0KstMuMu_V00_04_09/plugins' 
        ntoy = 100
        maincode = 'ExtractYield'

    elif label in ['ntoy50_mult100v2']:
        procdir = '/afs/cern.ch/user/x/xshi/work/cms/afb/rel/\
CMSSW_4_2_8_patch7/src/BphAna/B0KstMuMu_V00_04_09/plugins' 
        ntoy = 50
        maincode = 'ExtractYield'

    elif label in ['ntoy1000v8']:
        procdir = '/afs/cern.ch/user/x/xshi/work/cms/afb/rel/\
CMSSW_4_2_8_patch7/src/BphAna/B0KstMuMu_V00_06_03/plugins' 
        ntoy = 1000
        maincode = 'ExtractYield'
        parafile  = '../results/ParameterFile_Data.txt'

    elif label in ['ntoy1000v9']:
        procdir = '/afs/cern.ch/user/x/xshi/work/cms/afb/rel/\
CMSSW_4_2_8_patch7/src/BphAna/B0KstMuMu_V00_06_05/plugins' 
        ntoy = 1000
        maincode = 'ExtractYield'
        parafile  = '../results/ParameterFile_Data.txt'

    elif label in ['ntoy1000v10']:
        procdir = '/afs/cern.ch/user/x/xshi/work/cms/afb/rel/\
CMSSW_4_2_8_patch7/src/BphAna/B0KstMuMu_V00_06_06/plugins' 
        ntoy = 1000
        maincode = 'ExtractYield'
        parafile  = '../results/ParameterFile_Data.txt'

    elif label in ['ntoy1000v11']:
        procdir = '/afs/cern.ch/user/x/xshi/work/cms/afb/rel/\
CMSSW_4_2_8_patch7/src/BphAna/B0KstMuMu_V00_09_01/plugins' 
        ntoy = 1000
        maincode = 'ExtractYield'
        parafile  = '../results/ParameterFile_Data.txt'

    elif label in ['ntoy1000v12']:
        procdir = '/afs/cern.ch/user/x/xshi/work/cms/afb/rel/\
CMSSW_4_2_8_patch7/src/BphAna/B0KstMuMu_V00_09_02/plugins' 
        ntoy = 1000
        maincode = 'ExtractYield'
        parafile  = '../python/ParameterFile.txt'
	
    else:
        raise NameError(label)

    
    if fittype == 'BF':
        fitcode = 21
        if label == 'ntoy1000v1':
            parafile  = '../python/ParameterFile_Data'
        elif label == 'ntoy1000v3':
            parafile  = '../results/ParameterFile_Data'

    elif fittype == 'FL':
        fitcode = 23
        if label == 'ntoy1000v1':
            parafile  = '../python/ParameterFile_Data_FL'
        if label == 'ntoy1000v3':
            parafile  = '../results/ParameterFile_Data_FL'
        if label == 'ntoy1000v6':
            maincode = 'ExtractYield_FL'
        if label == 'ntoy1000_mult100v2':
            maincode = 'ExtractYield_FL_Mult100'
        
    elif fittype == 'AFB':
        fitcode = 24
        if label == 'ntoy1000v5':
            maincode = 'ExtractYield_AFB'
        if label == 'ntoy1000v1':
            parafile  = '../python/ParameterFile_Data_AFB'
        if label == 'ntoy1000v3':
            parafile = '../results/ParameterFile_Data_AFB'
        if label in ['ntoy1000_mult100v1', 'ntoy1000_mult100v2',
                     'ntoy100_mult100v1', 'ntoy100_mult100v2', 'ntoy50_mult100v2']:
            maincode = 'ExtractYield_AFB_Mult100'

    elif fittype == 'FL-AFB':
        fitcode = 26

    else:
        raise NameError(fittype)

    if batch:
        cmd = create_batch_cmd()

        if label in ['ntoy4000v1', 'ntoy1000v1', 'ntoy1000v3', 'ntoy1000v4',
                     'ntoy1000v5', 'ntoy1000_mult100v1', 'ntoy1000v6',
                     'ntoy1000_mult100v2', 'ntoy100_mult100v1', 'ntoy100_mult100v2',
                     'ntoy50_mult100v2', 'ntoy1000v8', 'ntoy1000v9', 'ntoy1000v10',
		     'ntoy1000v11', 'ntoy1000v12']:
            cmd = 'setroot\n' + cmd
        else:
            raise NameError(label)

        if label in ['ntoy1000_mult100v1', 'ntoy1000_mult100v2']:
            queue = '1nw'
            
        bashfile = set_file(atr.bashpath, label, figname, '.sh', test=test)
	update_bashfile_cmd(bashfile, cmd, test=test)
        logfile = set_file(atr.logpath, label, figname, '.log', test=test)
        if jobname and '[' in jobname and ']' in jobname:
            logfile += '.%I'

        bsub_jobs(logfile, jobname, bashfile, test=test, queue=queue)
        return


    cmd = './%s %s %s %s %s %s %s' %(
        maincode, fitcode, plotfile, effcorr, binnum, ntoy, parafile)

    if label in ['ntoy4000v1', 'ntoy1000v1', 'ntoy1000v3']:
        cmd = './%s %s %s %s --ntoy %s --plotname %s --parafile %s' %(
            maincode, fitcode, effcorr, binnum, ntoy, plotfile, parafile)

    proc_cmd(cmd, test, procdir=procdir)

def run_toy_mc(fittype, effcorr, binnum):
    # Now make TOY-MC for fit to B0 total invariant mass and cos(theta_L)
    #from ROOT import RooRealVar, RooGaussian
    B0MassArb = RooRealVar("B0MassArb","B0 mass arbitrated",
                           B0Mass - B0MassIntervalLeft,
                           B0Mass + B0MassIntervalRight, "GeV/c^{2}")
     
    instantiate_mass_angle_fit(B0MassArb)
    sys.exit()    
    #  (RooAbsPdf** TotalPDF,
    #  bool correctByPDF,
    #  RooRealVar* x,
    #  RooRealVar* y,
    #  string fitName,
    #  unsigned int FitType,
    #  unsigned int PsiIntrusion,
    #  bool use2GaussS,
    #  bool use2BExp,
    #  vector<double>* cosThetaKBins, vector<double>* cosThetaLBins,
    #  TF2* effFunc,
    #  bool isFirstDataset)

    # instantiate_mass_angle_fit(TotalPDFRejectPsi,
    #                            useEffPDF,
    #                            B0MassArb,
    #                            CosThetaMuArb,
    #                            "TotalPDFRejectPsi",
    #                            FitType,
    #                            configParam[0]->operator[](specBin),
    #                            configParam[1]->operator[](specBin),
    #                            configParam[2]->operator[](specBin),
    #                            cosThetaKBins,&cosThetaLBins, effFuncs1[specBin],true) 

    
def instantiate_mass_angle_fit(x):
    '''
    2D Model: 
    x: mass                              
    y: angle cos(theta_L) OR cos(theta_K)
    '''


    # Define mass fit variables and pdf for signal
    meanS = RooRealVar("meanS","Signal mean", B0Mass,"GeV/c^{2}")
    sigmaS1 = RooRealVar("sigmaS1","Signal sigma-1",SIGMAS1,"GeV/c^{2}")
    MassS1 = RooGaussian("MassS1","Signal Gaussian-1", x, meanS, sigmaS1)
    sigmaS2 = RooRealVar("sigmaS2","Signal sigma-2",SIGMAS2,"GeV/c^{2}")

    MassS2 = RooGaussian("MassS2","Signal Gaussian-2", x, meanS, sigmaS2)
    meanS.setConstant(True)
    sigmaS1.setConstant(True)
    sigmaS2.setConstant(True)

    fracMassS = RooRealVar("fracMassS","Fraction of signal Gaussian",
                           FRACMASSS,0.0,1.0)
    fracMassS.setConstant(True)

    print 'here' 
    sys.exit()

    
    
    
