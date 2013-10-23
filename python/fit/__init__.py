"""
Module for A_FB Fitting 

"""

__author__ = "Xin Shi <Xin.Shi@cern.ch>"
__copyright__ = "Copyright (c) Xin Shi"

import sys 
from tls import * 
import atr

def main(args):

    if args[0] == 'toy':
        import toy
        return toy.main(args[1:])

    if args[0] == 'b0mass':
        import b0mass
        return b0mass.main(args[1:])
    
    if args[0] == 'bmass':
        return bmass(args[1:])
        
    from ROOT import (gROOT, TFile, TTree, TCanvas, gDirectory, RooRealVar,
                      RooCategory, RooArgSet, RooDataSet, RooArgList, RooChebychev,
                      RooGaussian, RooAddPdf, RooFit, kTRUE, kFALSE, 
                      kDashed, kDot, kFullCircle, kFullSquare, RooAbsData)
    
    var = args[0]
    datatype = args[1]
    label = args[2]
    test = get_options(args, 'test')
    batch = get_options(args, 'batch')

    set_root_style(PadLeftMargin=0.13, TitleOffsetY=1.4)
    
    if batch:
        cmd = create_batch_cmd()
        bashname = 'fit_%s.sh' % var 
        bashfile = create_bashfile_cmd(cmd, bashname, label, test=test)
        logfile = set_logfile('fit', datatype, label, var)
        jobname = 'fit'
        bsub_jobs(logfile, jobname, bashfile, test)
        return

    selfile = get_selfile(datatype, label)
    t = root_chain(selfile, 'sel')

    massB0_low = 5.0
    massB0_high = 5.7
    massB0 = RooRealVar(var, "B^{0} mass(GeV/c^{2})", massB0_low, massB0_high)
    ar_massb0 = RooArgSet(massB0)
    dataset = RooDataSet("b0mass","B0 mass", t, ar_massb0)

    c = TCanvas('c', 'canvas', 300, 300)
        
    if var == 'b0massinjpsi':
        title = 'B^{0}#rightarrowJ/#psiK*'
        pars = fit_gau2_che(massB0, dataset, title, test=test)
    elif var == 'b0massout':
        title = 'B^{0}#rightarrowK*#mu^{+}#mu^{-}'
        if label == '5ifbv2.3.1':
            mean = 5.279
            sigma = 0.0586
            sigma1 = 0.03
            sigmaFraction = None
            
        elif label == '5ifbv2.3.2':
            mean = 5.2786
            sigma = 0.05927
            sigma1 = 0.0300
            sigmaFraction = 0.4154

        elif label == '5ifbv2.6.1':
            mean = 5.2788
            sigma = 0.0615
            sigma1 = 0.0300
            sigmaFraction = 0.4022
            
        else:
            raise Exception(label)
        
        pars = fit_gau2_che(massB0, dataset, title, test=test,
                            mean_=mean, sigma_=sigma, sigma1_=sigma1,
            sigmaFraction_ = sigmaFraction)
    else:
        raise Exception(var)

    figname = 'fit_%s_%s' % (var, datatype) 
        
    pdffile = set_figfile(figname, label, '.pdf', test=test)
    txtfile = set_figfile(figname, label, '.txt', test=test)

    c.SaveAs(pdffile)

    save_fit_results(pars, txtfile, verbose=1)


def fit_gau2_che(var, dataset, title='', print_pars=False, test=False,
                 mean_=None, sigma_=None, sigma1_=None, sigmaFraction_=None):
    # define background

    c0 = RooRealVar('c0', 'constant', 0.1, -1, 1)
    c1 = RooRealVar('c1', 'linear', 0.6, -1, 1)
    c2 = RooRealVar('c2', 'quadratic', 0.1, -1, 1)
    c3 = RooRealVar('c3', 'c3', 0.1, -1, 1)

    bkg = RooChebychev('bkg', 'background pdf', var,
                       RooArgList(c0, c1, c2, c3))
    
    # define signal
    val = 5.28
    dmean = 0.05 
    valL = val - dmean
    valR = val + dmean

    if mean_ is None:
        mean = RooRealVar("mean", "mean", val, valL, valR)
    else:
        mean = RooRealVar("mean", "mean", mean_)

    val = 0.05
    dmean = 0.02
    valL = val - dmean
    valR = val + dmean

    if sigma_ is None:
        sigma = RooRealVar('sigma', 'sigma', val, valL, valR)
    else:
        sigma = RooRealVar('sigma', 'sigma', sigma_)

    if sigma1_ is None:
        sigma1 = RooRealVar('sigma1', 'sigma1', val, valL, valR)
    else:
        sigma1 = RooRealVar('sigma1', 'sigma1', sigma1_)

    peakGaus = RooGaussian("peakGaus", "peakGaus", var, mean, sigma)
    peakGaus1 = RooGaussian("peakGaus1", "peakGaus1", var, mean, sigma1)    
    
    if sigmaFraction_ is None:
        sigmaFraction = RooRealVar("sigmaFraction", "Sigma Fraction", 0.5, 0., 1.)
    else:
        sigmaFraction = RooRealVar("sigmaFraction", "Sigma Fraction", sigmaFraction_)

    glist = RooArgList(peakGaus, peakGaus1)
    peakG = RooAddPdf("peakG","peakG", glist, RooArgList(sigmaFraction))
    
    listPeak = RooArgList("listPeak")
    
    listPeak.add(peakG)
    listPeak.add(bkg)
    
    fbkg = 0.45
    nEntries = dataset.numEntries()

    val=(1-fbkg)* nEntries
    listArea = RooArgList("listArea")
    
    areaPeak = RooRealVar("areaPeak", "areaPeak", val, 0.,nEntries)
    listArea.add(areaPeak)

    nBkg = fbkg*nEntries
    areaBkg = RooRealVar("areaBkg","areaBkg", nBkg, 0.,nEntries)
    
    listArea.add(areaBkg)
    model = RooAddPdf("model", "fit model", listPeak, listArea)

    if not test:
        fitres = model.fitTo(dataset, RooFit.Extended(kTRUE),
                             RooFit.Minos(kTRUE),RooFit.Save(kTRUE))

    nbins = 35
    frame = var.frame(nbins)

    frame.GetXaxis().SetTitle("B^{0} mass (GeV/c^{2})")
    frame.GetXaxis().CenterTitle()
    frame.GetYaxis().CenterTitle()
    frame.SetTitle(title)

    mk_size = RooFit.MarkerSize(0.3)
    mk_style = RooFit.MarkerStyle(kFullCircle)
    dataset.plotOn(frame, mk_size, mk_style)

    model.plotOn(frame)

    as_bkg = RooArgSet(bkg)
    cp_bkg = RooFit.Components(as_bkg)
    line_style = RooFit.LineStyle(kDashed)
    model.plotOn(frame, cp_bkg, line_style)

    if print_pars:
        fmt = RooFit.Format('NEU')
        lyt = RooFit.Layout(0.65, 0.95, 0.92)
        param = model.paramOn(frame, fmt, lyt)
        param.getAttText().SetTextSize(0.02)
        param.getAttText().SetTextFont(60)
    
    frame.Draw()

    pars = [areaBkg, areaPeak, c0, c1, c2, c3, mean, sigma, sigma1, sigmaFraction]
    return pars


def test_plot():
    c = TCanvas('c', 'canvas', 300, 300)

    x = RooRealVar("x","x",-10,10)
    mean = RooRealVar("mean","mean of gaussian",1,-10,10) 
    sigma = RooRealVar("sigma","width of gaussian",1,0.1,10)

    gauss = RooGaussian("gauss","gaussian PDF",x,mean,sigma) 

    xframe = x.frame(RooFit.Title("Gaussian p.d.f."))
    #gauss.plotOn(xframe)
    sigma.setVal(3)
    
    as_x = RooArgSet(x)
    data = gauss.generate(as_x,10000)
    data.plotOn(xframe, RooFit.MarkerSize(0.6), RooFit.MarkerStyle(20)) 
    xframe.Draw()

    pdffile = 'test.pdf'
    c.SaveAs(pdffile)
    

def bmass(args):
    figname = 'fit_bmass_%s' % set_figname(args) 
    test=option_exists(args, '-t')
    datatype = args[0] 
    label    = args[1] 
    cut      = args[2] 
    procdir = '../plugins'
    outfile = os.path.join(atr.figpath, figname)
    cmd = './fit bmass %s %s %s %s' %(datatype, label, cut, outfile)
    output = proc_cmd(cmd, procdir=procdir, test=test)
    if test: 
        return 

    print output 
    sys.stdout.write('Figure saved as:\n [[%s.pdf][%s]]\n' %(outfile, figname))

    
