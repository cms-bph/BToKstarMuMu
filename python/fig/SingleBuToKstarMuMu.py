"""
Module for Single Bu To Kstar Mu Mu  

"""

__author__ = "Xin Shi <Xin.Shi@cern.ch>"
__copyright__ = "Copyright (c) Xin Shi"

import sys 
from tls import * 
import atr 
from ROOT import TFile, TH1F, TCanvas, gStyle, TPaveText
import SingleBuToKstarMuMu

def main(args, figname):
    function = getattr(SingleBuToKstarMuMu, args[0])
    figfile = function(args[1:], figname)
    sys.stdout.write('[[%s][%s.pdf]]\n' %(figfile, figname))
    return

    # return function(args[1:], figname)
    # if args[0] == 'summary':
    #     summary(args[1:], figname)
    # if args[0] == 'summary':
    #     summary(args[1:], figname)
    # else:
    #     raise NameError(args)

def summary(args, figname):
    datatype = args[0]
    label = args[1]

    test = option_exists(args, '-t')
    force = option_exists(args, '-f')

    figfile = set_file(atr.figpath, label, figname, '.root', test=test)

    if os.access(figfile, os.F_OK) and not test and not force:
        s = raw_input('File exists!: \n\
%s\nDo you want to overwrite? (yes or no) ' %figfile)
        if s == 'no':
            return
    rootfile = atr.sel.rootfile(datatype, label)
    
    cmd = './SingleBuToKstarMuMuFigures %s "%s" %s' % (label, rootfile, figfile)
    procdir = atr.sel.procdir(label)

    output = proc_cmd(cmd, procdir=procdir, shell=True) 
    print output 
    pdffile = figfile.replace('.root', '.pdf')
    return pdffile

    obj = atr.sel.root_tree_obj(datatype, label)
    chain = root_chain(rootfile, obj)

    h_mumumass = TH1F('h_mumumass', ';#mu^{+} #mu^{-} mass (GeV)', 100, 2, 4)
    h_kstarpmass = TH1F('h_kstarpmass', ';K^{*+} mass (GeV)', 100, 0.6, 1.2)
    h_kstarmmass = TH1F('h_kstarmmass', ';K^{*-} mass (GeV)', 100, 0.6, 1.2)
    h_bpmass = TH1F('h_bpmass', '; B^{+} mass (GeV)', 100, 1, 10)
    h_bmmass = TH1F('h_bmmass', '; B^{-} mass (GeV)', 100, 1, 10)
    h_bvtxcl = TH1F('h_bvtxcl', '; B Vertex CL', 100, 0, 1)
    h_blxysig = TH1F('h_blxysig', '; B Lxy significace', 100, 0, 150)
    h_bcosalphabs = TH1F('h_bcosalphabs', '; B cos#alpha beam spot', 100, -1, 1)
    h_bppt = TH1F('h_bppt', '; B^{+} pT (GeV)', 100, 0, 60)
    h_bmpt = TH1F('h_bmpt', '; B^{-} pT (GeV)', 100, 0, 60)

    h_bctau = TH1F('h_bctau', '; B c#tau (cm)', 100, 0, 1)
    h_bpmass_res = TH1F('h_bpmass_res', 'B^{+} mass resonant; B^{+} mass (GeV)', 100, 1, 10)
    h_bpmass_res_zoom_5_6 = TH1F('h_bpmass_res_zoom_5_6', 'B^{+} mass resonant; B^{+} mass (GeV)', 100, 5, 6)
    h_bpmass_nonres = TH1F('h_bpmass_nonres', 'B^{+} mass non resonant; B^{+} mass (GeV)', 100, 1, 10)
    h_bpmass_nonres_below5 = TH1F('h_bpmass_nonres_below5', 'B^{+} mass; B^{+} mass (GeV)', 100, 3, 5)

    ntot = chain.GetEntries()
    if test:
        ntot = 1000

    sys.stdout.write('Processing %s events ...\n' %ntot)
    sys.stdout.flush()
    nfill = 0 
    for i in range(ntot):
        chain.LoadTree(i)
        chain.GetEntry(i)

        h_mumumass.Fill(chain.Mumumass)
        h_bvtxcl.Fill(chain.Bvtxcl)
        h_blxysig.Fill(chain.Blxysig)
        h_bcosalphabs.Fill(chain.Bcosalphabs)
        h_bctau.Fill(chain.Bctau)

        if chain.Bchg > 0:
            h_kstarpmass.Fill(chain.Kstarmass)
            h_bpmass.Fill(chain.Bmass)
            h_bppt.Fill(chain.Bpt)

            if sel_bmass_nonres(label, chain):
                h_bpmass_nonres.Fill(chain.Bmass)
                if chain.Bmass > 3 and chain.Bmass < 5: 
                    h_bpmass_nonres_below5.Fill(chain.Bmass)

            else:
                h_bpmass_res.Fill(chain.Bmass)
                if chain.Bmass > 5 and chain.Bmass < 6: 
                    h_bpmass_res_zoom_5_6.Fill(chain.Bmass)

        else:
            h_kstarmmass.Fill(chain.Kstarmass)
            h_bmmass.Fill(chain.Bmass)
            h_bmpt.Fill(chain.Bpt)

        nfill += 1
        
    sys.stdout.write('Filled %s events. \n' %nfill)

    f = TFile(figfile, "recreate")

    h_mumumass.Write()
    h_mumumass.Delete()

    h_kstarpmass.Write()
    h_kstarpmass.Delete()

    h_kstarmmass.Write()
    h_kstarmmass.Delete()
    
    h_bpmass.Write()
    h_bpmass.Delete()

    h_bmmass.Write()
    h_bmmass.Delete()
    
    h_bppt.Write()
    h_bppt.Delete()

    h_bmpt.Write()
    h_bmpt.Delete()
    
    h_bvtxcl.Write()
    h_bvtxcl.Delete()

    h_blxysig.Write()
    h_blxysig.Delete()

    h_bcosalphabs.Write()
    h_bcosalphabs.Delete()

    h_bctau.Write()
    h_bctau.Delete()

    h_bpmass_res.Write()
    h_bpmass_res.Delete()

    h_bpmass_res_zoom_5_6.Write()
    h_bpmass_res_zoom_5_6.Delete()

    h_bpmass_nonres.Write()
    h_bpmass_nonres.Delete()

    h_bpmass_nonres_below5.Write()
    h_bpmass_nonres_below5.Delete()

    f.Close()
    sys.stdout.write('Saved as: \n %s \n' %figfile)


def mumumass(args, figname):
    datatype = args[0]
    label = args[1]
    test = option_exists(args, '-t')
    rootname = figname.replace('mumumass', 'summary')
    rootfile = set_file(atr.figpath, label, rootname, '.root', test=test)
    figfile = set_file(atr.figpath, label, figname, '.pdf', test=test)
    
    f = TFile(rootfile)
    h = f.Get('h_mumumass')

    gStyle.SetPadTopMargin(0.1) 
    c = TCanvas("aCanvas", "Canvas", 600, 600)

    ytitle = get_y_title(h, 'GeV')
    h.GetYaxis().SetTitle(ytitle)
    h.Draw()

    # Fine tune in the very later state! 
    # pt = TPaveText()
    # if label == 'Run2011v10.1': 
    #     pt = TPaveText(3.4, 800, 3.8, 1000)
    # if label == 'Run2011v10.2': 
    #     pt = TPaveText(3.4, 350, 3.8, 400)
    
    # pt.SetBorderSize(0)
    # pt.SetFillColor(0) 
    # pt.AddText("Entries = %d" %h.GetEntries());  
    # pt.Draw()
    
    c.UseCurrentStyle() 
    c.Print(figfile)
    f.Close()
    return figfile

    
def kstarpmass(args, figname):
    datatype = args[0]
    label = args[1]
    test = option_exists(args, '-t')
    #rootname = 'SingleBuToKstarMuMu_summary_data_Run2011v10'
    rootname = figname.replace('kstarpmass', 'summary')
    rootfile = set_file(atr.figpath, label, rootname, '.root', test=test)
    figfile = set_file(atr.figpath, label, figname, '.pdf', test=test)
    c = TCanvas("aCanvas", "Canvas", 600, 600)
    gStyle.SetPadLeftMargin(0.15)
    gStyle.SetOptTitle(0) 
    c.UseCurrentStyle() 
    
    f = TFile(rootfile)
    h = f.Get('h_kstarpmass')
    h.GetYaxis().SetTitle(get_y_title(h, 'GeV'))
    h.GetYaxis().SetTitleOffset(1.5)
    h.Draw()
    
    pt = TPaveText()
    if label == 'Run2011v10.1': 
        pt = TPaveText(1, 70, 1.1, 80)
    if label == 'Run2011v10.2': 
        pt = TPaveText(1, 70, 1.1, 80)

    pt.SetBorderSize(0)
    pt.SetFillColor(0) 
    pt.AddText("Entries = %d" %h.GetEntries());  
    pt.Draw()

    c.Print(figfile)
    f.Close()
    
    
def bpmass(args, figname):
    datatype = args[0]
    label = args[1]
    test = option_exists(args, '-t')
    #rootname = 'SingleBuToKstarMuMu_summary_data_Run2011v10'
    rootname = figname.replace('bpmass', 'summary')
    rootfile = set_file(atr.figpath, label, rootname, '.root', test=test)
    figfile = set_file(atr.figpath, label, figname, '.pdf', test=test)
    c = TCanvas("aCanvas", "Canvas", 600, 600)
    gStyle.SetPadLeftMargin(0.15)
    gStyle.SetOptTitle(0) 

    
    f = TFile(rootfile)
    h = f.Get('h_bpmass')
    h.GetYaxis().SetTitle(get_y_title(h, 'GeV'))
    h.GetYaxis().SetTitleOffset(1.6)
    h.Draw()
    
    # pt = TPaveText()
    # if label == 'Run2011v10.1': 
    #     pt = TPaveText(7, 200, 9, 220)
    # if label == 'Run2011v10.2': 
    #     pt = TPaveText(7, 80, 9, 90)

    # pt.SetBorderSize(0)
    # pt.SetFillColor(0) 
    # pt.AddText("Entries = %d" %h.GetEntries());  
    # pt.Draw()
    c.UseCurrentStyle() 

    c.Print(figfile)
    f.Close()
    return figfile

    

def bvtxcl(args, figname):
    datatype = args[0]
    label = args[1]
    test = option_exists(args, '-t')
    rootname = figname.replace('bvtxcl', 'summary')
    rootfile = set_file(atr.figpath, label, rootname, '.root', test=test)
    figfile = set_file(atr.figpath, label, figname, '.pdf', test=test)
    c = TCanvas("aCanvas", "Canvas", 600, 600)
    gStyle.SetPadLeftMargin(0.15)
    gStyle.SetOptTitle(0) 
    c.UseCurrentStyle() 
    
    f = TFile(rootfile)
    h = f.Get('h_bvtxcl')
    h.GetYaxis().SetTitle(get_y_title(h, ''))
    h.GetYaxis().SetTitleOffset(1.6)
    h.Draw()
    
    pt = TPaveText()
    pt.SetBorderSize(0)
    pt.SetFillColor(0) 
    pt.AddText("Entries = %d" %h.GetEntries());  
    pt.Draw()

    c.Print(figfile)
    f.Close()
    
def blxysig(args, figname):
    datatype = args[0]
    label = args[1]
    test = option_exists(args, '-t')
    rootname = figname.replace('blxysig', 'summary')
    rootfile = set_file(atr.figpath, label, rootname, '.root', test=test)
    figfile = set_file(atr.figpath, label, figname, '.pdf', test=test)
    c = TCanvas("aCanvas", "Canvas", 600, 600)
    gStyle.SetPadLeftMargin(0.15)
    gStyle.SetOptTitle(0) 
    c.UseCurrentStyle() 
    
    f = TFile(rootfile)
    h = f.Get('h_blxysig')
    h.GetYaxis().SetTitle(get_y_title(h, ''))
    h.GetYaxis().SetTitleOffset(1.6)
    h.Draw()
    
    pt = TPaveText(1.4, 5e5, 1.8, 6e5)
    pt.SetBorderSize(0)
    pt.SetFillColor(0) 
    pt.AddText("Entries = %d" %h.GetEntries());  
    pt.Draw()

    c.Print(figfile)
    f.Close()
    
    
def bcosalphabs(args, figname):
    datatype = args[0]
    label = args[1]
    test = option_exists(args, '-t')
    rootname = figname.replace('bcosalphabs', 'summary')
    rootfile = set_file(atr.figpath, label, rootname, '.root', test=test)
    figfile = set_file(atr.figpath, label, figname, '.pdf', test=test)
    c = TCanvas("aCanvas", "Canvas", 600, 600)
    gStyle.SetPadLeftMargin(0.15)
    gStyle.SetOptTitle(0) 
    c.UseCurrentStyle() 
    
    f = TFile(rootfile)
    h = f.Get('h_bcosalphabs')
    h.GetYaxis().SetTitle(get_y_title(h, ''))
    h.GetYaxis().SetTitleOffset(1.6)
    h.Draw()
    
    pt = TPaveText(0.2, 1.6e5, 0.6, 1.8e5)
    pt.SetBorderSize(0)
    pt.SetFillColor(0) 
    pt.AddText("Entries = %d" %h.GetEntries());  
    pt.Draw()

    c.Print(figfile)
    f.Close()
    
   
def bctau(args, figname):
    datatype = args[0]
    label = args[1]
    test = option_exists(args, '-t')
    rootname = figname.replace('bctau', 'summary')
    rootfile = set_file(atr.figpath, label, rootname, '.root', test=test)
    figfile = set_file(atr.figpath, label, figname, '.pdf', test=test)
    c = TCanvas("aCanvas", "Canvas", 600, 600)
    gStyle.SetPadLeftMargin(0.15)
    gStyle.SetOptTitle(0) 
    c.UseCurrentStyle() 
    
    f = TFile(rootfile)
    h = f.Get('h_bctau')
    h.GetYaxis().SetTitle(get_y_title(h, 'cm'))
    h.GetYaxis().SetTitleOffset(1.6)
    h.Draw()
    
    pt = TPaveText(0.2, 1.6e5, 0.6, 1.8e5)
    pt.SetBorderSize(0)
    pt.SetFillColor(0) 
    pt.AddText("Entries = %d" %h.GetEntries());  
    pt.Draw()

    c.Print(figfile)
    f.Close()
    
def bppt(args, figname):
    datatype = args[0]
    label = args[1]
    test = option_exists(args, '-t')
    rootname = figname.replace('bppt', 'summary')
    rootfile = set_file(atr.figpath, label, rootname, '.root', test=test)
    figfile = set_file(atr.figpath, label, figname, '.pdf', test=test)
    c = TCanvas("aCanvas", "Canvas", 600, 600)
    gStyle.SetPadLeftMargin(0.15)
    gStyle.SetOptTitle(0) 
    c.UseCurrentStyle() 
    
    f = TFile(rootfile)
    h = f.Get('h_bppt')
    h.GetYaxis().SetTitle(get_y_title(h, 'GeV'))
    h.GetYaxis().SetTitleOffset(1.6)
    h.Draw()
    
    pt = TPaveText(40, 1.6e5, 50, 1.8e5)
    pt.SetBorderSize(0)
    pt.SetFillColor(0) 
    pt.AddText("Entries = %d" %h.GetEntries());  
    pt.Draw()

    c.Print(figfile)
    f.Close()
    
   
def sel_bmass_nonres(label, chain):
    jpsimass = 3.09692 
    psi2smass = 3.68609 

    if ( label == 'Run2011v10.2' or 
         label == 'Run2011v10.3' or
         label == 'Run2011v10.4' or
         label == 'Run2011v10.5') : 
        jpsimass_min = jpsimass - 5*0.03 # 2.94692
        jpsimass_max = jpsimass + 5*0.03 # 3.24692
        
        psi2smass_min = psi2smass - 5*0.03 # 3.53609 
        psi2smass_max = psi2smass + 5*0.03 # 3.83609

        mumumass = chain.Mumumass
        
        if ( (mumumass < jpsimass_min )
             or (mumumass > jpsimass_max and mumumass < psi2smass_min)
             or (mumumass > psi2smass_max)):
            return True

    else:
        raise NameError(label)

    return False


def bpmass_res(args, figname):
    datatype = args[0]
    label = args[1]
    test = option_exists(args, '-t')
    rootname = figname.replace('bpmass_res', 'summary')
    rootfile = set_file(atr.figpath, label, rootname, '.root', test=test)
    figfile = set_file(atr.figpath, label, figname, '.pdf', test=test)
    c = TCanvas("aCanvas", "Canvas", 600, 600)
    gStyle.SetPadLeftMargin(0.15)
    gStyle.SetOptTitle(0) 
    c.UseCurrentStyle() 
    
    f = TFile(rootfile)
    h = f.Get('h_bpmass_res')
    h.GetYaxis().SetTitle(get_y_title(h, 'GeV'))
    h.GetYaxis().SetTitleOffset(1.6)
    h.Draw()
    
    # pt = TPaveText()
    # if label == 'Run2011v10.2': 
    #     pt = TPaveText(7, 30, 9, 35)

    # pt.SetBorderSize(0)
    # pt.SetFillColor(0) 
    # pt.AddText("Entries = %d" %h.GetEntries());  
    # pt.Draw()

    c.UseCurrentStyle() 
    c.Print(figfile)
    f.Close()
    return figfile


def bpmass_res_zoom_5_6(args, figname):
    datatype = args[0]
    label = args[1]
    test = option_exists(args, '-t')
    rootname = figname.replace('bpmass_res_zoom_5_6', 'summary')
    rootfile = set_file(atr.figpath, label, rootname, '.root', test=test)
    figfile = set_file(atr.figpath, label, figname, '.pdf', test=test)
    c = TCanvas("aCanvas", "Canvas", 600, 600)
    gStyle.SetPadLeftMargin(0.15)
    gStyle.SetOptTitle(0) 
    c.UseCurrentStyle() 
    
    f = TFile(rootfile)
    h = f.Get('h_bpmass_res_zoom_5_6')
    h.GetYaxis().SetTitle(get_y_title(h, 'GeV'))
    h.GetYaxis().SetTitleOffset(1.6)
    h.Draw()
    
    # pt = TPaveText()
    # if label == 'Run2011v10.2': 
    #     pt = TPaveText(5.6, 10, 5.9, 12)

    # pt.SetBorderSize(0)
    # pt.SetFillColor(0) 
    # pt.AddText("Entries = %d" %h.GetEntries());  
    # pt.Draw()

    c.UseCurrentStyle() 
    c.Print(figfile)
    f.Close()
    return figfile
    

def bpmass_nonres(args, figname):
    datatype = args[0]
    label = args[1]
    test = option_exists(args, '-t')
    rootname = figname.replace('bpmass_nonres', 'summary')
    rootfile = set_file(atr.figpath, label, rootname, '.root', test=test)
    figfile = set_file(atr.figpath, label, figname, '.pdf', test=test)
    c = TCanvas("aCanvas", "Canvas", 600, 600)
    gStyle.SetPadLeftMargin(0.15)
    gStyle.SetOptTitle(0) 
    c.UseCurrentStyle() 
    
    f = TFile(rootfile)
    h = f.Get('h_bpmass_nonres')
    h.GetYaxis().SetTitle(get_y_title(h, 'GeV'))
    h.GetYaxis().SetTitleOffset(1.6)
    h.Draw()
    
    # pt = TPaveText()
    # if label == 'Run2011v10.2': 
    #     pt = TPaveText(7, 30, 9, 35)

    # pt.SetBorderSize(0)
    # pt.SetFillColor(0) 
    # pt.AddText("Entries = %d" %h.GetEntries());  
    # pt.Draw()

    c.UseCurrentStyle() 
    c.Print(figfile)
    f.Close()
    return figfile

    
def bpmass_nonres_below5(args, figname):
    datatype = args[0]
    label = args[1]
    test = option_exists(args, '-t')
    rootname = figname.replace('bpmass_nonres_below5', 'summary')
    rootfile = set_file(atr.figpath, label, rootname, '.root', test=test)
    figfile = set_file(atr.figpath, label, figname, '.pdf', test=test)
    c = TCanvas("aCanvas", "Canvas", 600, 600)
    gStyle.SetPadLeftMargin(0.15)
    gStyle.SetOptTitle(0) 
    c.UseCurrentStyle() 
    
    f = TFile(rootfile)
    h = f.Get('h_bpmass_nonres_below5')
    h.GetYaxis().SetTitle(get_y_title(h, 'GeV'))
    h.GetYaxis().SetTitleOffset(1.6)
    h.Draw()
    
    pt = TPaveText()
    if label == 'Run2011v10.2': 
        pt = TPaveText(4.2, 10, 4.8, 12)

    pt.SetBorderSize(0)
    pt.SetFillColor(0) 
    pt.AddText("Entries = %d" %h.GetEntries());  
    pt.Draw()

    c.Print(figfile)
    f.Close()

    
    
   

    
    

    
    

    
