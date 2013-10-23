"""
Module for B0 Mass fitting 

"""

__author__ = "Xin Shi <Xin.Shi@cern.ch>"
__copyright__ = "Copyright (c) Xin Shi"

import sys 
from tls import *
from ROOT import gROOT, kTRUE
import atr

def main(args):
    gROOT.SetBatch(kTRUE)
    name = args[0]
    datatype = args[1]
    label = args[2]

    test = get_options(args, 'test')
    batch = get_options(args, 'batch')

    #selfile = atr.sel.rootfile(datatype, label, name)
    sel_label = get_label_by_version(label, 'x.x')
    selfile = atr.sel.rootfile(datatype, sel_label, name)

    ccname = 'B0MassFitter.cc'
    ccpath = os.path.join(os.environ['HOME'], 'bak/src/cms/afb/cc')
    ccfile = os.path.join(ccpath, ccname)

    jpsimass = 3.09692 
    psi2smass = 3.68609 

    title = 'B^{0}#rightarrowK*#mu^{+}#mu^{-}'
    datacut = 'mumumass<(%s-5*mumumasserr) || (mumumass>(%s+3*mumumasserr) &&\
mumumass<(%s-3*mumumasserr)) || mumumass>(%s+3*mumumasserr)' %(
    jpsimass, jpsimass, psi2smass, psi2smass)
    
    figname = 'b0mass_kstmumu'
    if '-jpsi' in args: 
        title = 'B^{0}#rightarrowK*J/#psi'
        datacut = 'mumumass > (%s-5*mumumasserr) && mumumass < (%s+3*mumumasserr)' %(
            jpsimass, jpsimass)

        figname = 'b0mass_kstjpsi'

    pdffile = set_file(atr.figpath, label, figname, '.pdf', test=test)
    resfile = set_file(atr.figpath, label, figname, '.root', test=test)

    par_str = '( "%s", "%s", "%s", "%s", "%s")'
    par_tuple = (title, selfile, datacut, pdffile, resfile)
    print_paras(par_tuple)
    par = par_str % par_tuple

    sys.stdout.write('  ROOT Macro = %s \n' % ccfile)
    if test:
        return
    gROOT.Macro( ccfile + par)


def print_paras(paras):
    sys.stdout.write( '  title   = %s\n'     %paras[0])
    sys.stdout.write( '  selfile = %s\n'     %paras[1])
    sys.stdout.write( '  datacut = %s\n'     %paras[2])
    sys.stdout.write( '  pdffile = %s\n'     %paras[3])
    sys.stdout.write( '  resfile = %s\n'     %paras[4])

    

