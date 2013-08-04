"""
Module for Download 

"""

__author__ = "Xin Shi <Xin.Shi@cern.ch>"
__copyright__ = "Copyright (c) Xin Shi"

import sys 
from tls import * 

def main(args):
    if args[0] == 'fit':
        dld_fit(args[1:])
    else:
        raise NameError(args)

def dld_fit(args):
    var = args[0]
    datatype = args[1]
    label = args[2]
    figname = 'fit_%s' % var
    
    test = get_options(args, 'test')    
    if test:
        figname += '_test'
        
    rpdffile = set_figfile(figname, label, '.pdf', remote=True)
    pdffile = set_figfile(figname, label, '.pdf')

    rtxtfile = set_figfile(figname, label, '.txt', remote=True)
    txtfile = set_figfile(figname, label, '.txt')

    check_and_copy(rpdffile, pdffile)
    check_and_copy(rtxtfile, txtfile)


