"""
Module for Merging ROOT files for BuToKstarMuMu analysis. 

"""

__author__ = "Xin Shi <Xin.Shi@cern.ch>"
__copyright__ = "Copyright (c) Xin Shi"

import sys 
from tls import * 
import atr
import SingleBuToKstarMuMu 
import glob

def main(args):
    datatype = args[0]
    label = args[1]

    #com_name = get_name_from_label(label)

    run_label = label.split('/')[0]
    base = os.path.join(atr.datpath, 'sel', datatype, run_label)
    ntp_labels = atr.sel.ntp_labels(label)
    for ntp_label in ntp_labels:
        subdir = ntp_label.split('/')[1]

        dstfile = os.path.join(base, '%s.root' % subdir)

        srcfiles = ' '.join(glob.glob('%s/%s/SingleBuToKstarMuMu_*.root'
                                      % (base, subdir)))
        
        cmd = 'hadd %s %s' %(dstfile, srcfiles)
        if option_exists(args, '-f'): 
            cmd = 'hadd -f %s %s' %(dstfile, srcfiles)

        if option_exists(args, '-t'): 
            print cmd
            continue 

        output = proc_cmd(cmd)
        print output 
        
    #function = getattr(SingleBuToKstarMuMu, com_name)
    #return function(args)


def Run2011v10_1(args):
    subdirs = ['B_PromptReco_v1_10_1.1', 
               'PromptReco_v4_10.1', 
               'PromptReco_v6_10.1',
               'May10ReReco_v1.1',
               'PromptReco_v5_10.1']
    
    base = os.path.join(atr.datpath, 'sel/data/Run2011v10.1')
    for subdir in subdirs:
        dstfile = os.path.join(base, '%s.root' % subdir)

        srcfiles = ' '.join(glob.glob('%s/%s/SingleBuToKstarMuMu_*.root'
                                      % (base, subdir)))
        
        cmd = 'hadd %s %s' %(dstfile, srcfiles)
        if option_exists(args, '-f'): 
            cmd = 'hadd -f %s %s' %(dstfile, srcfiles)

        if option_exists(args, '-t'): 
            print cmd
            return 

        output = proc_cmd(cmd)
        print output 
        


