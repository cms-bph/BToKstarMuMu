"""
Module for SingleCandB0KstMuMu Selection

"""

__author__ = "Xin Shi <Xin.Shi@cern.ch>"
__copyright__ = "Copyright (c) Xin Shi"

import sys
import atr
from tls import *


def main(args):
    datatype = args[0]
    label = args[1]
    test = get_options(args, 'test')

    ntp_label = get_label_by_version(label, 'x.x')
    ntp_label_name = get_name_from_label(ntp_label)
    ntppath = os.path.join(atr.datpath, 'ntp', datatype, ntp_label_name)

    inputname = 'B0ToKstMuMu.root'
    output_version = label.split('.')[-1]
    outputname = 'SingleCandB0KstMuMu_v%s.root' % output_version

    infile = os.path.join(ntppath, inputname)
    outfile = os.path.join(ntppath, outputname)

    binname = 'SingleCandB0KstMuMu'
    
    if label in ['Run2012B-PromptReco-v1.5.1',
                 'Run2012B-PromptReco-v1.5.2',
                 ]:
        cmssw = 'CMSSW_5_3_2_patch4'
        procdir = os.path.join(atr.afbpath, 'rel', cmssw,
                               'src/BphAna/B0KstMuMu_V00_06_09/plugins')
    else:
        raise NameError(label)

    cmd = './%s singlecand %s %s' %(binname, infile, outfile)
    
    output = proc_cmd(cmd, procdir=procdir, test=test)

    print output 
