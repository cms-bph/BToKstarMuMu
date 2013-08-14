"""
Attributes for the A_FB analysis.

"""

__author__ = "Xin Shi <Xin.Shi@cern.ch>"
__copyright__ = "Copyright (c) Xin Shi"

import os
import sys

import dat
import ntp
import sel
import srm
import cfg

afbpath = os.environ['afb']
datpath = os.path.join(afbpath, 'dat')
figpath = os.path.join(afbpath, 'doc/fig')
logpath = os.path.join(afbpath, 'log')
bashpath = os.path.join(afbpath, 'src/sh')

castor_afbpath = '/castor/cern.ch/user/x/xshi/afb'
castor_datpath = os.path.join(castor_afbpath, 'dat')

rafbpath = os.environ['rafb']
rfigpath = os.path.join(rafbpath, 'doc', 'fig')

muon_mass = 0.10565837 # GeV 

def get_afb_from_label(label):
    afb = 0 
    return afb 

