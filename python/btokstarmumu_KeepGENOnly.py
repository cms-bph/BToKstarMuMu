################
#  variables
################
runMC               = True
run2012not2011      = True  # True for 2012 and False for 2011

print "\n@@@ CMSSW run configuration flags @@@"

if (runMC == True):
    print " running MC : ", runMC
else:
    print " not running MC, please xCheck !!"

if (run2012not2011 == True):
    print "---> 2012 MC running"
else:
    print "---> 2011 MC running"

print "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"

##################
# cmssw configs
##################

import FWCore.ParameterSet.Config as cms
from btokstarmumu_cfi import process 

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
'file:/wk2/pchen/work/BToKstarMuMu/localResources/dpm/BuToKstarMuMu_NoFilter_8TeV_GENSIM_v1/BuToKstarMuMu_NoFilter_8TeV_GENSIM_v1/fb03d90927ce761d45f353f0dd279fce/BuToKstarMuMu_NoFilter_8TeV-pythia6-evtgen_GENSIM_100_1_iUG.root',
                            )
                        )
if (run2012not2011 == True):
    process.GlobalTag.globaltag = cms.string('START53_V19F::All')
    print "\nGlobalTag : START53_V19F::All\n"
else:
    process.GlobalTag.globaltag = cms.string('START42_V14B::All')
    print "\nGlobalTag : START42_V14B::All\n"

process.ntuple.IsMonteCarlo = cms.untracked.bool(True)
process.ntuple.KeepGENOnly  = cms.untracked.bool(True)
process.p = cms.Path(process.ntuple)
 
