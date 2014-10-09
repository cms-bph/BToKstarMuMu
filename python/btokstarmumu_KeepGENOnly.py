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
        fileNames = cms.untracked.vstring('root://eoscms//store/user/pchen/BToKstarMuMu/dat/GEN/PYTHIA6_BuToKstarMuMu_TuneZ2star_8TeV_GEN_NoFilter_999_2_3Yn.root' )
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
