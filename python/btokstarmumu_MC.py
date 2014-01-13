import FWCore.ParameterSet.Config as cms
from btokstarmumu_cfi import process 

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                               #'file:/afs/cern.ch/user/x/xshi/work/cms/afb/dat/aod/mc/BuToKstarMuMu_7TeV_PYTHIA6.root', 
# 'file:/afs/cern.ch/user/x/xshi/work/cms/afb/dat/aod/mc/BuToKstarMuMu_7TeV_PYTHIA6_1_1_96m.root', 
# 'file:/afs/cern.ch/user/x/xshi/work/cms/afb/dat/aod/mc/BuToKstarMuMu_7TeV_PYTHIA6_2_1_dfI.root', 
# 'file:/afs/cern.ch/user/x/xshi/work/cms/afb/dat/aod/mc/BuToKstarMuMu_7TeV_PYTHIA6_3_1_ud5.root', 
# 'file:/afs/cern.ch/user/x/xshi/work/cms/afb/dat/aod/mc/BuToKstarMuMu_7TeV_PYTHIA6_4_1_neu.root', 
# 'file:/afs/cern.ch/user/x/xshi/work/cms/afb/dat/aod/mc/BuToKstarMuMu_7TeV_PYTHIA6_5_1_8Kc.root', 
#                                'file:/afs/cern.ch/work/n/nsahoo/BPH-ANALYSIS/afb/SE/aod/mc/2012_test_BuToKstarJPsi_500.root',
#                                'file:/afs/cern.ch/work/n/nsahoo/BPH-ANALYSIS/afb/SE/aod/mc/2012_test_BuToKstarMuMu_500.root',
                                'file:/afs/cern.ch/work/n/nsahoo/BPH-ANALYSIS/afb/SE/aod/mc/2012_test_BuToKstarPsi2S_500.root', 
                            )
                        )
#process.GlobalTag.globaltag = cms.string('GR_R_42_V25::All')r
#process.GlobalTag.globaltag = cms.string('START42_V14A::All')
process.GlobalTag.globaltag = cms.string('START53_V19F::All')

# do trigger matching for muons
triggerProcessName = 'HLT'

process.cleanMuonTriggerMatchHLT0 = cms.EDProducer(
    # match by DeltaR only (best match by DeltaR)
    'PATTriggerMatcherDRLessByR',                         
    src                   = cms.InputTag('cleanPatMuons'),
    # default producer label as defined in
    # PhysicsTools/PatAlgos/python/triggerLayer1/triggerProducer_cfi.py
    matched               = cms.InputTag('patTrigger'),
    matchedCuts           = cms.string('path("HLT_DoubleMu3p5_LowMass_Displaced*",0,0)'),
    maxDeltaR             = cms.double(0.1),
    # only one match per trigger object
    resolveAmbiguities    = cms.bool(True),
    # take best match found per reco object (by DeltaR here, see above)       
    resolveByMatchQuality = cms.bool(False))

from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTriggerMatchEmbedding(process, triggerMatchers = ['cleanMuonTriggerMatchHLT0'],
                              hltProcess = triggerProcessName, outputModule = '')

g_TriggerNames_LastFilterNames = [
    ('HLT_DoubleMu3p5_LowMass_Displaced',  'hltDisplacedmumuFilterDoubleMu3p5LowMass')  #5E32, v8.1-v8.3
    ]

g_TriggerNames = [i[0] for i in g_TriggerNames_LastFilterNames]
g_LastFilterNames = [i[1] for i in g_TriggerNames_LastFilterNames]

process.ntuple.TriggerNames = cms.vstring(g_TriggerNames)
process.ntuple.LastFilterNames = cms.vstring(g_LastFilterNames)
process.ntuple.IsMonteCarlo = cms.untracked.bool(True)
 
