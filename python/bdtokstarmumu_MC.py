import FWCore.ParameterSet.Config as cms
from btokstarmumu_cfi import process 

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                                'file:/afs/cern.ch/user/x/xshi/work/cms/afb/dat/aod/mc/BuToKstarJPsi_7TeV_PYTHIA6_100.root', 
                            )
                        )
process.GlobalTag.globaltag = cms.string('START42_V14A::All')

# do trigger matching for muons
triggerProcessName = 'HLT'

process.cleanMuonTriggerMatchHLT0 = cms.EDProducer(
    # match by DeltaR only (best match by DeltaR)
    'PATTriggerMatcherDRLessByR',                         
    src                   = cms.InputTag('cleanPatMuons'),
    # default producer label as defined in
    # PhysicsTools/PatAlgos/python/triggerLayer1/triggerProducer_cfi.py
    matched               = cms.InputTag('patTrigger'),
    matchedCuts           = cms.string('path("HLT_Dimuon6p5_LowMass_Displaced_v*",0,0)'),
    maxDeltaR             = cms.double(0.1),
    # only one match per trigger object
    resolveAmbiguities    = cms.bool(True),
    # take best match found per reco object (by DeltaR here, see above)       
    resolveByMatchQuality = cms.bool(False))

from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTriggerMatchEmbedding(process, triggerMatchers = ['cleanMuonTriggerMatchHLT0'],
                              hltProcess = triggerProcessName, outputModule = '')

g_TriggerNames_LastFilterNames = [
    ('HLT_Dimuon6p5_LowMass_Displaced',  'hltDisplacedmumuFilterLowMass')  #5E32, v8.1-v8.3
    ]

g_TriggerNames = [i[0] for i in g_TriggerNames_LastFilterNames]
g_LastFilterNames = [i[1] for i in g_TriggerNames_LastFilterNames]

process.ntuple.BuildBuToKstarMuMu = cms.untracked.bool(False)
process.ntuple.BuildBdToKstarMuMu = cms.untracked.bool(True)

process.ntuple.TriggerNames = cms.vstring(g_TriggerNames)
process.ntuple.LastFilterNames = cms.vstring(g_LastFilterNames)
process.ntuple.IsMonteCarlo = cms.untracked.bool(True)

#process.p = cms.Path(process.patDefaultSequence * process.ntuple)

