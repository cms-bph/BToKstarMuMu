import FWCore.ParameterSet.Config as cms
from btokstarmumu_cfi import process 

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2000))

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    '/store/data/Run2011A/MuOnia/AOD/12Oct2013-v1/00000/002104A5-2A44-E311-B69F-00304867918A.root')
    )
process.GlobalTag.globaltag = cms.string('FT_53_LV5_AN1::All')  # GT for 2011AB legacy ReReco in cmssw_5_3_X

print "\n=> Global Tag: FT_53_LV5_AN1::All \n"

# do trigger matching for muons
triggerProcessName = 'HLT'


process.cleanMuonTriggerMatchHLT1 = cms.EDProducer(
    'PATTriggerMatcherDRLessByR',
    src                   = cms.InputTag('cleanPatMuons'),
    # default producer label as defined in 
    # PhysicsTools/PatAlgos/python/triggerLayer1/triggerProducer_cfi.py
    matched               = cms.InputTag('patTrigger'),
    matchedCuts           = cms.string('path("HLT_Dimuon7_LowMass_Displaced_v*")'),
    maxDeltaR             = cms.double(0.1),
    # only one match per trigger object
    resolveAmbiguities    = cms.bool(True),
    # take best match found per reco object (by DeltaR here, see above)
    resolveByMatchQuality = cms.bool(False)) #####

process.cleanMuonTriggerMatchHLT2 = cms.EDProducer(
    'PATTriggerMatcherDRLessByR',
    src                   = cms.InputTag('cleanPatMuons'),
    matched               = cms.InputTag('patTrigger'),
    matchedCuts           = cms.string('path("HLT_DoubleMu4_LowMass_Displaced_v*")'),
    maxDeltaR             = cms.double(0.1),
    resolveAmbiguities    = cms.bool(True),
    resolveByMatchQuality = cms.bool(False)) #####

process.cleanMuonTriggerMatchHLT3 = cms.EDProducer(
    'PATTriggerMatcherDRLessByR',
    src                   = cms.InputTag('cleanPatMuons'),
    matched               = cms.InputTag('patTrigger'),
    matchedCuts           = cms.string('path("HLT_DoubleMu4p5_LowMass_Displaced_v*")'),
    maxDeltaR             = cms.double(0.1),
    resolveAmbiguities    = cms.bool(True),
    resolveByMatchQuality = cms.bool(False)) #####

process.cleanMuonTriggerMatchHLT4 = cms.EDProducer(
    'PATTriggerMatcherDRLessByR',
    src                   = cms.InputTag('cleanPatMuons'),
    matched               = cms.InputTag('patTrigger'),
    matchedCuts           = cms.string('path("HLT_DoubleMu5_LowMass_Displaced_v*")'),
    maxDeltaR             = cms.double(0.1),
    resolveAmbiguities    = cms.bool(True),
    resolveByMatchQuality = cms.bool(False)) #####

from PhysicsTools.PatAlgos.tools.trigTools import *

switchOnTriggerMatchEmbedding(process, triggerMatchers = ['cleanMuonTriggerMatchHLT1','cleanMuonTriggerMatchHLT2','cleanMuonTriggerMatchHLT3','cleanMuonTriggerMatchHLT4'], hltProcess = triggerProcessName, outputModule = '')

g_TriggerNames_LastFilterNames = [
    ('HLT_Dimuon7_LowMass_Displaced',
      'hltDisplacedmumuFilterLowMass'), #1E33, 1.4E33
      
    ('HLT_DoubleMu4_LowMass_Displaced',
     'hltDisplacedmumuFilterLowMass'), #2E33

    ('HLT_DoubleMu4p5_LowMass_Displaced',
     'hltDisplacedmumuFilterDoubleMu4p5LowMass'), #3E33, 5E33
     
    ('HLT_DoubleMu5_LowMass_Displaced',
     'hltDisplacedmumuFilterDoubleMu5LowMass'), #3E33, 5E33
    ]

g_TriggerNames = [i[0] for i in g_TriggerNames_LastFilterNames]
g_LastFilterNames = [i[1] for i in g_TriggerNames_LastFilterNames]

process.ntuple.TriggerNames = cms.vstring(g_TriggerNames)
process.ntuple.LastFilterNames = cms.vstring(g_LastFilterNames)

