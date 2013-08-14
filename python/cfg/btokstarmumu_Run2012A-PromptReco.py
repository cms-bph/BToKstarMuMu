import FWCore.ParameterSet.Config as cms
from btokstarmumu_2012_cfi import process 

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('file:Run2012B.root' )
                            #fileNames = cms.untracked.vstring('/store/data/Run2012B/ParkingMonitor/AOD/PromptReco-v1/000/194/305/B2F35FE7-B5A1-E111-9B27-E0CB4E55367F.root')
                            
    )
    #process.GlobalTag.globaltag = cms.string('GR_R_42_V25::All')
process.GlobalTag.globaltag = cms.string('GR_R_52_V9D::All')
   
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
    ('HLT_DoubleMu3p5_LowMass_Displaced',  'hltDisplacedmumuFilterDoubleMu3p5LowMass') 
    ]

g_TriggerNames = [i[0] for i in g_TriggerNames_LastFilterNames]
g_LastFilterNames = [i[1] for i in g_TriggerNames_LastFilterNames]

process.ntuple.TriggerNames = cms.vstring(g_TriggerNames)
process.ntuple.LastFilterNames = cms.vstring(g_LastFilterNames)

test = False 
if test: 
    process.ntuple.FileName = cms.string("test.root")
    process.ntuple.KstarMinMass = cms.untracked.double(0.2)
    process.ntuple.KstarMaxMass = cms.untracked.double(100.29)
    process.ntuple.BMaxMass = cms.untracked.double(100.0)
    process.ntuple.MuonMinPt = cms.untracked.double(0.0) # [GeV]
    process.ntuple.MuonMaxEta = cms.untracked.double(5)  
    process.ntuple.TrkMaxDcaBs = cms.untracked.double(10.0) # [cm]
    process.ntuple.TrkMaxR = cms.untracked.double(1110.0) # [cm]
    process.ntuple.TrkMaxZ = cms.untracked.double(1280.0) # [cm]
    process.ntuple.MuMuMaxDca = cms.untracked.double(10.5) # [cm]
    process.ntuple.MuMuMinVtxCl = cms.untracked.double(10.05) 
    process.ntuple.MuMuMinPt = cms.untracked.double(16.9) # [GeV/c]
    process.ntuple.MuMuMinInvMass = cms.untracked.double(111.0) # [GeV/c2]
    process.ntuple.MuMuMaxInvMass = cms.untracked.double(114.8) # [GeV/c2]
    process.ntuple.MuMuMinLxySigmaBs = cms.untracked.double(13.0) 
    process.ntuple.MuMuMinCosAlphaBs = cms.untracked.double(0.0)

    
