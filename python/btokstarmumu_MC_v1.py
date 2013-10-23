import FWCore.ParameterSet.Config as cms
from btokstarmumu_cfi import process 

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )


process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
'file:/afs/cern.ch/user/x/xshi/work/cms/afb/dat/aod/mc/BuToKstarMuMu_7TeV_PYTHIA6.root', 
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

process.ntuple.TriggerNames = cms.vstring(g_TriggerNames)
process.ntuple.LastFilterNames = cms.vstring(g_LastFilterNames)
process.ntuple.SaveGenInfo = cms.untracked.bool(True)

# use very loose cuts 
process.ntuple.MuonMinPt = cms.untracked.double(0.) # 3.0 
process.ntuple.MuonMinPt = cms.untracked.double(0.) #3.0), # [GeV]
process.ntuple.MuonMaxEta = cms.untracked.double(10) # 2.2),  
process.ntuple.TrkMaxDcaBs = cms.untracked.double(100) #2.0), # [cm]
process.ntuple.TrkMaxR = cms.untracked.double(300) #110.0), # [cm]
process.ntuple.TrkMaxZ = cms.untracked.double(900) # 280.0), # [cm]
process.ntuple.MuMuMaxDca = cms.untracked.double(10) # 0.5), # [cm]
process.ntuple.MuMuMinVtxCl = cms.untracked.double(0) # 0.05), 
process.ntuple.MuMuMinPt = cms.untracked.double(0) # 6.9), # [GeV/c]
process.ntuple.MuMuMinInvMass = cms.untracked.double(0) #1.0), # [GeV/c2]
process.ntuple.MuMuMaxInvMass = cms.untracked.double(10) #4.8), # [GeV/c2]
process.ntuple.MuMuMinLxySigmaBs = cms.untracked.double(0) #3.0), 
process.ntuple.MuMuMinCosAlphaBs = cms.untracked.double(0) #0.9),

process.ntuple.KstarChargedTrackMinPt = cms.untracked.double(0.0) # no change, # [GeV/c]
process.ntuple.KstarMinMass = cms.untracked.double(0) # 0.49), # [GeV/c2] K*+ mass = 891.66 +- 0.26 MeV 
process.ntuple.KstarMaxMass = cms.untracked.double(50) #1.29), # [GeV/c2] K*0 mass = 895.94 +- 0.26 MeV

process.ntuple.BMaxMass = cms.untracked.double(100 ) #8.0), # [GeV/c2]n B+ mass = 5279 MeV 


    
