"""
Configuration File Initializer for B2KstarMuMu analysis

"""

import FWCore.ParameterSet.Config as cms
import sys 
from PhysicsTools.PatAlgos.tools.trackTools import makeTrackCandidates
from b2kstarmumu_cfi import HLTBitNames_DoubleMu_v1, HLTLastFilterNames_DoubleMu_v1, MuonSelection_v1
from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import addMCinfo, useExistingPATMuons, useL1MatchingWindowForSinglets, changeTriggerProcessName, switchOffAmbiguityResolution, addDiMuonTriggers

process = cms.Process("NTUPLE")

# Input source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('file:RECO.root' ))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100))


# initial configurations
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.MessageLogger.cerr.FwkReport.reportEvery = 10

# set global environment
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'START42_V14A::All'


# load default PAT process
process.load("PhysicsTools.PatAlgos.patSequences_cff")


# discard overlaps for muons and electrons
process.load("PhysicsTools.PatAlgos.cleaningLayer1.genericTrackCleaner_cfi")
process.cleanPatTracks.checkOverlaps.muons.requireNoOverlaps = cms.bool(False)
process.cleanPatTracks.checkOverlaps.electrons.requireNoOverlaps = cms.bool(False)

# muon MC match settings
process.load('PhysicsTools.PatAlgos.mcMatchLayer0.muonMatch_cfi')
process.muonMatch.maxDeltaR = cms.double(0.02)
process.muonMatch.resolveByMatchQuality = cms.bool(True)

# make track candidates
makeTrackCandidates(process,    # patAODTrackCands
    label='TrackCands',         # output collection 
    tracks=cms.InputTag('generalTracks'),  # input track collection
    particleType='pi+',         # particle type (for assigning a mass)
    preselection='pt > 0.1',    # preselection cut on candidates
    selection='pt > 0.1',       # Selection on PAT Layer 1 objects 
    isolation={},               # Isolations to use ('source':deltaR; set to {} for None)
    isoDeposits=[],     
    mcAs=None                   # Replicate MC match as the one used for Muons
);

# add the muon associators

process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
addDiMuonTriggers(process)    
useExistingPATMuons(process,'cleanPatMuons' , addL1Info=False)
changeTriggerProcessName(process, 'HLT')
switchOffAmbiguityResolution(process) # Switch off ambiguity resolution: allow multiple reco muons to match to the same trigger muon
useL1MatchingWindowForSinglets(process)
process.muonL1Info.maxDeltaR     =0.3 
process.muonL1Info.fallbackToME1 = True
process.muonMatchHLTL1.maxDeltaR     = 0.3
process.muonMatchHLTL1.fallbackToME1 = True
process.muonMatchHLTL2.maxDeltaR = 0.3
process.muonMatchHLTL2.maxDPtRel = 10.0
process.muonMatchHLTL3.maxDeltaR = 0.1
process.muonMatchHLTL3.maxDPtRel = 10.0
process.muonMatchHLTCtfTrack.maxDeltaR = 0.1
process.muonMatchHLTCtfTrack.maxDPtRel = 10.0
process.muonMatchHLTTrackMu.maxDeltaR = 0.1
process.muonMatchHLTTrackMu.maxDPtRel = 10.0


# define the main process

process.mkcands = cms.EDAnalyzer(
    'B2KstarMuMu',
    treeFileName=cms.string("tree_B2KstarMuMu.root"),
    AnalyzeMC=cms.untracked.bool(True), 

    StoreReco=cms.untracked.bool(True),
    reco_muons=cms.InputTag("patMuons"),

    src=cms.InputTag("patMuonsWithTrigger"),
    #src=cms.InputTag("patMuons"),
    applyCuts= cms.bool(True),
    TrackLabel=cms.InputTag("cleanPatTrackCands"),
    MinTrPt = cms.untracked.double(0.400),
    OniaPiPiMaxDR = cms.untracked.double(1.5),

    # At least one muon must pass this selection
    higherPuritySelection = MuonSelection_v1, 

    # BOTH muons must pass this selection
    lowerPuritySelection = MuonSelection_v1,

    # The dimuon must pass this selection before vertexing
    dimuonSelection  = cms.string("(mass>0.8 && mass<6.0)"), 

    OniaDiTrackRECO=cms.untracked.bool(True),
    # Configuration of trigger matching                           
    triggerResultsLabel = cms.InputTag("TriggerResults","","HLT"),
    HLTBitNames_DoubleMu = HLTBitNames_DoubleMu_v1, 
    HLTLastFilterNames_DoubleMu = HLTLastFilterNames_DoubleMu_v1
    )


# remove unused pat sequences
process.patDefaultSequence.remove(process.patJetCorrFactors)
process.patDefaultSequence.remove(process.patJetCharge)
process.patDefaultSequence.remove(process.patJetPartonMatch)
process.patDefaultSequence.remove(process.patJetGenJetMatch)
process.patDefaultSequence.remove(process.patJetPartons)
process.patDefaultSequence.remove(process.patJetPartonAssociation)
process.patDefaultSequence.remove(process.patJetFlavourAssociation)
process.patDefaultSequence.remove(process.patJets)
process.patDefaultSequence.remove(process.metJESCorAK5CaloJet)
process.patDefaultSequence.remove(process.metJESCorAK5CaloJetMuons)
process.patDefaultSequence.remove(process.patMETs)
process.patDefaultSequence.remove(process.selectedPatJets)
process.patDefaultSequence.remove(process.cleanPatJets)
process.patDefaultSequence.remove(process.countPatJets)


# schedule paths and run 
process.p = cms.Path(process.patDefaultSequence *
                     process.patMuonsWithTriggerSequence * 
                     process.mkcands)

process.schedule = cms.Schedule(process.p)


