import FWCore.ParameterSet.Config as cms

process = cms.Process('NTUPLE')

process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(False)
    )
# import of standard configurations
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.suppressInfo = cms.untracked.vstring( "mkcands" )
process.MessageLogger.suppressWarning = cms.untracked.vstring( "mkcands" )
process.MessageLogger.cerr.FwkReport.reportEvery = 10
MC=True
# Input source
process.source = cms.Source("PoolSource",
                            #skipEvents = cms.untracked.uint32(0),
                            fileNames = cms.untracked.vstring(
                                'file:RECO.root'
                                )
    )

process.source.inputCommands = cms.untracked.vstring(
    "keep *",
    "drop L1GlobalTriggerObjectMapRecord_hltL1GtObjectMap__RECO",
    "drop *_MEtoEDMConverter_*_*"
    )

process.maxEvents = cms.untracked.PSet(  input = cms.untracked.int32(10) )


process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")


process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
## GlobalTag prompt processing of 2011 data GR_P_V16
#process.GlobalTag.globaltag = 'GR_P_V22A::All' 
#process.GlobalTag.globaltag = 'GR_P_V22::All' 
##process.GlobalTag.globaltag = 'GR_R_311_V4::All'
process.GlobalTag.globaltag = 'START42_V14A::All' 


process.load('Configuration/EventContent/EventContent_cff')
#
#  Load common sequences
#
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskAlgoTrigConfig_cff')
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')



####################################################################################
##################################good collisions############################################

############
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
	 																				 vertexCollection = cms.InputTag('offlinePrimaryVertices'),
																					 minimumNDOF = cms.uint32(4) ,
																					 maxAbsZ = cms.double(24),	
																					 maxd0 = cms.double(2)	
)

process.noscraping = cms.EDFilter("FilterOutScraping",
                                  applyfilter = cms.untracked.bool(True),
																 	debugOn = cms.untracked.bool(False),
																 	numtrack = cms.untracked.uint32(10),
																 	thresh = cms.untracked.double(0.25)
)

# PAT Layer 0+1
process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.load("PhysicsTools.PatAlgos.cleaningLayer1.genericTrackCleaner_cfi")
process.cleanPatTracks.checkOverlaps.muons.requireNoOverlaps = cms.bool(False)
process.cleanPatTracks.checkOverlaps.electrons.requireNoOverlaps = cms.bool(False)
from PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi import *
patMuons.embedTrack = cms.bool(True)
patMuons.embedPickyMuon = cms.bool(False)
patMuons.embedTpfmsMuon = cms.bool(False)

# Prune generated particles to muons and their parents
process.genMuons = cms.EDProducer("GenParticlePruner",
        src = cms.InputTag("genParticles"),
        select = cms.vstring(
            "drop  *  ",                     # this is the default
            "++keep abs(pdgId) = 13",        # keep muons and their parents
            "++keep abs(pdgId) = 321",        # keep Kaons and their parents
            "++keep abs(pdgId) = 211",        # keep Pions and their parents
            "drop pdgId == 21 && status = 2" # remove intermediate qcd spam carrying no flavour info
      )
 )



process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import  addMCinfo, useExistingPATMuons, useL1MatchingWindowForSinglets, changeTriggerProcessName, switchOffAmbiguityResolution, addDiMuonTriggers
    # with some customization
if MC:
        addMCinfo(process)
        # since we match inner tracks, keep the matching tight and make it one-to-one
        process.muonMatch.maxDeltaR = 0.02
        process.muonMatch.resolveByMatchQuality = True
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


    
from PhysicsTools.PatAlgos.tools.trackTools import *
makeTrackCandidates(process,                                        #         patAODTrackCands
        label='TrackCands',                   # output collection will be 'allLayer0TrackCands', 'allLayer1TrackCands', 'selectedLayer1TrackCands'
        tracks=cms.InputTag('generalTracks'), # input track collection
        particleType='pi+',                   # particle type (for assigning a mass)
        preselection='pt > 0.1',              # preselection cut on candidates. Only methods of 'reco::Candidate' are available
        selection='pt > 0.1',                 # Selection on PAT Layer 1 objects ('selectedLayer1TrackCands')
        isolation={},                         # Isolations to use ('source':deltaR; set to {} for None)
        isoDeposits=[],
        mcAs=None                           # Replicate MC match as the one used for Muons
        );                                    # you can specify more than one collection for this

l1cands = getattr(process, 'patTrackCands')
l1cands.addGenMatch = False

process.mkcands = cms.EDAnalyzer('B2KstarMuMu',
      treeFileName=cms.string("tree_B2KstarMuMu.root"),
      AnalyzeMC=cms.untracked.bool(True),


      # -----------------------------
      # For StoreReco 
      # -----------------------------
      StoreReco=cms.untracked.bool(True),
      reco_muons=cms.InputTag("patMuons"),
      # -----------------------------
      
      src=cms.InputTag("patMuonsWithTrigger"),
      applyCuts= cms.bool(True),
      TrackLabel=cms.InputTag("cleanPatTrackCands"),
      MinTrPt = cms.untracked.double(0.400),
      OniaPiPiMaxDR = cms.untracked.double(1.5),
      # At least one muon must pass this selection
      higherPuritySelection = cms.string("(isGlobalMuon || isTrackerMuon || (innerTrack.isNonnull && genParticleRef(0).isNonnull)) && abs(innerTrack.dxy)<4 && abs(innerTrack.dz)<35 && muonID('TrackerMuonArbitrated')"),
      # BOTH muons must pass this selection
      lowerPuritySelection  = cms.string("(isGlobalMuon || isTrackerMuon || (innerTrack.isNonnull && genParticleRef(0).isNonnull)) && abs(innerTrack.dxy)<4 && abs(innerTrack.dz)<35 && muonID('TrackerMuonArbitrated')"),
      
      dimuonSelection  = cms.string("(mass>0.8 && mass<6.0)"), ## The dimuon must pass this selection before vertexing
		 	OniaDiTrackRECO=cms.untracked.bool(True), # Onia plus two track selection
	  # Configuration of trigger matching                           
      triggerResultsLabel = cms.InputTag("TriggerResults","","HLT"),

      ####### 2011, 5E32 1E33,1.4E33 ########
      HLTBitNames_DoubleMu = cms.vstring(
         "HLT_DoubleMu3_Jpsi_v1", # 5E32, v4.2-v5.3
         "HLT_DoubleMu3_Jpsi_v2", # 5E32, v6.1-v6.2
         "HLT_Dimuon6p5_Jpsi_v1", #5E32, v8.1-v8.3
         "HLT_Dimuon6p5_Jpsi_Displaced_v1", #5E32, v8.1-v8.3
         "HLT_Dimuon6p5_Barrel_Jpsi_v1", #5E32, v8.1-v8.3
         "HLT_DoubleMu3_Quarkonium_v1", #5E32, v4.2-v5.3
         "HLT_DoubleMu3_Quarkonium_v2", #5E32, v6.1-v6.2
         "HLT_DoubleMu3_LowMass_v1", #5E32, v6.1-v6.2
         "HLT_Dimuon6p5_LowMass_v1", #5E32, v8.1-v8.3
         "HLT_Dimuon6p5_LowMass_Displaced_v1", #5E32, v8.1-v8.3
         "HLT_DoubleMu3_Bs_v1", #5E32, v4.2-v5.3
         "HLT_DoubleMu2_Bs_v1", #5E32, v6.1-v6.2
         "HLT_DoubleMu2_Bs_v2", #5E32, v8.1-v8.3
         "HLT_DoubleMu2_Bs_v3", #1E33, v1.2-v2.4
         "HLT_DoubleMu2_Bs_v4", #1E33, v2.5
         "HLT_DoubleMu2_Bs_v5", #1.4E33, v1.2
         "HLT_Dimuon0_Jpsi_v1", #1E33, v1.2-v2.4
         "HLT_Dimuon0_Jpsi_v2", #1E33, v2.5
         "HLT_Dimuon0_Jpsi_v3", #1.4E33, v1.2
         "HLT_Dimuon0_Jpsi_v5", #2E33, v1.1-v1.2
         "HLT_Dimuon0_Jpsi_v6", #3E33,  v1.1-v5.0
         "HLT_Dimuon0_Jpsi_v9", #5E33,  v1.4
         "HLT_Dimuon0_Jpsi_NoVertexing_v2", #2E33, v1.1-v1.2
         "HLT_Dimuon0_Jpsi_NoVertexing_v3", #3E33, v1.1-v5.0
         "HLT_Dimuon0_Jpsi_NoVertexing_v6", #5E33, v1.4
         "HLT_Dimuon0_Upsilon_v1", #1E33, v1.2-v2.4
         "HLT_Dimuon0_Upsilon_v2",#1E33, v2.5
         "HLT_Dimuon0_Upsilon_v3",#1.4E33, v1.2
         "HLT_Dimuon0_Upsilon_v5",#2E33, v1.1-v1.2
         "HLT_Dimuon0_Upsilon_v6",#3E33, v1.1-v5.0
         "HLT_Dimuon0_Upsilon_v9",#5E33, v1.4
         "HLT_Dimuon4_Bs_Barrel_v2",#1E33, v1.2-v1.3
         "HLT_Dimuon4_Bs_Barrel_v3", #1E33, v2.2-v2.4
         "HLT_Dimuon4_Bs_Barrel_v4", #1E33, v2.5
         "HLT_Dimuon4_Bs_Barrel_v5",#1.4E33, v1.2
         "HLT_Dimuon4_Bs_Barrel_v7",#2E33, v1.1-v1.2
         "HLT_DoubleMu4_Dimuon4_Bs_Barrel_v1", #3E33, v1.1-v5.0
         "HLT_DoubleMu4_Dimuon4_Bs_Barrel_v4",#5E33, v1.4
         "HLT_Dimuon5_Upsilon_Barrel_v1", #1E33, v1.2-v2.4
         "HLT_Dimuon5_Upsilon_Barrel_v2",#1E33, v2.5
         "HLT_Dimuon5_Upsilon_Barrel_v3",#1.4E33, v1.2
         "HLT_Dimuon5_Upsilon_Barrel_v5",#2E33, v1.1-v1.2
         "HLT_Dimuon7_Upsilon_Barrel_v1",#3E33, v1.1-v5.0
         "HLT_Dimuon7_Upsilon_Barrel_v4",#5E33, v1.4
         "HLT_Dimuon9_Upsilon_Barrel_v1",#3E33, v1.1-v5.0
         "HLT_Dimuon9_Upsilon_Barrel_v4",#5E33, v1.4
         "HLT_Dimuon6_Bs_v1",#1E33, v1.2-v1.3
         "HLT_Dimuon6_Bs_v2",#1E33, v2.2-v2.4
         "HLT_Dimuon6_Bs_v3",#1E33, v2.5
         "HLT_Dimuon6_Bs_v4",#1.4E33, v1.2
         "HLT_Dimuon6_Bs_v6",#2E33, v1.1-v1.2
         "HLT_DoubleMu4_Dimuon6_Bs_v1",#3E33, v1.1-v5.0         
         "HLT_DoubleMu4_Dimuon6_Bs_v4",#5E33, v1.4         
         "HLT_Dimuon7_LowMass_Displaced_v1",#1E33, v1.2-v1.3
         "HLT_Dimuon7_LowMass_Displaced_v2", #1E33, v2.2-v2.4
         "HLT_Dimuon7_LowMass_Displaced_v3",#1E33, v2.5
         "HLT_Dimuon7_LowMass_Displaced_v4",#1.4E33, v1.2
         "HLT_DoubleMu4_LowMass_Displaced_v2",#2E33, v1.1-v1.2
         "HLT_DoubleMu4p5_LowMass_Displaced_v1",#3E33, v1.1-v5.0
         "HLT_DoubleMu4p5_LowMass_Displaced_v4",#5E33, v1.4
         "HLT_DoubleMu5_LowMass_Displaced_v1",#3E33, v1.1-v5.0
         "HLT_DoubleMu5_LowMass_Displaced_v4",#5E33, v1.4
         "HLT_Dimuon7_Jpsi_Displaced_v1", #1E33, v1.2-v2.4
         "HLT_Dimuon7_Jpsi_Displaced_v2",#1E33, v2.5
         "HLT_Dimuon7_Jpsi_Displaced_v3",#1.4E33, v1.2
         "HLT_DoubleMu3p5_Jpsi_Displaced_v2", #2E33, v1.1-v1.2
         "HLT_DoubleMu4_Jpsi_Displaced_v1", #3E33, v1.1-v5.0
         "HLT_DoubleMu4_Jpsi_Displaced_v4", #5E33, v1.4
         "HLT_DoubleMu5_Jpsi_Displaced_v1", #3E33, v1.1-v3.1
         "HLT_DoubleMu5_Jpsi_Displaced_v2", #3E33, v4.0-v5.0
         "HLT_DoubleMu5_Jpsi_Displaced_v4", #5E33, v1.4
         "HLT_Dimuon7_Jpsi_X_Barrel_v1", #1E33, v1.2-v2.4
         "HLT_Dimuon7_Jpsi_X_Barrel_v2",#1E33, v2.5
         "HLT_Dimuon7_Jpsi_X_Barrel_v3",#1.4E33, v1.2
         "HLT_Dimuon7_Jpsi_X_Barrel_v5",#2E33, v1.1-v1.2
         "HLT_Dimuon7_PsiPrime_v1", #1E33, v1.2-v2.4
         "HLT_Dimuon7_PsiPrime_v2",#1E33, v2.5
         "HLT_Dimuon7_PsiPrime_v3",#1.4E33, v1.2
         "HLT_Dimuon7_PsiPrime_v5",#2E33, v1.1-v1.2
         "HLT_Dimuon9_PsiPrime_v1",#3E33, v1.1-v1.2
         "HLT_Dimuon9_PsiPrime_v4",#5E33, v1.4
         "HLT_Dimuon11_PsiPrime_v1",#3E33, v1.1-v1.2
         "HLT_Dimuon11_PsiPrime_v4",#5E33, v1.4
         "HLT_Dimuon10_Jpsi_Barrel_v1", #1E33, v1.2-v2.4
         "HLT_Dimuon10_Jpsi_Barrel_v2",#1E33, v2.5
         "HLT_Dimuon10_Jpsi_Barrel_v3",#1.4E33, v1.2
         "HLT_Dimuon10_Jpsi_Barrel_v5",#2E33, v1.1-v1.2
         "HLT_Dimuon10_Jpsi_Barrel_v6",#3E33, v1.1-v5.0
         "HLT_Dimuon10_Jpsi_Barrel_v9",#5E33, v1.4
         "HLT_Dimuon13_Jpsi_Barrel_v1",#3E33, v1.1-v5.0
         "HLT_Dimuon13_Jpsi_Barrel_v4",#5E33, v1.4
         "HLT_Dimuon0_Jpsi_Muon_v1", #1E33, v1.2-v1.3
         "HLT_Dimuon0_Jpsi_Muon_v2", #1E33, v2.2-v2.4
         "HLT_Dimuon0_Jpsi_Muon_v3",#1E33, v2.5
         "HLT_Dimuon0_Jpsi_Muon_v4",#1.4E33, v1.2
         "HLT_Dimuon0_Jpsi_Muon_v6", #2E33, v1.1-v1.2
         "HLT_Dimuon0_Jpsi_Muon_v7", #3E33, v1.1-v5.0
         "HLT_Dimuon0_Jpsi_Muon_v10", #5E33, v1.4
         "HLT_Dimuon0_Upsilon_Muon_v1",#1E33, v1.2-v1.3
         "HLT_Dimuon0_Upsilon_Muon_v2", #1E33, v2.2-v2.4
         "HLT_Dimuon0_Upsilon_Muon_v3",#1E33, v2.5
         "HLT_Dimuon0_Upsilon_Muon_v4",#1.4E33, v1.2
         "HLT_Dimuon0_Upsilon_Muon_v6",#2E33, v1.1-v1.2
         "HLT_Dimuon0_Upsilon_Muon_v7",#3E33, v1.1-v5.0
         "HLT_Dimuon0_Upsilon_Muon_v10",#5E33, v1.4
        ),
   # ONE FILTER NAME PER PATH    
   HLTLastFilterNames_DoubleMu = cms.vstring(
         "hltDoubleMu3JpsiL3Filtered", #HLT_DoubleMu3_Jpsi_v1
         "hltDoubleMu3JpsiL3Filtered", #HLT_DoubleMu3_Jpsi_v2
         "hltDimuon6p5JpsiL3Filtered", #HLT_Dimuon6p5_Jpsi_v1
         "hltDimuon6p5JpsiDisplacedL3Filtered", #HLT_Dimuon6p5_Jpsi_Displaced_v1
         "hltDimuon6p5BarrelJpsiL3Filtered", #HLT_Dimuon6p5_Barrel_Jpsi_v1
         "hltDoubleMu3QuarkoniumL3Filtered", #HLT_DoubleMu3_Quarkonium_v1	      
         "hltDoubleMu3QuarkoniumL3Filtered", #HLT_DoubleMu3_Quarkonium_v2
         "hltDoubleMu3LowMassL3Filtered", #HLT_DoubleMu3_LowMass_v1
         "hltDimuon6p5LowMassL3Filtered", #HLT_Dimuon6p5_LowMass_v1
         "hltDimuon6p5LowMassL3FilteredDisplaced", #HLT_Dimuon6p5_LowMass_Displaced_v1
         "hltDoubleMu3BsL3Filtered", #HLT_DoubleMu3_Bs_v1
         "hltDoubleMu2BsL3Filtered", #HLT_DoubleMu2_Bs_v1
         "hltDoubleMu2BsL3Filtered", #HLT_DoubleMu2_Bs_v2
         "hltDoubleMu2BsL3Filtered", #HLT_DoubleMu2_Bs_v3
         "hltDoubleMu2BsL3Filtered", #HLT_DoubleMu2_Bs_v4
         "hltDoubleMu2BsL3Filtered", #HLT_DoubleMu2_Bs_v5
         "hltVertexmumuFilterJpsi", #HLT_Dimuon0_Jpsi_v1 was "hltJpsiL3Filtered"
         "hltVertexmumuFilterJpsi", #HLT_Dimuon0_Jpsi_v2 was "hltJpsiL3Filtered"
         "hltVertexmumuFilterJpsi", #HLT_Dimuon0_Jpsi_v3 was "hltJpsiL3Filtered"
         "hltVertexmumuFilterJpsi", #HLT_Dimuon0_Jpsi_v5
         "hltVertexmumuFilterJpsi", #HLT_Dimuon0_Jpsi_v6
         "hltVertexmumuFilterJpsi", #HLT_Dimuon0_Jpsi_v9
         "hltJpsiNoVertexingL3Filtered", #HLT_Dimuon0_Jpsi_NoVertexing_v2
         "hltJpsiNoVertexingL3Filtered", #HLT_Dimuon0_Jpsi_NoVertexing_v3
         "hltJpsiNoVertexingL3Filtered", #HLT_Dimuon0_Jpsi_NoVertexing_v6
         "hltVertexmumuFilterUpsilon", #HLT_Dimuon0_Upsilon_v1 was "hltUpsilonL3Filtered"
         "hltVertexmumuFilterUpsilon", #HLT_Dimuon0_Upsilon_v2 was "hltUpsilonL3Filtered"
         "hltVertexmumuFilterUpsilon", #HLT_Dimuon0_Upsilon_v3 was "hltUpsilonL3Filtered"
         "hltVertexmumuFilterUpsilon", #HLT_Dimuon0_Upsilon_v5
         "hltVertexmumuFilterUpsilon", #HLT_Dimuon0_Upsilon_v6
         "hltVertexmumuFilterUpsilon", #HLT_Dimuon0_Upsilon_v9
         "hltDoubleMu2BarrelBsL3Filtered", #HLT_Dimuon4_Bs_Barrel_v2
         "hltDoubleMu2BarrelBsL3Filtered", #HLT_Dimuon4_Bs_Barrel_v3
         "hltDoubleMu2BarrelBsL3Filtered", #HLT_Dimuon4_Bs_Barrel_v4
         "hltDoubleMu2BarrelBsL3Filtered", #HLT_Dimuon4_Bs_Barrel_v5
         "hltDoubleMu2BarrelBsL3Filtered", #HLT_Dimuon4_Bs_Barrel_v7
         "hltVertexmumuFilterBs4", #HLT_DoubleMu4_Dimuon4_Bs_Barrel_v1
         "hltVertexmumuFilterBs4", #HLT_DoubleMu4_Dimuon4_Bs_Barrel_v4
         "hltVertexmumuFilterUpsilonBarrel", #HLT_Dimuon5_Upsilon_Barrel_v1 was "hltBarrelUpsilonL3Filtered"
         "hltVertexmumuFilterUpsilonBarrel", #HLT_Dimuon5_Upsilon_Barrel_v2 was "hltBarrelUpsilonL3Filtered"
         "hltVertexmumuFilterUpsilonBarrel", #HLT_Dimuon5_Upsilon_Barrel_v3 was "hltBarrelUpsilonL3Filtered"
         "hltVertexmumuFilterUpsilonBarrel", #HLT_Dimuon5_Upsilon_Barrel_v5
         "hltVertexmumuFilterDimuon7UpsilonBarrel", #HLT_Dimuon7_Upsilon_Barrel_v1
         "hltVertexmumuFilterDimuon7UpsilonBarrel", #HLT_Dimuon7_Upsilon_Barrel_v4
         "hltVertexmumuFilterDimuon9UpsilonBarrel", #HLT_Dimuon9_Upsilon_Barrel_v1
         "hltVertexmumuFilterDimuon9UpsilonBarrel", #HLT_Dimuon9_Upsilon_Barrel_v4
         "hltDoubleMu2Dimuon6BsL3Filtered", #HLT_Dimuon6_Bs_v1
         "hltDoubleMu2Dimuon6BsL3Filtered", #HLT_Dimuon6_Bs_v2
         "hltDoubleMu2Dimuon6BsL3Filtered", #HLT_Dimuon6_Bs_v3
         "hltDoubleMu2Dimuon6BsL3Filtered", #HLT_Dimuon6_Bs_v4
         "hltDoubleMu2Dimuon6BsL3Filtered", #HLT_Dimuon6_Bs_v6
         "hltVertexmumuFilterBs6", #HLT_DoubleMu4_Dimuon6_Bs_v1 
         "hltVertexmumuFilterBs6", #HLT_DoubleMu4_Dimuon6_Bs_v4 
         "hltDisplacedmumuFilterLowMass", #HLT_Dimuon7_LowMass_Displaced_v1 was "hltLowMassDisplacedL3Filtered"
         "hltDisplacedmumuFilterLowMass", #HLT_Dimuon7_LowMass_Displaced_v2 was "hltLowMassDisplacedL3Filtered"
         "hltDisplacedmumuFilterLowMass", #HLT_Dimuon7_LowMass_Displaced_v3 was "hltLowMassDisplacedL3Filtered"
         "hltDisplacedmumuFilterLowMass", #HLT_Dimuon7_LowMass_Displaced_v4 was "hltLowMassDisplacedL3Filtered"
         "hltDisplacedmumuFilterLowMass", #HLT_DoubleMu4_LowMass_Displaced_v2
         "hltDisplacedmumuFilterDoubleMu4p5LowMass",#HLT_DoubleMu4p5_LowMass_Displaced_v1
         "hltDisplacedmumuFilterDoubleMu4p5LowMass",#HLT_DoubleMu4p5_LowMass_Displaced_v4
         "hltDisplacedmumuFilterDoubleMu5LowMass",#HLT_DoubleMu5_LowMass_Displaced_v1
         "hltDisplacedmumuFilterDoubleMu5LowMass",#HLT_DoubleMu5_LowMass_Displaced_v4
         "hltDisplacedmumuFilterJpsi", #HLT_Dimuon7_Jpsi_Displaced_v1 was "hltJpsiDisplacedL3Filtered"
         "hltDisplacedmumuFilterJpsi", #HLT_Dimuon7_Jpsi_Displaced_v2 was "hltJpsiDisplacedL3Filtered"
         "hltDisplacedmumuFilterJpsi", #HLT_Dimuon7_Jpsi_Displaced_v3 was "hltJpsiDisplacedL3Filtered"
         "hltDisplacedmumuFilterJpsi", #HLT_DoubleMu3p5_Jpsi_Displaced_v2
         "hltDisplacedmumuFilterDoubleMu4Jpsi", #HLT_DoubleMu4_Jpsi_Displaced_v1
         "hltDisplacedmumuFilterDoubleMu4Jpsi", #HLT_DoubleMu4_Jpsi_Displaced_v4
         "hltDisplacedmumuFilterDoubleMu5Jpsi", #HLT_DoubleMu5_Jpsi_Displaced_v1
         "hltDisplacedmumuFilterDoubleMu5Jpsi", #HLT_DoubleMu5_Jpsi_Displaced_v2
         "hltDisplacedmumuFilterDoubleMu5Jpsi", #HLT_DoubleMu5_Jpsi_Displaced_v4          
         "hltVertexmumuFilterJpsiXBarrel", #HLT_Dimuon7_Jpsi_X_Barrel_v1 was "hltJpsiXBarrelL3Filtered"
         "hltVertexmumuFilterJpsiXBarrel", #HLT_Dimuon7_Jpsi_X_Barrel_v2 was "hltJpsiXBarrelL3Filtered"
         "hltVertexmumuFilterJpsiXBarrel", #HLT_Dimuon7_Jpsi_X_Barrel_v3 was "hltJpsiXBarrelL3Filtered"
         "hltVertexmumuFilterJpsiXBarrel", #HLT_Dimuon7_Jpsi_X_Barrel_v5
         "hltVertexmumuFilterPsiPrime", #HLT_Dimuon7_PsiPrime_v1 was "hltPsiPrimeL3Filtered"
         "hltVertexmumuFilterPsiPrime", #HLT_Dimuon7_PsiPrime_v2 was "hltPsiPrimeL3Filtered"
         "hltVertexmumuFilterPsiPrime", #HLT_Dimuon7_PsiPrime_v3 was "hltPsiPrimeL3Filtered"
         "hltVertexmumuFilterPsiPrime", #HLT_Dimuon7_PsiPrime_v5
         "hltVertexmumuFilterDimuon9PsiPrime", #HLT_Dimuon9_PsiPrime_v1
         "hltVertexmumuFilterDimuon9PsiPrime", #HLT_Dimuon9_PsiPrime_v4
         "hltVertexmumuFilterDimuon11PsiPrime", #HLT_Dimuon11_PsiPrime_v1
         "hltVertexmumuFilterDimuon11PsiPrime", #HLT_Dimuon11_PsiPrime_v4
         "hltVertexmumuFilterJpsiBarrel", #HLT_Dimuon10_Jpsi_Barrel_v1 was "hltBarrelJpsiL3Filtered"
         "hltVertexmumuFilterJpsiBarrel", #HLT_Dimuon10_Jpsi_Barrel_v2 was "hltBarrelJpsiL3Filtered"
         "hltVertexmumuFilterJpsiBarrel", #HLT_Dimuon10_Jpsi_Barrel_v3 was "hltBarrelJpsiL3Filtered"
         "hltVertexmumuFilterJpsiBarrel", #HLT_Dimuon10_Jpsi_Barrel_v5
         "hltVertexmumuFilterDimuon10JpsiBarrel", #HLT_Dimuon10_Jpsi_Barrel_v6
         "hltVertexmumuFilterDimuon10JpsiBarrel", #HLT_Dimuon10_Jpsi_Barrel_v9
         "hltVertexmumuFilterDimuon13JpsiBarrel", #HLT_Dimuon13_Jpsi_Barrel_v1
         "hltVertexmumuFilterDimuon13JpsiBarrel", #HLT_Dimuon13_Jpsi_Barrel_v4
         "hltVertexmumuFilterJpsiMuon", #HLT_Dimuon0_Jpsi_Muon_v1 was "hltJpsiMuonL3Filtered"
         "hltVertexmumuFilterJpsiMuon", #HLT_Dimuon0_Jpsi_Muon_v2 was "hltJpsiMuonL3Filtered"
         "hltVertexmumuFilterJpsiMuon", #HLT_Dimuon0_Jpsi_Muon_v3 was "hltJpsiMuonL3Filtered"
         "hltVertexmumuFilterJpsiMuon", #HLT_Dimuon0_Jpsi_Muon_v4 was "hltJpsiMuonL3Filtered"
         "hltVertexmumuFilterJpsiMuon", #HLT_Dimuon0_Jpsi_Muon_v6
         "hltVertexmumuFilterJpsiMuon", #HLT_Dimuon0_Jpsi_Muon_v7
         "hltVertexmumuFilterJpsiMuon", #HLT_Dimuon0_Jpsi_Muon_v10
         "hltVertexmumuFilterUpsilonMuon", #HLT_Dimuon0_Upsilon_Muon_v1 was "hltUpsilonMuonL3Filtered"
         "hltVertexmumuFilterUpsilonMuon", #HLT_Dimuon0_Upsilon_Muon_v2 was "hltUpsilonMuonL3Filtered"
         "hltVertexmumuFilterUpsilonMuon", #HLT_Dimuon0_Upsilon_Muon_v3 was "hltUpsilonMuonL3Filtered"
         "hltVertexmumuFilterUpsilonMuon", #HLT_Dimuon0_Upsilon_Muon_v4 was "hltUpsilonMuonL3Filtered"
         "hltVertexmumuFilterUpsilonMuon", #HLT_Dimuon0_Upsilon_Muon_v6
         "hltVertexmumuFilterUpsilonMuon", #HLT_Dimuon0_Upsilon_Muon_v7
         "hltVertexmumuFilterUpsilonMuon", #HLT_Dimuon0_Upsilon_Muon_v10
),

)

# turn off MC matching for the process
from PhysicsTools.PatAlgos.tools.coreTools import *
removeMCMatching(process, ['All'], outputInProcess = False)

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



process.filter = cms.Sequence(process.noscraping)
#process.filter = cms.Sequence(process.primaryVertexFilter+process.noscraping)

## AF: replace EndPath with path
process.ntup = cms.Path(process.filter*process.patDefaultSequence*process.patMuonsWithTriggerSequence*process.mkcands)


process.schedule = cms.Schedule(process.ntup)
