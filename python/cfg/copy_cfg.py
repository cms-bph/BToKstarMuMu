import FWCore.ParameterSet.Config as cms

# Give the process a name
process = cms.Process("PickEvent")

# Tell the process which files to use as the sourdce
process.source = cms.Source ("PoolSource",
#          fileNames = cms.untracked.vstring ("/store/mc/Summer12_DR53X/BuToKstarJPsi_EtaPtFilter_8TeV-pythia6-evtgen/AODSIM/PU_RD2_START53_V19F-v1/00000/0010B7CC-994C-E311-BB5E-002590D0AF74.root")
#          fileNames = cms.untracked.vstring ("/store/mc/Summer12_DR53X/BuToKstarMuMu_EtaPtFilter_8TeV-pythia6-evtgen/AODSIM/PU_RD2_START53_V19F-v1/00000/006C6D81-8B4D-E311-BEBC-E0CB4E1A1194.root")
          fileNames = cms.untracked.vstring ("/store/mc/Summer12_DR53X/BuToKstarPsi2S_EtaPtFilter_8TeV-pythia6-evtgen/AODSIM/PU_RD2_START53_V19F-v1/20000/001FC228-B948-E311-93B2-003048FFD7C2.root")                   
)

# tell the process to only run over 500 events (-1 would mean run over
#  everything
process.maxEvents = cms.untracked.PSet(
            input = cms.untracked.int32 (500)

)

# Tell the process what filename to use to save the output
process.Out = cms.OutputModule("PoolOutputModule",
#         fileName = cms.untracked.string ("/afs/cern.ch/work/n/nsahoo/BPH-ANALYSIS/afb/SE/aod/mc/2012_test_BuToKstarJPsi_500.root")
#         fileName = cms.untracked.string ("/afs/cern.ch/work/n/nsahoo/BPH-ANALYSIS/afb/SE/aod/mc/2012_test_BuToKstarMuMu_500.root")
         fileName = cms.untracked.string ("/afs/cern.ch/work/n/nsahoo/BPH-ANALYSIS/afb/SE/aod/mc/2012_test_BuToKstarPsi2S_500.root")

)

# make sure everything is hooked up
process.end = cms.EndPath(process.Out)
