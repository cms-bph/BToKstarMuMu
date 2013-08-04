import FWCore.ParameterSet.Config as cms

process = cms.Process("Ntuple")


# Message Service
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger = cms.Service(
    "MessageLogger",
    debugModules = cms.untracked.vstring('*'),
    destinations = cms.untracked.vstring(
        'cerr', 'info'),
        info = cms.untracked.PSet(
            threshold  = cms.untracked.string('INFO')),
        cerr = cms.untracked.PSet(
            threshold  = cms.untracked.string('WARNING'))
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:RECO.root'
    )
)

process.ntuple = cms.EDAnalyzer('BToKMuMu'
)


process.p = cms.Path(process.ntuple)
