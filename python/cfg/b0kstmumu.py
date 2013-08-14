####################
### B0KstMuMu.py ###
####################




#################
### Variables ###
#################
runDataMC = 2 # 1 = Data; 2 = MC (Reco + Gen); 3 = MC (Gen)
useJSON   = False
PrintMsg  = False
triggerProcessName = 'HLT' # 'GEN' or 'HLT' or 'RECO' or 'TEST' or 'REDIGI36X'




#####################
### CMSSW configs ###
#####################
import FWCore.ParameterSet.Config as cms

process = cms.Process('NTUPLIZER')

process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.MessageLogger.suppressWarning = cms.untracked.vstring('B0Cand')


#################
### GlobalTag ###
#################
if (runDataMC != 1):
    process.GlobalTag.globaltag = cms.string('START42_V17::All')
else:
    process.GlobalTag.globaltag = cms.string('GR_R_42_V25::All')




#######################
### Read input file ###
#######################
import sys
if (len(sys.argv) > 2):
    readFiles = sys.argv[2]
# RECO MC:
#    path = 'file:/bestman/storage/cms/store/mc/Summer11/B0ToPsiMuMu_2MuPEtaFilter_7TeV-pythia6-evtgen/GEN-SIM-RECO/PU_S3_START42_V11-v2/'
#    path = 'file:/bestman/storage/cms/store/user/drell/B0ToJPsiKstar_Prod3/B0ToJPsiKstar_Prod1_RECO/623615e3655a310c49f9760910496a19/'
#    path = 'file:/bestman/storage/cms/store/user/drell/B0ToKstarMuMu_Prod4/B0toKstarMuMu_RECO/623615e3655a310c49f9760910496a19/'
#    path = 'file:/nfs/data36/cms/dinardo/B0ToKstMuMu_MyMCRECO001_NTuples/'
    path = 'file:/bestman/storage/cms/store/user/drell/B0dKstMuMuKPi_MuMuEtaFilter_pythia6_evtgen_prod2_1/B0dKstMuMuKPi_MuMuEtaFilter_pythia6_evtgen_RECO_prod2_1/d80aca2d30d1865a7fb9254b7a4518c6/'
# GEN MC:
#    path = 'file:/bestman/storage/cms/store/user/dinardo/PYTHIA6_B0dToKstMuMuKPi_7TeV_cff_GEN_Filter/'
#    path = 'file:/bestman/storage/cms/store/user/dinardo/PYTHIA6_B0dToKstMuMuKPi_7TeV_cff_GEN_NoFilter_01/'
#    path = 'file:/bestman/storage/cms/store/user/dinardo/PYTHIA6_B0dToJPsiKstMuMuKPi_7TeV_cff_GEN_Filter/'
#    path = 'file:/bestman/storage/cms/store/user/dinardo/PYTHIA6_B0dToJPsiKstMuMuKPi_7TeV_cff_GEN_NoFilter/'
#    path = 'file:/bestman/storage/cms/store/user/dinardo/PYTHIA6_B0dToPsi2SKstMuMuKPi_7TeV_cff_GEN_Filter/'
#    path = 'file:/bestman/storage/cms/store/user/dinardo/PYTHIA6_B0dToPsi2SKstMuMuKPi_7TeV_cff_GEN_NoFilter/'
    file = readFiles.replace(path, '')
else:
# RECO MC:
#    from B0ToJPsiMuMu_MC_cff import readFiles
#    from B0ToJPsiKst_MC_cff import readFiles
#    from B0ToKstMuMu_MC_cff import readFiles
#    from B0ToKstMuMu_MyMCRECO001_cff import readFiles
#from B0ToKstMuMu_BrMC1_cff import readFiles
    readFiles =  [
        'root://eoscms//eos/cms/store/user/xshi/B2KstarMuMu_GEN_SIM_HLT_100M_v1/B2KstarMuMu_RECO_100M_v1/74adfdc7814766ad34dc5589b8cb6432/7TeV_PYTHIA6_97_1_sJR.root']
    
# GEN MC:
#    from B0ToKstMuMu_GEN_Filter_MC_cff import readFiles
#    from B0ToKstMuMu_GEN_NoFilter_01_MC_cff import readFiles
#    from B0ToJPsiKst_GEN_Filter_MC_cff import readFiles
#    from B0ToJPsiKst_GEN_NoFilter_MC_cff import readFiles
#    from B0ToPsi2SKst_GEN_Filter_MC_cff import readFiles
#    from B0ToPsi2SKst_GEN_NoFilter_MC_cff import readFiles
process.source = cms.Source('PoolSource',
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            skipEvents = cms.untracked.uint32(0),
                            fileNames = cms.untracked.vstring(
#    'rfio:/castor/cern.ch/user/d/dinardo/MyJPsiKstpMC/MyJPsiKstpSkim_1_1_yR4.root'
    readFiles
    ))


##################################
### Use JSON file interatively ###
##################################
if (runDataMC == 1 and useJSON == True):
    import PhysicsTools.PythonAnalysis.LumiList as LumiList
    import FWCore.ParameterSet.Types as CfgTypes
    myLumis = LumiList.LumiList(filename = 'Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON_MuonPhys.txt').getCMSSWString().split(',')
    process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
    process.source.lumisToProcess.extend(myLumis)


###################
### Output file ###
###################
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))
process.TFileService = cms.Service('TFileService', fileName = cms.string(
    'B0ToKstMuMu.root'
# RECO MC:
#    '/nfs/data36/cms/dinardo/B0ToJPsiMuMu_MC_NTuples/' + readFiles[len(readFiles) - 41 : len(readFiles)]
#    '/nfs/data36/cms/dinardo/B0ToJPsiKst_MC_NTuples/' + file
#    '/nfs/data36/cms/dinardo/B0ToKstMuMu_MC_NTuples/' + file
#    '/nfs/data36/cms/dinardo/B0ToKstMuMu_MyMCANALYSIS001_NTuples/' + file
#    '/nfs/data36/cms/dinardo/ManyTests/' + file
# GEN MC:
#    '/nfs/data36/cms/dinardo/B0ToKstMuMu_GEN_Filter_MC_NTuples/' + file
#    '/nfs/data36/cms/dinardo/B0ToKstMuMu_GEN_NoFilter_01_MC_NTuples/' + file
#    '/nfs/data36/cms/dinardo/B0ToJPsiKst_GEN_Filter_MC_NTuples/' + file
#    '/nfs/data36/cms/dinardo/B0ToJPsiKst_GEN_NoFilter_MC_NTuples/' + file
#    '/nfs/data36/cms/dinardo/B0ToPsi2SKst_GEN_Filter_MC_NTuples/' + file
#    '/nfs/data36/cms/dinardo/B0ToPsi2SKst_GEN_NoFilter_MC_NTuples/' + file
    ), closeFileFast = cms.untracked.bool(True))




###################################
### Import PAT skeleton process ###
###################################
process.load("PhysicsTools.PatAlgos.patSequences_cff")


############################
### Add track candidates ###
############################
process.load("PhysicsTools.PatAlgos.cleaningLayer1.genericTrackCleaner_cfi")
process.cleanPatTracks.checkOverlaps.muons.requireNoOverlaps = cms.bool(False)
process.cleanPatTracks.checkOverlaps.electrons.requireNoOverlaps = cms.bool(False)

process.load('PhysicsTools.PatAlgos.mcMatchLayer0.muonMatch_cfi')
process.muonMatch.maxDeltaR = cms.double(0.02)
process.muonMatch.resolveByMatchQuality = cms.bool(True)

from PhysicsTools.PatAlgos.tools.trackTools import *
if (runDataMC != 1):
    makeTrackCandidates(process,
                        label        = 'TrackCands',                  # output collection
                        tracks       = cms.InputTag('generalTracks'), # input track collection
                        particleType = 'pi+',                         # particle type (for assigning a mass)
                        preselection = 'pt > 0.1',                    # preselection cut on candidates
                        selection    = 'pt > 0.1',                    # selection on PAT Layer 1 objects
                        isolation    = {},                            # isolations to use (set to {} for None)
                        isoDeposits  = [],
                        mcAs         = 'muon'                         # replicate MC match as the one used for Muons
                        )
    
    process.patTrackCandsMCMatch.mcPdgId               = cms.vint32(211) # = pi+
    process.patTrackCandsMCMatch.mcStatus              = cms.vint32(1)
    process.patTrackCandsMCMatch.maxDeltaR             = cms.double(0.02)
    process.patTrackCandsMCMatch.resolveByMatchQuality = cms.bool(True)
else:
    makeTrackCandidates(process,
                        label        = 'TrackCands',                  # output collection
                        tracks       = cms.InputTag('generalTracks'), # input track collection
                        particleType = 'pi+',                         # particle type (for assigning a mass)
                        preselection = 'pt > 0.1',                    # preselection cut on candidates
                        selection    = 'pt > 0.1',                    # selection on PAT Layer 1 objects
                        isolation    = {},                            # isolations to use (set to {} for None)
                        isoDeposits  = [],
                        mcAs         = None                           # replicate MC match as the one used for Muons
                        )
    
    removeMCMatching(process, ['All'], outputInProcess = False)


#####################################
### Do trigger matching for muons ###
#####################################
process.cleanMuonTriggerMatchHLT = cms.EDProducer(
    'PATTriggerMatcherDRLessByR',                          # match by DeltaR only (best match by DeltaR)
    src                   = cms.InputTag('cleanPatMuons'),
    matched               = cms.InputTag('patTrigger'),    # default producer label as defined in PhysicsTools/PatAlgos/python/triggerLayer1/triggerProducer_cfi.py

    # For May10 ReReco only
#    matchedCuts           = cms.string('path("HLT_Dimuon6p5_LowMass_Displaced_v*",0,0)'),

    # After May10 ReReco
    matchedCuts           = cms.string('path("HLT_Dimuon7_LowMass_Displaced_v*") || path("HLT_DoubleMu4_LowMass_Displaced_v*") || path("HLT_DoubleMu4p5_LowMass_Displaced_v*") || path("HLT_DoubleMu5_LowMass_Displaced_v*")'),

    maxDeltaR             = cms.double(0.1),
    resolveAmbiguities    = cms.bool(False),               # only one match per trigger object
    resolveByMatchQuality = cms.bool(False))               # take best match found per reco object (by DeltaR here, see above))


##################################
### Switch on PAT trigger info ###
##################################
from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger(process, hltProcess = triggerProcessName, outputModule = '')
switchOnTriggerMatching(process, triggerMatchers = ['cleanMuonTriggerMatchHLT'], hltProcess = triggerProcessName, outputModule = '')
switchOnTriggerMatchEmbedding(process, triggerMatchers = ['cleanMuonTriggerMatchHLT'], hltProcess = triggerProcessName, outputModule = '')


################################
### Remove not used from PAT ###
################################
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




##########################
### B0 --> K*0 mu+ mu- ###
##########################
process.B0Cand = cms.EDAnalyzer('B0KstMuMu',
                                HLTriggerResults = cms.untracked.string(triggerProcessName),
                                VtxSample        = cms.untracked.string('offlinePrimaryVertices'),
                                BeamSpot         = cms.untracked.string('offlineBeamSpot'),
                                GenParticles     = cms.untracked.string('genParticles'),
                                MuonType         = cms.untracked.string('cleanPatMuonsTriggerMatch'),
                                TrackType        = cms.untracked.string('cleanPatTrackCands'),
                                ParameterFile    = cms.untracked.string('ParameterFile.txt'),
                                doGenReco        = cms.untracked.uint32(runDataMC),
                                printMsg         = cms.untracked.bool(PrintMsg))


###########
### RUN ###
###########
process.patPath  = cms.Path(process.patDefaultSequence)
process.ntupPath = cms.Path(process.B0Cand)

if (runDataMC == 3):
    process.schedule = cms.Schedule(process.ntupPath)
else:
    process.schedule = cms.Schedule(process.patPath, process.ntupPath)
