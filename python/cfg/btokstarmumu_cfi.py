import FWCore.ParameterSet.Config as cms

process = cms.Process("Ntuple")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.MessageLogger = cms.Service(
    "MessageLogger",
    destinations = cms.untracked.vstring('cerr', 'cout', 'message'),
    categories = cms.untracked.vstring('myHLT', 'myVertex', 'myDimuon', 'myKshort',
                                       'myMuonMatch', 'myKstarCharged',
                                       'myKshortMatch', 'myBu'),
    cerr = cms.untracked.PSet(threshold = cms.untracked.string('WARNING')),
    cout = cms.untracked.PSet(
        threshold = cms.untracked.string('INFO'),
        INFO = cms.untracked.PSet(limit = cms.untracked.int32(0)), 
        myHLT = cms.untracked.PSet(limit = cms.untracked.int32(0)), 
        myVertex = cms.untracked.PSet(limit = cms.untracked.int32(0)), 
        myDimuon = cms.untracked.PSet(limit = cms.untracked.int32(-1)), 
        myKshort = cms.untracked.PSet(limit = cms.untracked.int32(-1)), 
        myMuonMatch = cms.untracked.PSet(limit = cms.untracked.int32(0)), 
        myKstarCharged = cms.untracked.PSet(limit = cms.untracked.int32(-1)), 
        myKshortMatch = cms.untracked.PSet(limit = cms.untracked.int32(0)), 
        myBu = cms.untracked.PSet(limit = cms.untracked.int32(-1)), 
        ), 
     message = cms.untracked.PSet(
        threshold = cms.untracked.string('INFO'),
        INFO = cms.untracked.PSet(limit = cms.untracked.int32(0)), 
        myHLT = cms.untracked.PSet(limit = cms.untracked.int32(0)), 
        myVertex = cms.untracked.PSet(limit = cms.untracked.int32(0)), 
        myDimuon = cms.untracked.PSet(limit = cms.untracked.int32(-1)), 
        myKshort = cms.untracked.PSet(limit = cms.untracked.int32(-1)), 
        myMuonMatch = cms.untracked.PSet(limit = cms.untracked.int32(0)), 
        myKstarCharged = cms.untracked.PSet(limit = cms.untracked.int32(-1)), 
        myKshortMatch = cms.untracked.PSet(limit = cms.untracked.int32(0)), 
        myBu = cms.untracked.PSet(limit = cms.untracked.int32(-1)), 
        )
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load("PhysicsTools.PatAlgos.patSequences_cff")

# add track candidates
from PhysicsTools.PatAlgos.tools.trackTools import *
makeTrackCandidates(process,
                    label        = 'TrackCands',                  
                    tracks       = cms.InputTag('generalTracks'), 
                    particleType = 'pi+',                         
                    preselection = 'pt > 0.1',                     
                    selection    = 'pt > 0.1',                     
                    isolation    = {},                            
                    isoDeposits  = [],                            
                    mcAs         = None          
    )    

from PhysicsTools.PatAlgos.tools.coreTools import *
removeMCMatching(process, ['All'], outputInProcess = False)


process.localV0Candidates = cms.EDProducer(
    "V0Producer",
                                     
    # InputTag that tells which TrackCollection to use for vertexing
    trackRecoAlgorithm = cms.InputTag('generalTracks'),

    # These bools decide whether or not to reconstruct
    #  specific V0 particles
    selectKshorts = cms.bool(True),
    selectLambdas = cms.bool(True),

    # Recommend leaving this one as is.
    vertexFitter = cms.InputTag('KalmanVertexFitter'),

    # set to true, uses tracks refit by the KVF for V0Candidate kinematics
    #  NOTE: useSmoothing is automatically set to FALSE
    #  if using the AdaptiveVertexFitter (which is NOT recommended)
    useSmoothing = cms.bool(True),
                                     
    # Select tracks using TrackBase::TrackQuality.
    # Select ALL tracks by leaving this vstring empty, which
    #   is equivalent to using 'loose'
    #trackQualities = cms.vstring('highPurity', 'goodIterative'),
    trackQualities = cms.vstring('loose'),
                                     
    # The next parameters are cut values.
    # All distances are in cm, all energies in GeV, as usual.

    # --Track quality/compatibility cuts--
    #   Normalized track Chi2 <
    tkChi2Cut = cms.double(5.0),
    #   Number of valid hits on track >=
    tkNhitsCut = cms.int32(6),
    #   Track impact parameter significance >
    #impactParameterSigCut = cms.double(2.),
    impactParameterSigCut = cms.double(0.5),
    # We calculate the PCA of the tracks quickly in RPhi, extrapolating
    # the z position as well, before vertexing.  Used in the following 2 cuts:
    #   m_pipi calculated at PCA of tracks <
    mPiPiCut = cms.double(0.6),
    #   PCA distance between tracks <
    tkDCACut = cms.double(1.),

    # --V0 Vertex cuts--
    #   Vertex chi2 < 
    vtxChi2Cut = cms.double(7.0),
    #   Lambda collinearity cut
    #   (UNUSED)
    collinearityCut = cms.double(0.02),
    #   Vertex radius cut >
    #   (UNUSED)
    rVtxCut = cms.double(0.0),
    #   V0 decay length from primary cut >
    #   (UNUSED)
    lVtxCut = cms.double(0.0),
    #   Radial vertex significance >
    #vtxSignificance2DCut = cms.double(15.0),
    vtxSignificance2DCut = cms.double(5.0),
    #   3D vertex significance using primary vertex
    #   (UNUSED)
    vtxSignificance3DCut = cms.double(0.0),
    #   V0 mass window, Candidate mass must be within these values of
    #     the PDG mass to be stored in the collection
    # kShortMassCut = cms.double(0.07),
    kShortMassCut = cms.double(0.08),
    lambdaMassCut = cms.double(0.05),
    #   Mass window cut using normalized mass (mass / massError)
    #   (UNUSED)
    kShortNormalizedMassCut = cms.double(0.0),
    lambdaNormalizedMassCut = cms.double(0.0),
    # We check if either track has a hit inside (radially) the vertex position
    #  minus this number times the sigma of the vertex fit
    #  NOTE: Set this to -1 to disable this cut, which MUST be done
    #  if you want to run V0Producer on the AOD track collection!
    #innerHitPosCut = cms.double(4.)
    innerHitPosCut = cms.double(-1)
)

 
process.ntuple = cms.EDAnalyzer(
    'BToKstarMuMu',
    FileName = cms.string("BToKstarMuMu.root"),
    SaveGenInfo = cms.untracked.bool(False),
    GenParticlesLabel = cms.InputTag("genParticles"),
    TriggerResultsLabel = cms.InputTag("TriggerResults","", 'HLT'),
    TriggerNames = cms.vstring([]),
    LastFilterNames = cms.vstring([]),
    BeamSpotLabel = cms.InputTag('offlineBeamSpot'),
    VertexLabel = cms.InputTag('offlinePrimaryVertices'),
    MuonLabel = cms.InputTag('cleanPatMuonsTriggerMatch'),

    MuonMinPt = cms.untracked.double(3.0), # [GeV]
    MuonMaxEta = cms.untracked.double(2.2),  
    TrkMaxDcaBs = cms.untracked.double(2.0), # [cm]
    TrkMaxR = cms.untracked.double(110.0), # [cm]
    TrkMaxZ = cms.untracked.double(280.0), # [cm]
    MuMuMaxDca = cms.untracked.double(0.5), # [cm]
    MuMuMinVtxCl = cms.untracked.double(0.05), 
    MuMuMinPt = cms.untracked.double(6.9), # [GeV/c]
    MuMuMinInvMass = cms.untracked.double(1.0), # [GeV/c2]
    MuMuMaxInvMass = cms.untracked.double(4.8), # [GeV/c2]
    MuMuMinLxySigmaBs = cms.untracked.double(3.0), 
    MuMuMinCosAlphaBs = cms.untracked.double(0.9),

    MuonMass = cms.untracked.double(0.10565837), 
    MuonMassErr = cms.untracked.double(0.10565837*1e-6), 

    #KshortLabel = cms.InputTag('generalV0Candidates:Kshort'),
    KshortLabel = cms.InputTag('localV0Candidates:Kshort'),
    TrackLabel = cms.InputTag('cleanPatTrackCands'), 
    PionMass = cms.untracked.double(0.13957018), 
    PionMassErr = cms.untracked.double(0.13957018*1e-6),
    KshortMass = cms.untracked.double(0.497614), 
    KshortMassErr = cms.untracked.double(0.000024),
    KstarChargedTrackMinPt = cms.untracked.double(0.0), # [GeV/c]
    KstarMinMass = cms.untracked.double(0.49), # [GeV/c2] K*+ mass = 891.66 +- 0.26 MeV 
    KstarMaxMass = cms.untracked.double(1.29), # [GeV/c2] K*0 mass = 895.94 +- 0.26 MeV

    BMaxMass = cms.untracked.double(10.0), # [GeV/c2]n B+ mass = 5279 MeV 
    BuMass = cms.untracked.double(5.27925),
    
    B3MinMass = cms.untracked.double(4.5), 
    B3MaxMass = cms.untracked.double(6.0), 
    B3MinLsBs = cms.untracked.double(0.5),
    )


# Remove not used from PAT 
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

process.p = cms.Path(process.patDefaultSequence * process.localV0Candidates * process.ntuple)
