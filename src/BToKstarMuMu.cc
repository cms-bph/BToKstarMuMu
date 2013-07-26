// -*- C++ -*-
//
// Package:    BToKstarMuMu
// Class:      BToKstarMuMu
// 
/**\class BToKstarMuMu BToKstarMuMu.cc BphAna/BToKstarMuMu/src/BToKstarMuMu.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Xin Shi <Xin.Shi@cern.ch>
//         Created:  Mon Jan 28 11:06:56 CET 2013
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TH1.h>

using namespace std;

// 
// Constans
// 

const int MUONMINUS_PDG_ID = 13; 
const int PIONPLUS_PDG_ID = 211; 
const int KSHORTZERO_PDG_ID = 310; 
const int KSTARPLUS_PDG_ID = 323; 
const int BPLUS_PDG_ID = 521; 
const double PI = 3.141592653589793;

// Structures 
struct HistArgs{
  char name[128];
  char title[128];
  int n_bins;
  double x_min;
  double x_max;
};

enum HistName{
  h_mumumass,
  h_kshortmass,
  kHistNameSize
};

// Global hist args

HistArgs hist_args[kHistNameSize] = {
  // name, title, n_bins, x_min, x_max  
  {"h_mumumass", "#mu^{+}#mu^{-} invariant mass; M(#mu^{+}#mu^{-}) [GeV]", 100, 2, 4},
  {"h_kshortmass", "Kshort mass; M(Kshort) [GeV]", 100, 0.2, 0.8},
};

// Define histograms 
TH1F *BToKstarMuMuFigures[kHistNameSize];

//
// class declaration
//

class BToKstarMuMu : public edm::EDAnalyzer {
public:
  explicit BToKstarMuMu(const edm::ParameterSet&);
  ~BToKstarMuMu();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& );
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  
  void buildBuToKstarMuMu(const edm::Event &); 

  void computeLS (double, double, double, double, double, double, double, double, double, 
		  double, double, double, double, double, double, double, double, double, 
		  double* deltaD, double* deltaDErr); 

  void computeCosAlpha (double Vx, double Vy, double Vz, double Wx, double Wy, double Wz,
			double VxErr2, double VyErr2, double VzErr2, double VxyCov,
			double VxzCov, double VyzCov, double WxErr2, double WyErr2,
			double WzErr2, double WxyCov, double WxzCov, double WyzCov,
			double* cosAlpha, double* cosAlphaErr); 
  
  void computeCtau(RefCountedKinematicTree, double &, double &);
  double computeEta (double, double, double); 
  double computePhi (double, double, double); 
  double computeEtaPhiDistance (double, double, double,	double, double, double); 
  void clearVariables(); 

  bool hasBeamSpot(const edm::Event&);

  bool hasGoodClosestApproachTracks (const reco::TransientTrack,
				     const reco::TransientTrack, double &);

  bool hasGoodKshortVertex(const vector<reco::TrackRef>, RefCountedKinematicTree &); 
  bool hasGoodKshortVertexMKC(const vector<reco::TrackRef>, RefCountedKinematicTree &); 

  bool hasGoodDimuonKshortMass(const reco::TrackRef, const reco::TrackRef, 
			       const reco::VertexCompositeCandidate, double &); 
  bool hasGoodTrackDcaBs (const reco::TransientTrack, double &, double &); 
  bool hasGoodTrackDcaPoint (const reco::TransientTrack, const GlobalPoint, 
			     double, double &, double &);

  bool hasGoodKstarChargedMass(RefCountedKinematicTree); 

  bool hasGoodBuMass(RefCountedKinematicTree); 

  bool hasGoodBuVertex(const reco::TrackRef, const reco::TrackRef, 
		       const vector<reco::TrackRef>, const reco::TrackRef, 
		       RefCountedKinematicTree &, RefCountedKinematicTree &); 
  
  bool hasGoodMuMuVertex (const reco::TransientTrack, const reco::TransientTrack,
			  reco::TransientTrack &, reco::TransientTrack &, 
			  double &, double &, double &, double &, double &);

  bool hasGoodPionTrack(const edm::Event&, const pat::GenericParticle);

  bool hasPrimaryVertex(const edm::Event &); 

  void hltReport(const edm::Event&);

  bool isGenKstarCharged(const reco::Candidate *);
  bool isGenKshort(const reco::Candidate *);
  bool isGenMuonP(const reco::Candidate *);

  bool matchMuonTrack (const edm::Event&, const reco::TrackRef);
  bool matchMuonTracks (const edm::Event&, const vector<reco::TrackRef>);
  bool matchKshortTrack (const edm::Event&, const reco::TrackRef);
  bool matchPrimaryVertexTracks (); 

  void saveBuToKstarMuMu(RefCountedKinematicTree); 
  void saveBuVertex(RefCountedKinematicTree); 
  void saveBuCosAlpha(RefCountedKinematicTree); 
  void saveBuLsig(RefCountedKinematicTree);
  void saveBuCtau(RefCountedKinematicTree); 

  void saveGenInfo(const edm::Event&); 
  void saveKshortVariables(RefCountedKinematicTree); 

  void saveSoftMuonVariables(pat::Muon, pat::Muon, reco::TrackRef, reco::TrackRef); 
  void saveDimuVariables(double, double, double, double, double, double, double,
			 double, double, double, double);
  void saveMuonTriggerMatches(const pat::Muon, const pat::Muon); 
  void saveTruthMatch(const edm::Event& iEvent); 

  // ----------member data ---------------------------

  // --- input from python file --- 
  string FileName_; 
  bool SaveGenInfo_; 
  edm::InputTag GenParticlesLabel_;
  double TruthMatchMuonMaxR_; 
  double TruthMatchPionMaxR_; 
  double TruthMatchKsMaxVtx_; 

  edm::InputTag TriggerResultsLabel_;
  vector<string> TriggerNames_, LastFilterNames_;
  map<string, string> mapTriggerToLastFilter_;
  edm::InputTag BeamSpotLabel_;
  edm::InputTag VertexLabel_;

  // Muon 
  edm::InputTag MuonLabel_;
  double MuonMinPt_; 
  double MuonMaxEta_; 
  double TrkMaxDcaBs_; 
  double TrkMaxR_;
  double TrkMaxZ_; 
  double MuMuMaxDca_; 
  double MuMuMinVtxCl_; 
  double MuMuMinPt_; 
  double MuMuMinInvMass_; 
  double MuMuMaxInvMass_; 
  double MuMuMinLxySigmaBs_; 
  double MuMuMinCosAlphaBs_; 
  ParticleMass MuonMass_; 
  float MuonMassErr_; 
  // pat::CompositeCandidateCollection DimuonCandidates_; 

  // Kshort 
  edm::InputTag KshortLabel_;
  edm::InputTag TrackLabel_;
  ParticleMass PionMass_; 
  float PionMassErr_; 
  ParticleMass KshortMass_; 
  float KshortMassErr_; 
  // pat::CompositeCandidateCollection KshortCandidates_; 

  double DimukshortMinMass_, DimukshortMaxMass_; 
  
  // Kstar
  double KstarChargedTrackMinPt_; 
  double KstarMinMass_, KstarMaxMass_; 
  // pat::CompositeCandidateCollection KstarChargedCandidates_; 

  // B meson
  // pat::CompositeCandidateCollection BuCandidates_;  
  //   BuToPiMuMuCandidates_; 
  double BMaxMass_; 
  double BuMass_; 
  

  // Across the event 
  reco::BeamSpot beamSpot_;  
  edm::ESHandle<MagneticField> bFieldHandle_;
  reco::Vertex primaryVertex_;


  // ---- Root Variables ---- 
  TFile* fout_;
  TTree* tree_;
  
  unsigned int run, event, lumiblock, nprivtx; 
  vector<string> *triggernames;
  vector<int> *triggerprescales;
  
  // dimuon 
  vector<double> *mumdcabs, *mumdcabserr, *mumpx, *mumpy, *mumpz; 
  vector<double> *mummass, *mummasserr; 
  vector<double> *mupdcabs, *mupdcabserr, *muppx, *muppy, *muppz; 
  vector<double> *mupmass, *mupmasserr; 
  vector<double> *mumudca, *mumuvtxcl, *mumulsbs, *mumulsbserr;  
  vector<double> *mumucosalphabs, *mumucosalphabserr; 

  // soft muon variables 
  vector<bool> *mumisgoodmuon, *mupisgoodmuon ; 
  vector<int> *mumnpixhits, *mupnpixhits, *mumnpixlayers, *mupnpixlayers; 
  vector<int> *mumntrkhits, *mupntrkhits, *mumntrklayers, *mupntrklayers; 
  vector<double> *mumnormchi2, *mupnormchi2; 
  vector<double> *mumdxyvtx, *mupdxyvtx, *mumdzvtx, *mupdzvtx; 
  vector<string> *mumtriglastfilter, *muptriglastfilter; 
  vector<double> *mumpt, *muppt, *mumeta, *mupeta; 
  
  // pion track 
  vector<int> *trkchg; // +1 for pi+, -1 for pi-
  vector<double> *trkpx, *trkpy, *trkpz, *trkmass, *trkmasserr; 


  // kshort 
  vector<double> *pimpx, *pimpy, *pimpz, *pimmass, *pimd0, *pimd0err; 
  vector<double> *pippx, *pippy, *pippz, *pipmass, *pipd0, *pipd0err; 
  vector<double> *kspx, *kspy, *kspz, *ksmass, *ksmasserr; 
  vector<double> *ksvtxx, *ksvtxy, *ksvtxz, *ksvtxcl, *kslsbs, *kslsbserr; 
  vector<double> *dimuksmass; 

  // B+ and B- 
  vector<int>    *bchg; // +1 for b+, -1 for b-
  vector<double> *bpx, *bpxerr, *bpy, *bpyerr, *bpz, *bpzerr, *bmass, *bmasserr; 
  vector<double> *bvtxcl, *bvtxx, *bvtxxerr, *bvtxy, *bvtxyerr, *bvtxz, *bvtxzerr; 
  vector<double> *bcosalphabs, *bcosalphabserr, *blsbs, *blsbserr, *bctau, *bctauerr; 
  

  // For MC
  int genbchg; // +1 for b+, -1 for b-
  double genbpx, genbpy, genbpz;
  double genkstpx, genkstpy, genkstpz;
  double genkspx, genkspy, genkspz;
  double genksvtxx, genksvtxy, genksvtxz; 

  int gentrkchg; 
  double gentrkpx, gentrkpy, gentrkpz;
  double genmumpx, genmumpy, genmumpz;
  double genmuppx, genmuppy, genmuppz;
  double genpippx, genpippy, genpippz;
  double genpimpx, genpimpy, genpimpz;

  vector<bool> *istruemum, *istruemup, *istrueks, *istruetrk, *istruebu;  
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
BToKstarMuMu::BToKstarMuMu(const edm::ParameterSet& iConfig):
  FileName_(iConfig.getParameter<string>("FileName")),
  SaveGenInfo_(iConfig.getUntrackedParameter<bool>("SaveGenInfo")),
  GenParticlesLabel_(iConfig.getParameter<edm::InputTag>("GenParticlesLabel")),
  TruthMatchMuonMaxR_(iConfig.getUntrackedParameter<double>("TruthMatchMuonMaxR")),
  TruthMatchPionMaxR_(iConfig.getUntrackedParameter<double>("TruthMatchPionMaxR")),
  TruthMatchKsMaxVtx_(iConfig.getUntrackedParameter<double>("TruthMatchKsMaxVtx")),
  TriggerResultsLabel_(iConfig.getParameter<edm::InputTag>("TriggerResultsLabel")),
  TriggerNames_(iConfig.getParameter< vector<string> >("TriggerNames")),
  LastFilterNames_(iConfig.getParameter< vector<string> >("LastFilterNames")),
  BeamSpotLabel_(iConfig.getParameter<edm::InputTag>("BeamSpotLabel")),
  VertexLabel_(iConfig.getParameter<edm::InputTag>("VertexLabel")),
  MuonLabel_(iConfig.getParameter<edm::InputTag>("MuonLabel")),
  MuonMinPt_(iConfig.getUntrackedParameter<double>("MuonMinPt")),
  MuonMaxEta_(iConfig.getUntrackedParameter<double>("MuonMaxEta")),
  TrkMaxDcaBs_(iConfig.getUntrackedParameter<double>("TrkMaxDcaBs")),
  TrkMaxR_(iConfig.getUntrackedParameter<double>("TrkMaxR")),
  TrkMaxZ_(iConfig.getUntrackedParameter<double>("TrkMaxZ")),
  MuMuMaxDca_(iConfig.getUntrackedParameter<double>("MuMuMaxDca")),
  MuMuMinVtxCl_(iConfig.getUntrackedParameter<double>("MuMuMinVtxCl")),
  MuMuMinPt_(iConfig.getUntrackedParameter<double>("MuMuMinPt")),
  MuMuMinInvMass_(iConfig.getUntrackedParameter<double>("MuMuMinInvMass")),
  MuMuMaxInvMass_(iConfig.getUntrackedParameter<double>("MuMuMaxInvMass")),
  MuMuMinLxySigmaBs_(iConfig.getUntrackedParameter<double>("MuMuMinLxySigmaBs")), 
  MuMuMinCosAlphaBs_(iConfig.getUntrackedParameter<double>("MuMuMinCosAlphaBs")), 
  MuonMass_(iConfig.getUntrackedParameter<double>("MuonMass")),
  MuonMassErr_(iConfig.getUntrackedParameter<double>("MuonMassErr")),
  KshortLabel_(iConfig.getParameter<edm::InputTag>("KshortLabel")),
  TrackLabel_(iConfig.getParameter<edm::InputTag>("TrackLabel")),
  PionMass_(iConfig.getUntrackedParameter<double>("PionMass")),
  PionMassErr_(iConfig.getUntrackedParameter<double>("PionMassErr")),
  KshortMass_(iConfig.getUntrackedParameter<double>("KshortMass")),
  KshortMassErr_(iConfig.getUntrackedParameter<double>("KshortMassErr")),
  DimukshortMinMass_(iConfig.getUntrackedParameter<double>("DimukshortMinMass")),
  DimukshortMaxMass_(iConfig.getUntrackedParameter<double>("DimukshortMaxMass")),

  KstarChargedTrackMinPt_(iConfig.getUntrackedParameter<double>("KstarChargedTrackMinPt")),
  KstarMinMass_(iConfig.getUntrackedParameter<double>("KstarMinMass")),
  KstarMaxMass_(iConfig.getUntrackedParameter<double>("KstarMaxMass")),
  BMaxMass_(iConfig.getUntrackedParameter<double>("BMaxMass")),
  BuMass_(iConfig.getUntrackedParameter<double>("BuMass")),

  tree_(0), 
  triggernames(0), triggerprescales(0), 
  mumdcabs(0), mumdcabserr(0), mumpx(0), mumpy(0), mumpz(0), mummass(0), mummasserr(0), 
  mupdcabs(0),  mupdcabserr(0), muppx(0),  muppy(0),  muppz(0), mupmass(0), mupmasserr(0), 
  mumudca(0),  mumuvtxcl(0),  mumulsbs(0),  mumulsbserr(0), 
  mumucosalphabs(0),  mumucosalphabserr(0), 
  mumisgoodmuon(0), mupisgoodmuon(0), 
  mumnpixhits(0), mupnpixhits(0), mumnpixlayers(0), mupnpixlayers(0), 
  mumntrkhits(0), mupntrkhits(0), mumntrklayers(0), mupntrklayers(0), 
  mumnormchi2(0), mupnormchi2(0), mumdxyvtx(0), mupdxyvtx(0), mumdzvtx(0), mupdzvtx(0), 
  mumtriglastfilter(0), muptriglastfilter(0), 
  mumpt(0), muppt(0), mumeta(0), mupeta(0), 

  trkchg(0), 
  trkpx(0), trkpy(0), trkpz(0), trkmass(0), trkmasserr(0), 

  pimpx(0), pimpy(0), pimpz(0), pimmass(0), pimd0(0), pimd0err(0), 
  pippx(0), pippy(0), pippz(0), pipmass(0), pipd0(0), pipd0err(0), 
  kspx(0), kspy(0), kspz(0), ksmass(0), ksmasserr(0), 
  ksvtxx(0), ksvtxy(0), ksvtxz(0), ksvtxcl(0), kslsbs(0), kslsbserr(0), 
  dimuksmass(0), 
  bchg(0), bpx(0), bpxerr(0), bpy(0), bpyerr(0), bpz(0), bpzerr(0), bmass(0), bmasserr(0), 
  bvtxcl(0), bvtxx(0), bvtxxerr(0), bvtxy(0), bvtxyerr(0), bvtxz(0), bvtxzerr(0), 
  bcosalphabs(0), bcosalphabserr(0), blsbs(0), blsbserr(0), bctau(0), bctauerr(0), 


  genbchg(0), 
  genbpx(0), genbpy(0), genbpz(0), 
  genkstpx(0), genkstpy(0), genkstpz(0), 
  genkspx(0), genkspy(0), genkspz(0), 
  genksvtxx(0), genksvtxy(0), genksvtxz(0), 
  gentrkchg(0), gentrkpx(0), gentrkpy(0), gentrkpz(0), 
  genmumpx(0), genmumpy(0), genmumpz(0),
  genmuppx(0), genmuppy(0), genmuppz(0), 
  genpippx(0), genpippy(0), genpippz(0), 
  genpimpx(0), genpimpy(0), genpimpz(0), 
  istruemum(0), istruemup(0), istrueks(0), istruetrk(0), istruebu(0) 
  
{ 
  //now do what ever initialization is needed
  
  assert(TriggerNames_.size() == LastFilterNames_.size());
  for (size_t i = 0; i < TriggerNames_.size(); ++i)
    mapTriggerToLastFilter_[TriggerNames_[i]] = LastFilterNames_[i];
  
}


BToKstarMuMu::~BToKstarMuMu()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
BToKstarMuMu::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  clearVariables(); 

  run = iEvent.id().run() ;
  event = iEvent.id().event() ;
  lumiblock = iEvent.luminosityBlock(); 

  hltReport(iEvent);

  if ( hasBeamSpot(iEvent) ) {
    iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle_);      
    
    if ( bFieldHandle_.isValid() && hasPrimaryVertex(iEvent) ) { 
      
      buildBuToKstarMuMu(iEvent) ;  

      if (SaveGenInfo_) {
	saveGenInfo(iEvent);
	saveTruthMatch(iEvent);
      }
      
      tree_->Fill();
    }

  }
  
  clearVariables(); 
}


// ------------ method called once each job just before starting event loop  ------------
void 
BToKstarMuMu::beginJob()
{
  fout_ = new TFile(FileName_.c_str(), "RECREATE");
  fout_->cd();

  for(int i=0; i<kHistNameSize; i++) {
    BToKstarMuMuFigures[i] = new TH1F(hist_args[i].name, hist_args[i].title,
				       hist_args[i].n_bins,
				       hist_args[i].x_min, hist_args[i].x_max);
  }
  

  tree_ = new TTree ("tree", "BToKstarMuMu");

  tree_->Branch("run", &run, "run/i");
  tree_->Branch("event", &event, "event/i");
  tree_->Branch("lumiblock", &lumiblock, "lumiblock/i");
  tree_->Branch("nprivtx", &nprivtx, "nprivtx/i");
  tree_->Branch("triggernames", &triggernames);
  tree_->Branch("triggerprescales", &triggerprescales);

  tree_->Branch("mumdcabs", &mumdcabs);
  tree_->Branch("mumdcabserr", &mumdcabserr);
  tree_->Branch("mumpx", &mumpx);
  tree_->Branch("mumpy", &mumpy);
  tree_->Branch("mumpz", &mumpz);
  tree_->Branch("mummass", &mummass); 
  tree_->Branch("mummasserr", &mummasserr); 
  tree_->Branch("mupdcabs", &mupdcabs);
  tree_->Branch("mupdcabserr", &mupdcabserr);
  tree_->Branch("muppx", &muppx);
  tree_->Branch("muppy", &muppy);
  tree_->Branch("muppz", &muppz);
  tree_->Branch("mupmass", &mupmass); 
  tree_->Branch("mupmasserr", &mupmasserr); 
  tree_->Branch("mumudca", &mumudca);
  tree_->Branch("mumuvtxcl", &mumuvtxcl);
  tree_->Branch("mumulsbs", &mumulsbs);
  tree_->Branch("mumulsbserr", &mumulsbserr);
  tree_->Branch("mumucosalphabs", &mumucosalphabs);
  tree_->Branch("mumucosalphabserr", &mumucosalphabserr);
  tree_->Branch("mumisgoodmuon", &mumisgoodmuon);
  tree_->Branch("mupisgoodmuon", &mupisgoodmuon);
  tree_->Branch("mumnpixhits", &mumnpixhits);
  tree_->Branch("mupnpixhits", &mupnpixhits);
  tree_->Branch("mumnpixlayers", &mumnpixlayers);
  tree_->Branch("mupnpixlayers", &mupnpixlayers);
  tree_->Branch("mumntrkhits", &mumntrkhits);
  tree_->Branch("mupntrkhits", &mupntrkhits);
  tree_->Branch("mumntrklayers", &mumntrklayers);
  tree_->Branch("mupntrklayers", &mupntrklayers);
  tree_->Branch("mumnormchi2", &mumnormchi2);
  tree_->Branch("mupnormchi2", &mupnormchi2);
  tree_->Branch("mumdxyvtx", &mumdxyvtx);
  tree_->Branch("mupdxyvtx", &mupdxyvtx);
  tree_->Branch("mumdzvtx", &mumdzvtx);
  tree_->Branch("mupdzvtx", &mupdzvtx);
  tree_->Branch("mumtriglastfilter", &mumtriglastfilter);
  tree_->Branch("muptriglastfilter", &muptriglastfilter);
  tree_->Branch("mumpt", &mumpt);
  tree_->Branch("muppt", &muppt);
  tree_->Branch("mumeta", &mumeta);
  tree_->Branch("mupeta", &mupeta);

  tree_->Branch("trkchg", &trkchg);
  tree_->Branch("trkpx", &trkpx);
  tree_->Branch("trkpy", &trkpy);
  tree_->Branch("trkpz", &trkpz);
  tree_->Branch("trkmass", &trkmass);
  tree_->Branch("trkmasserr", &trkmasserr);

  tree_->Branch("pimpx", &pimpx);
  tree_->Branch("pimpy", &pimpy);
  tree_->Branch("pimpz", &pimpz);
  tree_->Branch("pimmass", &pimmass);
  tree_->Branch("pimd0", &pimd0);
  tree_->Branch("pimd0err", &pimd0err);
  tree_->Branch("pippx", &pippx);
  tree_->Branch("pippy", &pippy);
  tree_->Branch("pippz", &pippz);
  tree_->Branch("pipmass", &pipmass);
  tree_->Branch("pipd0", &pipd0);
  tree_->Branch("pipd0err", &pipd0err);
  tree_->Branch("kspx", &kspx);
  tree_->Branch("kspy", &kspy);
  tree_->Branch("kspz", &kspz);
  tree_->Branch("ksmass", &ksmass);
  tree_->Branch("ksmasserr", &ksmasserr);
  tree_->Branch("ksvtxx", &ksvtxx);
  tree_->Branch("ksvtxy", &ksvtxy);
  tree_->Branch("ksvtxz", &ksvtxz);
  tree_->Branch("ksvtxcl", &ksvtxcl);
  tree_->Branch("kslsbs", &kslsbs);
  tree_->Branch("kslsbserr", &kslsbserr);
  tree_->Branch("dimuksmass", &dimuksmass);
    
  tree_->Branch("bchg", &bchg);
  tree_->Branch("bpx", &bpx);
  tree_->Branch("bpxerr", &bpxerr);
  tree_->Branch("bpy", &bpy);
  tree_->Branch("bpyerr", &bpyerr);
  tree_->Branch("bpz", &bpz);
  tree_->Branch("bpzerr", &bpzerr);
  tree_->Branch("bmass", &bmass);
  tree_->Branch("bmasserr", &bmasserr);
  tree_->Branch("bvtxcl", &bvtxcl);
  tree_->Branch("bvtxx", &bvtxx);
  tree_->Branch("bvtxxerr", &bvtxxerr);
  tree_->Branch("bvtxy", &bvtxy);
  tree_->Branch("bvtxyerr", &bvtxyerr);
  tree_->Branch("bvtxz", &bvtxz);
  tree_->Branch("bvtxzerr", &bvtxzerr);
  tree_->Branch("bcosalphabs", &bcosalphabs);
  tree_->Branch("bcosalphabserr", &bcosalphabserr);
  tree_->Branch("blsbs", &blsbs);
  tree_->Branch("blsbserr", &blsbserr);
  tree_->Branch("bctau", &bctau);
  tree_->Branch("bctauerr", &bctauerr);


  if (SaveGenInfo_) {
    tree_->Branch("genbchg",     &genbchg    , "genbchg/I"   );
    tree_->Branch("genbpx",      &genbpx     , "genbpx/D"    );
    tree_->Branch("genbpy",      &genbpy     , "genbpy/D"    );
    tree_->Branch("genbpz",      &genbpz     , "genbpz/D"    );
    tree_->Branch("genkstpx",    &genkstpx   , "genkstpx/D"  );
    tree_->Branch("genkstpy",    &genkstpy   , "genkstpy/D"  );
    tree_->Branch("genkstpz",    &genkstpz   , "genkstpz/D"  );
    tree_->Branch("genkspx",     &genkspx    , "genkspx/D"   );
    tree_->Branch("genkspy",     &genkspy    , "genkspy/D"   );
    tree_->Branch("genkspz",     &genkspz    , "genkspz/D"   );
    tree_->Branch("genksvtxx",     &genksvtxx    , "genksvtxx/D"   );
    tree_->Branch("genksvtxy",     &genksvtxy    , "genksvtxy/D"   );
    tree_->Branch("genksvtxz",     &genksvtxz    , "genksvtxz/D"   );
    tree_->Branch("gentrkchg",     &gentrkchg    , "gentrkchg/I"  );
    tree_->Branch("gentrkpx",     &gentrkpx    , "gentrkpx/D"   );
    tree_->Branch("gentrkpy",     &gentrkpy    , "gentrkpy/D"   );
    tree_->Branch("gentrkpz",     &gentrkpz    , "gentrkpz/D"   );
    tree_->Branch("genmumpx",    &genmumpx   , "genmumpx/D"  );
    tree_->Branch("genmumpy",    &genmumpy   , "genmumpy/D"  );
    tree_->Branch("genmumpz",    &genmumpz   , "genmumpz/D"  );
    tree_->Branch("genmuppx",    &genmuppx   , "genmuppx/D"  );
    tree_->Branch("genmuppy",    &genmuppy   , "genmuppy/D"  );
    tree_->Branch("genmuppz",    &genmuppz   , "genmuppz/D"  );
    tree_->Branch("genpippx",  &genpippx , "genpippx/D");
    tree_->Branch("genpippy",  &genpippy , "genpippy/D");
    tree_->Branch("genpippz",  &genpippz , "genpippz/D");
    tree_->Branch("genpimpx",  &genpimpx , "genpimpx/D");
    tree_->Branch("genpimpy",  &genpimpy , "genpimpy/D");
    tree_->Branch("genpimpz",  &genpimpz , "genpimpz/D");
    tree_->Branch("istruemum",  &istruemum );
    tree_->Branch("istruemup",  &istruemup );
    tree_->Branch("istrueks",   &istrueks  );
    tree_->Branch("istruetrk",   &istruetrk  );
    tree_->Branch("istruebu",   &istruebu  );

  } 


}

// ------------ method called once each job just after ending the event loop  ------------
void 
BToKstarMuMu::endJob() 
{
  fout_->cd();
  tree_->Write();
  
  for(int i = 0; i < kHistNameSize; i++) {
    BToKstarMuMuFigures[i]->Write();
    BToKstarMuMuFigures[i]->Delete();
  }
  fout_->Close();
}

// ------------ method called when starting to processes a run  ------------
void 
BToKstarMuMu::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
BToKstarMuMu::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
BToKstarMuMu::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
BToKstarMuMu::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
BToKstarMuMu::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}



void 
BToKstarMuMu::clearVariables(){
  run = 0; 
  event = 0;
  lumiblock = 0;
  nprivtx = 0; 
  triggernames->clear();
  triggerprescales->clear();
  mumdcabs->clear();  mumdcabserr->clear();  mumpx->clear();   mumpy->clear();   mumpz->clear(); 
  mummass->clear(); mummasserr->clear(); 
  mupdcabs->clear();  mupdcabserr->clear();  muppx->clear();   muppy->clear();   muppz->clear(); 
  mupmass->clear(); mupmasserr->clear(); 
  mumudca->clear();  mumuvtxcl->clear();   mumulsbs->clear();  mumulsbserr->clear(); 
  mumucosalphabs->clear();  mumucosalphabserr->clear();  
  mumisgoodmuon->clear();  mupisgoodmuon->clear(); 
  mumnpixhits->clear();  mupnpixhits->clear();  mumnpixlayers->clear();  mupnpixlayers->clear(); 
  mumntrkhits->clear();  mupntrkhits->clear();  mumntrklayers->clear();  mupntrklayers->clear(); 
  mumnormchi2->clear(); mupnormchi2->clear(); mumdxyvtx->clear(); mupdxyvtx->clear(); 
  mumdzvtx->clear(); mupdzvtx->clear();  mumtriglastfilter->clear(); muptriglastfilter->clear(); 
  mumpt->clear(); muppt->clear(); 
  mumeta->clear(); mupeta->clear(); 

  trkchg->clear(); trkpx->clear(); trkpy->clear(); trkpz->clear(); 
  trkmass->clear(); trkmasserr->clear(); 

  pimpx->clear(); pimpy->clear(); pimpz->clear(); pimmass->clear(); pimd0->clear(); pimd0err->clear();
  pippx->clear(); pippy->clear(); pippz->clear(); pipmass->clear(); pipd0->clear(); pipd0err->clear();
  kspx->clear(); kspy->clear(); kspz->clear(); ksmass->clear(); ksmasserr->clear(); 
  ksvtxx->clear(); ksvtxy->clear(); ksvtxz->clear(); ksvtxcl->clear();
  kslsbs->clear(); kslsbserr->clear(); dimuksmass->clear();  
  // DimuonCandidates_.clear();  KshortCandidates_.clear(); 
  // KstarChargedCandidates_.clear(); BuCandidates_.clear();

  bpx->clear(); bpxerr->clear(); bpy->clear();  bpyerr->clear(); bpz->clear(); bpzerr->clear(); 
  bchg->clear(); bmass->clear(); bmasserr->clear(); 
  bvtxcl->clear(); bvtxx->clear(); bvtxxerr->clear(); bvtxy->clear(); bvtxyerr->clear();
  bvtxz->clear(); bvtxzerr->clear(); bcosalphabs->clear(); bcosalphabserr->clear(); 
  blsbs->clear(); blsbserr->clear();  bctau->clear(); bctauerr->clear(); 

  if (SaveGenInfo_) {
    genbchg = 0; genbpx = 0;    genbpy = 0;    genbpz = 0; 
    genkstpx = 0;  genkstpy = 0;  genkstpz = 0; 
    genkspx = 0;   genkspy = 0;   genkspz = 0; 
    genksvtxx = 0;   genksvtxy = 0;   genksvtxz = 0; 
    gentrkchg = 0; gentrkpx = 0;   gentrkpy = 0;   gentrkpz = 0; 
    genmumpx = 0;  genmumpy = 0;  genmumpz = 0; 
    genmuppx = 0;  genmuppy = 0;  genmuppz = 0; 
    genpippx = 0;  genpippy = 0;  genpippz = 0; 
    genpimpx = 0;  genpimpy = 0;  genpimpz = 0; 

    istruemum->clear(); istruemup->clear(); istrueks->clear();
    istruetrk->clear(); istruebu->clear(); 
  }

}


void
BToKstarMuMu::hltReport(const edm::Event& iEvent)
{
  edm::Handle<edm::TriggerResults> hltTriggerResults;
  try {iEvent.getByLabel( TriggerResultsLabel_, hltTriggerResults ); }
  catch ( ... ) { edm::LogInfo("myHLT") << __LINE__ << " : couldn't get handle on HLT Trigger" ; }
  
  // Save trigger info
  HLTConfigProvider hltConfig_;
  if (hltTriggerResults.isValid()) {
    const edm::TriggerNames& triggerNames_ = iEvent.triggerNames(*hltTriggerResults);

    for (unsigned int itrig = 0; itrig < hltTriggerResults->size(); itrig++){

      // Only consider the triggered case. 
      if ((*hltTriggerResults)[itrig].accept() == 1){ 

	string triggername = triggerNames_.triggerName(itrig);	
	int triggerprescale = hltConfig_.prescaleValue(itrig, triggername);

	// Loop over our interested HLT trigger names to find if this event contains. 
	for (unsigned int it=0; it<TriggerNames_.size(); it++){
	  //if (TriggerNames_[it] == triggername)	{
	  if (triggername.find(TriggerNames_[it]) != string::npos) {
	    // triggernames->push_back(triggername); 
	    // save the no versioned case 
	    triggernames->push_back(TriggerNames_[it]); 
	    triggerprescales->push_back(triggerprescale); 

	  }}}}}
}


bool
BToKstarMuMu::hasBeamSpot(const edm::Event& iEvent)
{
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByLabel(BeamSpotLabel_, beamSpotHandle);
  
  if ( ! beamSpotHandle.isValid() ) {
    edm::LogError("myBeam") << "No beam spot available from EventSetup" ; 
    return false; 
  }
  
  beamSpot_ = *beamSpotHandle; 
  return true; 
}


bool
BToKstarMuMu::hasPrimaryVertex(const edm::Event& iEvent)
{
  edm::Handle<reco::VertexCollection> recVtxs;
  iEvent.getByLabel(VertexLabel_, recVtxs);
  nprivtx = recVtxs->size(); 

  for (std::vector<reco::Vertex>::const_iterator iVertex = recVtxs->begin();
       iVertex != recVtxs->end(); iVertex++) { 
    primaryVertex_ = *(iVertex); 
    if (primaryVertex_.isValid()) break; 
  }

  if (!primaryVertex_.isValid()) return false; 
 
  return true; 
}



bool 
BToKstarMuMu::hasGoodPionTrack(const edm::Event& iEvent, const pat::GenericParticle iTrack)
{
   reco::TrackRef theTrackRef = iTrack.track(); 
   if ( theTrackRef.isNull() ) return false; 

   // veto muon tracks
   if ( matchMuonTrack(iEvent, theTrackRef) ) return false; 
   
   // veto pion tracks from Kshort
   if ( matchKshortTrack(iEvent, theTrackRef) ) return false; 
   
   // check the track kinematics
   if ( theTrackRef->pt() < KstarChargedTrackMinPt_ ) return false; 

   // compute track DCA to beam spot 
   const reco::TransientTrack theTrackTT(theTrackRef, &(*bFieldHandle_));   
   double DCAKstTrkBS, DCAKstTrkBSErr; 
   if ( ! hasGoodTrackDcaBs(theTrackTT, DCAKstTrkBS, DCAKstTrkBSErr)) return false; 

   return true;
}



void 
BToKstarMuMu::buildBuToKstarMuMu(const edm::Event& iEvent)
{
  // init variables 
  edm::Handle< vector<pat::Muon> > patMuonHandle;
  iEvent.getByLabel(MuonLabel_, patMuonHandle);
  if( patMuonHandle->size() < 2 ) return;

  edm::Handle<reco::VertexCompositeCandidateCollection> theKshorts;
  iEvent.getByLabel(KshortLabel_, theKshorts);
  if ( theKshorts->size() <= 0) return;

  edm::Handle< vector<pat::GenericParticle> >thePATTrackHandle;
  iEvent.getByLabel(TrackLabel_, thePATTrackHandle);

  double DCAmumBS, DCAmumBSErr, DCAmupBS, DCAmupBSErr, DCAmumu; 
  reco::TransientTrack refitMupTT, refitMumTT; 
  double mu_mu_vtx_cl, MuMuLSBS, MuMuLSBSErr, MuMuCosAlphaBS, MuMuCosAlphaBSErr;

  double dimu_ks_mass; 

  vector<reco::TrackRef> kshortDaughterTracks;
  int nBu = 0; 
  RefCountedKinematicTree vertexFitTree, ksVertexFitTree;


  // ---------------------------------
  // loop 1: mu-
  // ---------------------------------
  for (vector<pat::Muon>::const_iterator iMuonM = patMuonHandle->begin(); 
       iMuonM != patMuonHandle->end(); iMuonM++){
    
    reco::TrackRef muTrackm = iMuonM->innerTrack(); 
    if ( muTrackm.isNull() || 
	 (muTrackm->charge() != -1) ||
	 (muTrackm->pt() < MuonMinPt_) ||
	 (fabs(muTrackm->eta()) > MuonMaxEta_)) continue;
    
    // check mu- DCA to beam spot 
    const reco::TransientTrack muTrackmTT(muTrackm, &(*bFieldHandle_));   
    if ( ! hasGoodTrackDcaBs(muTrackmTT, DCAmumBS, DCAmumBSErr)) continue; 

    // ---------------------------------
    // loop 2: mu+ 
    // ---------------------------------
    for (vector<pat::Muon>::const_iterator iMuonP = patMuonHandle->begin(); 
	 iMuonP != patMuonHandle->end(); iMuonP++){

      reco::TrackRef muTrackp = iMuonP->innerTrack(); 
      if ( muTrackp.isNull() || 
	   (muTrackp->charge() != 1) ||
	   (muTrackp->pt() < MuonMinPt_) ||
	   (fabs(muTrackp->eta()) > MuonMaxEta_)) continue;
      
      // check mu+ DCA to beam spot 
      const reco::TransientTrack muTrackpTT(muTrackp, &(*bFieldHandle_));   
      if ( ! hasGoodTrackDcaBs(muTrackpTT, DCAmupBS, DCAmupBSErr)) continue; 
      
      // check goodness of muons closest approach and the 3D-DCA
      if ( !hasGoodClosestApproachTracks(muTrackpTT, muTrackmTT, DCAmumu) ) continue; 
      
      // check dimuon vertex 
      if ( !hasGoodMuMuVertex(muTrackpTT, muTrackmTT, refitMupTT, refitMumTT, 
			      mu_mu_vtx_cl,MuMuLSBS, MuMuLSBSErr,
			      MuMuCosAlphaBS, MuMuCosAlphaBSErr) ) continue; 
      
      // ---------------------------------
      // loop 3: kshort 
      // ---------------------------------
      for ( reco::VertexCompositeCandidateCollection::const_iterator iKshort 
	      = theKshorts->begin(); iKshort != theKshorts->end(); ++iKshort) {
	
	// check dimuon + kshort mass 
	if ( ! hasGoodDimuonKshortMass(muTrackm, muTrackp, *iKshort, dimu_ks_mass ) )
	  continue; 
	
	// check the daughter pions from Kshort is overlap with muons
	kshortDaughterTracks.push_back((dynamic_cast<const reco::RecoChargedCandidate *>
					(iKshort->daughter(0)))->track()); 
	kshortDaughterTracks.push_back((dynamic_cast<const reco::RecoChargedCandidate *>
					(iKshort->daughter(1)))->track()); 
	if ( matchMuonTracks(iEvent, kshortDaughterTracks) ) continue; 

	// ---------------------------------
	// loop 4: track 
	// ---------------------------------
	for ( vector<pat::GenericParticle>::const_iterator iTrack = thePATTrackHandle->begin();
	    iTrack != thePATTrackHandle->end(); ++iTrack ) {

	  if ( ! hasGoodPionTrack(iEvent, *iTrack)) continue; 
	  
	  reco::TrackRef pionTrack = iTrack->track(); 
	  
	  if ( ! hasGoodBuVertex(muTrackm, muTrackp, kshortDaughterTracks, pionTrack,
				 vertexFitTree, ksVertexFitTree)) continue; 
	  
	  if ( ! hasGoodKstarChargedMass(vertexFitTree) ) continue; 
	  
	  if ( ! hasGoodBuMass(vertexFitTree) ) continue; 

	  nBu++; 
	  // save the tree variables 
	  saveDimuVariables(DCAmumBS, DCAmumBSErr, DCAmupBS, DCAmupBSErr,
			    DCAmumu,  mu_mu_vtx_cl,  MuMuLSBS, MuMuLSBSErr, 
			    MuMuCosAlphaBS, MuMuCosAlphaBSErr, dimu_ks_mass); 
	  
	  saveSoftMuonVariables(*iMuonM, *iMuonP, muTrackm, muTrackp); 

	  saveKshortVariables(ksVertexFitTree); 

	  bchg->push_back(iTrack->charge()); 
	  saveBuToKstarMuMu(vertexFitTree); 
	  saveBuVertex(vertexFitTree); 
	  saveBuCosAlpha(vertexFitTree); 
	  saveBuLsig(vertexFitTree); 
	  saveBuCtau(vertexFitTree); 

	} // close track loop
      } // close kshort loop
    } // close mu+ loop
  } // close mu- loop 

  if ( nBu > 0) edm::LogInfo("myBu") << "Found " << nBu  << " Bu -> K* mu mu."; 
    
}


void 
BToKstarMuMu::computeLS (double Vx, double Vy, double Vz,
			 double Wx, double Wy, double Wz,
			 double VxErr2, double VyErr2, double VzErr2,
			 double VxyCov, double VxzCov, double VyzCov,
			 double WxErr2, double WyErr2, double WzErr2,
			 double WxyCov, double WxzCov, double WyzCov,
			 double* deltaD, double* deltaDErr)
{
  *deltaD = sqrt((Vx-Wx) * (Vx-Wx) + (Vy-Wy) * (Vy-Wy) + (Vz-Wz) * (Vz-Wz));
  if (*deltaD > 0.)
    *deltaDErr = sqrt((Vx-Wx) * (Vx-Wx) * VxErr2 +
		      (Vy-Wy) * (Vy-Wy) * VyErr2 +
		      (Vz-Wz) * (Vz-Wz) * VzErr2 +
		      
		      (Vx-Wx) * (Vy-Wy) * 2.*VxyCov +
		      (Vx-Wx) * (Vz-Wz) * 2.*VxzCov +
		      (Vy-Wy) * (Vz-Wz) * 2.*VyzCov +
		      
		      (Vx-Wx) * (Vx-Wx) * WxErr2 +
		      (Vy-Wy) * (Vy-Wy) * WyErr2 +
		      (Vz-Wz) * (Vz-Wz) * WzErr2 +
		      
		      (Vx-Wx) * (Vy-Wy) * 2.*WxyCov +
		      (Vx-Wx) * (Vz-Wz) * 2.*WxzCov +
		      (Vy-Wy) * (Vz-Wz) * 2.*WyzCov) / *deltaD;
  else *deltaDErr = 0.;
}


void 
BToKstarMuMu::computeCosAlpha (double Vx, double Vy, double Vz,
			       double Wx, double Wy, double Wz,
			       double VxErr2, double VyErr2, double VzErr2,
			       double VxyCov, double VxzCov, double VyzCov,
			       double WxErr2, double WyErr2, double WzErr2,
			       double WxyCov, double WxzCov, double WyzCov,
			       double* cosAlpha, double* cosAlphaErr)
{
  double Vnorm = sqrt(Vx*Vx + Vy*Vy + Vz*Vz);
  double Wnorm = sqrt(Wx*Wx + Wy*Wy + Wz*Wz);
  double VdotW = Vx*Wx + Vy*Wy + Vz*Wz;
  
  if ((Vnorm > 0.) && (Wnorm > 0.)) {
      *cosAlpha = VdotW / (Vnorm * Wnorm);
      *cosAlphaErr = sqrt( ((Vx*Wnorm - VdotW*Wx) * (Vx*Wnorm - VdotW*Wx) * WxErr2 +
			    (Vy*Wnorm - VdotW*Wy) * (Vy*Wnorm - VdotW*Wy) * WyErr2 +
			    (Vz*Wnorm - VdotW*Wz) * (Vz*Wnorm - VdotW*Wz) * WzErr2 +
			   
			    (Vx*Wnorm - VdotW*Wx) * (Vy*Wnorm - VdotW*Wy) * 2.*WxyCov +
			    (Vx*Wnorm - VdotW*Wx) * (Vz*Wnorm - VdotW*Wz) * 2.*WxzCov +
			    (Vy*Wnorm - VdotW*Wy) * (Vz*Wnorm - VdotW*Wz) * 2.*WyzCov) /
			   (Wnorm*Wnorm*Wnorm*Wnorm) +
			   
			   ((Wx*Vnorm - VdotW*Vx) * (Wx*Vnorm - VdotW*Vx) * VxErr2 +
			    (Wy*Vnorm - VdotW*Vy) * (Wy*Vnorm - VdotW*Vy) * VyErr2 +
			    (Wz*Vnorm - VdotW*Vz) * (Wz*Vnorm - VdotW*Vz) * VzErr2 +
			    
			    (Wx*Vnorm - VdotW*Vx) * (Wy*Vnorm - VdotW*Vy) * 2.*VxyCov +
			    (Wx*Vnorm - VdotW*Vx) * (Wz*Vnorm - VdotW*Vz) * 2.*VxzCov +
			    (Wy*Vnorm - VdotW*Vy) * (Wz*Vnorm - VdotW*Vz) * 2.*VyzCov) /
			   (Vnorm*Vnorm*Vnorm*Vnorm) ) / (Wnorm*Vnorm);
    }  else {
      *cosAlpha = 0.;
      *cosAlphaErr = 0.;
    }
}


bool 
BToKstarMuMu::hasGoodTrackDcaBs (const reco::TransientTrack muTrackTT, 
				 double &muDcaBs, double &muDcaBsErr)
{
  TrajectoryStateClosestToPoint theDCAXBS = muTrackTT.trajectoryStateClosestToPoint(GlobalPoint(beamSpot_.position().x(),beamSpot_.position().y(),beamSpot_.position().z()));
  
  if ( !theDCAXBS.isValid() )  return false; 

  muDcaBs = theDCAXBS.perigeeParameters().transverseImpactParameter();
  muDcaBsErr = theDCAXBS.perigeeError().transverseImpactParameterError();
  if ( muDcaBs > TrkMaxDcaBs_ )   return false; 
  return true; 
}


bool 
BToKstarMuMu::hasGoodTrackDcaPoint (const reco::TransientTrack track, const GlobalPoint p, 
				    double maxdca, double &dca, double &dcaerr)
{
  TrajectoryStateClosestToPoint theDCAX = track.trajectoryStateClosestToPoint(p); 
  if ( !theDCAX.isValid() ) return false; 
  
  dca = theDCAX.perigeeParameters().transverseImpactParameter();
  dcaerr = theDCAX.perigeeError().transverseImpactParameterError();
  if ( dca > maxdca ) return false; 

  return true; 
}


bool
BToKstarMuMu::hasGoodClosestApproachTracks (const reco::TransientTrack muTrackpTT, 
					    const reco::TransientTrack muTrackmTT,
					    double & DCAmumu)
{
  ClosestApproachInRPhi ClosestApp; 
  ClosestApp.calculate(muTrackpTT.initialFreeState(),muTrackmTT.initialFreeState());
  if (! ClosestApp.status())  return false;
  
  GlobalPoint XingPoint = ClosestApp.crossingPoint();
  
  if ((sqrt(XingPoint.x()*XingPoint.x() + XingPoint.y()*XingPoint.y()) > 
       TrkMaxR_) || (fabs(XingPoint.z()) > TrkMaxZ_))  return false;

  DCAmumu = ClosestApp.distance();
  if (DCAmumu > MuMuMaxDca_) return false;
 
  return true; 
}


bool 
BToKstarMuMu::hasGoodMuMuVertex (const reco::TransientTrack muTrackpTT, 
				 const reco::TransientTrack muTrackmTT, 
				 reco::TransientTrack &refitMupTT, 
				 reco::TransientTrack &refitMumTT, 
				 double & mu_mu_vtx_cl,
				 double & MuMuLSBS, double & MuMuLSBSErr, 
				 double & MuMuCosAlphaBS, double & MuMuCosAlphaBSErr)
{
  KinematicParticleFactoryFromTransientTrack partFactory;
  KinematicParticleVertexFitter PartVtxFitter;   

  vector<RefCountedKinematicParticle> muonParticles;
  double chi = 0.;
  double ndf = 0.;
  muonParticles.push_back(partFactory.particle(muTrackmTT, MuonMass_,chi,ndf,MuonMassErr_));
  muonParticles.push_back(partFactory.particle(muTrackpTT,MuonMass_,chi,ndf,MuonMassErr_));
  
  RefCountedKinematicTree mumuVertexFitTree = PartVtxFitter.fit(muonParticles);
  
  if ( !mumuVertexFitTree->isValid())  {
    // edm::LogInfo("myDimuon") << "Invalid vertex from the mu+ mu- vertex fit"; 
    return false;
  }
  
  mumuVertexFitTree->movePointerToTheTop();
  RefCountedKinematicParticle mumu_KP = mumuVertexFitTree->currentParticle();
  RefCountedKinematicVertex mumu_KV = mumuVertexFitTree->currentDecayVertex();
  
  if ( !mumu_KV->vertexIsValid()) {
    // edm::LogInfo("myDimuon") << "Invalid vertex from the mu+ mu- vertex fit";
    return false;
  }
  
  mu_mu_vtx_cl = TMath::Prob((double)mumu_KV->chiSquared(), int(rint(mumu_KV->degreesOfFreedom())));
  
  if (mu_mu_vtx_cl < MuMuMinVtxCl_)   {
    // edm::LogInfo("myDimuon") << "Bad vtx CL from mu+ mu- fit: " << mu_mu_vtx_cl; 
    return false;
  }
  // extract the re-fitted tracks
  mumuVertexFitTree->movePointerToTheTop();
  
  mumuVertexFitTree->movePointerToTheFirstChild();
  RefCountedKinematicParticle refitMum = mumuVertexFitTree->currentParticle();
  refitMumTT = refitMum->refittedTransientTrack();
  
  mumuVertexFitTree->movePointerToTheNextChild();
  RefCountedKinematicParticle refitMup = mumuVertexFitTree->currentParticle();
  refitMupTT = refitMup->refittedTransientTrack();
  
  TLorentzVector mymum, mymup, mydimu; 
  
  mymum.SetXYZM(refitMumTT.track().momentum().x(), refitMumTT.track().momentum().y(),
		refitMumTT.track().momentum().z(), MuonMass_); 
  mymup.SetXYZM(refitMupTT.track().momentum().x(), refitMupTT.track().momentum().y(),
		refitMupTT.track().momentum().z(), MuonMass_); 
  
  mydimu = mymum + mymup; 
  
  if ((mydimu.Perp() < MuMuMinPt_) || 
      (mydimu.M() < MuMuMinInvMass_) ||
      (mydimu.M() > MuMuMaxInvMass_)) {
    return false;
  }
  // compute the distance between mumu vtx and beam spot 
  computeLS (mumu_KV->position().x(),mumu_KV->position().y(),0.0,
	     beamSpot_.position().x(),beamSpot_.position().y(),0.0,
	     mumu_KV->error().cxx(),mumu_KV->error().cyy(),0.0,
	     mumu_KV->error().matrix()(0,1),0.0,0.0,
	     beamSpot_.covariance()(0,0),beamSpot_.covariance()(1,1),0.0,
	     beamSpot_.covariance()(0,1),0.0,0.0,
	     &MuMuLSBS,&MuMuLSBSErr);
  
  if (MuMuLSBS/MuMuLSBSErr < MuMuMinLxySigmaBs_)  return false;

  computeCosAlpha(mumu_KP->currentState().globalMomentum().x(),
		  mumu_KP->currentState().globalMomentum().y(), 
		  0.0,
		  mumu_KV->position().x() - beamSpot_.position().x(),
		  mumu_KV->position().y() - beamSpot_.position().y(),
		  0.0,
		  mumu_KP->currentState().kinematicParametersError().matrix()(3,3),
		  mumu_KP->currentState().kinematicParametersError().matrix()(4,4),
		  0.0,
		  mumu_KP->currentState().kinematicParametersError().matrix()(3,4),
		  0.0,
		  0.0,
		  mumu_KV->error().cxx() + beamSpot_.covariance()(0,0),
		  mumu_KV->error().cyy() + beamSpot_.covariance()(1,1),
		  0.0,
		  mumu_KV->error().matrix()(0,1) + beamSpot_.covariance()(0,1),
		  0.0,
		  0.0,
		  &MuMuCosAlphaBS,&MuMuCosAlphaBSErr);	  
  
  if (MuMuCosAlphaBS < MuMuMinCosAlphaBs_)  return false;
  
  return true; 
}


  
bool
BToKstarMuMu::matchMuonTrack (const edm::Event& iEvent, 
			      const reco::TrackRef theTrackRef)
{
  if ( theTrackRef.isNull() ) return false;

  edm::Handle< vector<pat::Muon> > thePATMuonHandle; 
  iEvent.getByLabel(MuonLabel_, thePATMuonHandle);
  
  reco::TrackRef muTrackRef; 
  for (vector<pat::Muon>::const_iterator iMuon = thePATMuonHandle->begin();
       iMuon != thePATMuonHandle->end(); iMuon++){

    muTrackRef = iMuon->innerTrack();
    if ( muTrackRef.isNull() ) continue; 

    if (muTrackRef == theTrackRef) return true; 
  }
  
  return false; 
}


bool
BToKstarMuMu::matchMuonTracks (const edm::Event& iEvent, 
			       const vector<reco::TrackRef> theTracks)
{
  reco::TrackRef theTrackRef; 
  for(unsigned int j = 0; j < theTracks.size(); ++j) {
    theTrackRef = theTracks[j];
    if ( matchMuonTrack(iEvent, theTrackRef) ) return true; 
  }
  return false; 
}


bool 
BToKstarMuMu::hasGoodKshortVertex(const vector<reco::TrackRef> theDaughterTracks, 
				  RefCountedKinematicTree &ksVertexFitTree)
{
  reco::TransientTrack pion1TT(theDaughterTracks[0], &(*bFieldHandle_) );
  reco::TransientTrack pion2TT(theDaughterTracks[1], &(*bFieldHandle_) );
  
  KinematicParticleFactoryFromTransientTrack pFactory;

  float chi = 0.;
  float ndf = 0.;
  vector<RefCountedKinematicParticle> pionParticles;
  pionParticles.push_back(pFactory.particle(pion1TT,PionMass_,chi,ndf,PionMassErr_));
  pionParticles.push_back(pFactory.particle(pion2TT,PionMass_,chi,ndf,PionMassErr_));

  KinematicParticleVertexFitter fitter;   
  ksVertexFitTree = fitter.fit(pionParticles); 
  if (!ksVertexFitTree->isValid()) return false;

  return true; 
}

bool 
BToKstarMuMu::hasGoodKshortVertexMKC(const vector<reco::TrackRef> theDaughterTracks, 
				     RefCountedKinematicTree &ksVertexFitTree)
{
  if ( !hasGoodKshortVertex(theDaughterTracks, ksVertexFitTree) ) return false; 
  KinematicParticleFitter csFitterKs; 
  KinematicConstraint * ks_c = new MassKinematicConstraint(KshortMass_, KshortMassErr_);
  ksVertexFitTree = csFitterKs.fit(ks_c,ksVertexFitTree);
  
  delete ks_c; 
  if (!ksVertexFitTree->isValid()) return false; 
  return true;
}


bool 
BToKstarMuMu::matchKshortTrack (const edm::Event& iEvent, 
				const reco::TrackRef theTrackRef)
{
  if ( theTrackRef.isNull() ) {
    // edm::LogInfo("myKshortMatch") << "Null track ref!"; 
    return false;
  }

  edm::Handle<reco::VertexCompositeCandidateCollection> theKshorts;
  iEvent.getByLabel(KshortLabel_, theKshorts);

  reco::TrackRef thePionTrackRef;

  for ( reco::VertexCompositeCandidateCollection::const_iterator iKshort 
	  = theKshorts->begin(); iKshort != theKshorts->end(); ++iKshort) {

    thePionTrackRef = (*(dynamic_cast<const reco::RecoChargedCandidate *>
			 (iKshort->daughter(0)))).track(); 

    if ( ! thePionTrackRef.isNull() && thePionTrackRef == theTrackRef) return true; 

    thePionTrackRef = (*(dynamic_cast<const reco::RecoChargedCandidate *>
			 (iKshort->daughter(1)))).track(); 
    
    if ( ! thePionTrackRef.isNull() && thePionTrackRef == theTrackRef) return true; 
  }
  
  return false; 
}



bool 
BToKstarMuMu::hasGoodKstarChargedMass(RefCountedKinematicTree vertexFitTree)
{
  vertexFitTree->movePointerToTheTop();
  TLorentzVector myks, mypi, mykstar; 

  vertexFitTree->movePointerToTheFirstChild(); // mu1
  vertexFitTree->movePointerToTheNextChild(); // mu2
  vertexFitTree->movePointerToTheNextChild(); // pion 

  RefCountedKinematicParticle pion_KP = vertexFitTree->currentParticle();

  mypi.SetXYZM(pion_KP->currentState().globalMomentum().x(), 
	       pion_KP->currentState().globalMomentum().y(), 
	       pion_KP->currentState().globalMomentum().z(), 
	       pion_KP->currentState().mass()); 

  vertexFitTree->movePointerToTheNextChild(); // ks 
  RefCountedKinematicParticle ks_KP = vertexFitTree->currentParticle();

  myks.SetXYZM(ks_KP->currentState().globalMomentum().x(), 
	       ks_KP->currentState().globalMomentum().y(), 
	       ks_KP->currentState().globalMomentum().z(), 
	       ks_KP->currentState().mass()); 
 
  mykstar = myks + mypi; 
  if ( mykstar.M() < KstarMinMass_  
       ||  mykstar.M() > KstarMaxMass_ )  return false; 

  return true; 
}


bool 
BToKstarMuMu::hasGoodBuMass(RefCountedKinematicTree vertexFitTree)
{
  vertexFitTree->movePointerToTheTop();
  RefCountedKinematicParticle b_KP = vertexFitTree->currentParticle();
  if ( b_KP->currentState().mass() > BMaxMass_ ) return false;  
  return true; 
}


bool 
BToKstarMuMu::hasGoodBuVertex(const reco::TrackRef mu1Track,
			      const reco::TrackRef mu2Track,
			      const vector<reco::TrackRef> kshortDaughterTracks, 
			      const reco::TrackRef pionTrack, 
			      RefCountedKinematicTree &vertexFitTree, 
			      RefCountedKinematicTree &ksVertexFitTree)
{
  if ( ! hasGoodKshortVertexMKC(kshortDaughterTracks, ksVertexFitTree) ) return false; 

  ksVertexFitTree->movePointerToTheTop();
  RefCountedKinematicParticle ks_KP = ksVertexFitTree->currentParticle();

  // do vertex fit for Bu
  KinematicParticleFactoryFromTransientTrack pFactory;
  reco::TransientTrack mu1TT(mu1Track, &(*bFieldHandle_) );
  reco::TransientTrack mu2TT(mu2Track, &(*bFieldHandle_) );
  reco::TransientTrack pionTT(pionTrack, &(*bFieldHandle_) );

  float chi = 0.;
  float ndf = 0.;
  vector<RefCountedKinematicParticle> vFitMCParticles;
  vFitMCParticles.push_back(pFactory.particle(mu1TT,MuonMass_,chi,ndf,MuonMassErr_));
  vFitMCParticles.push_back(pFactory.particle(mu2TT,MuonMass_,chi,ndf,MuonMassErr_));
  vFitMCParticles.push_back(pFactory.particle(pionTT, PionMass_, chi, ndf, PionMassErr_));
  vFitMCParticles.push_back(ks_KP);

  KinematicParticleVertexFitter fitter;   
  vertexFitTree = fitter.fit(vFitMCParticles);
  if (!vertexFitTree->isValid()) return false; 

  vertexFitTree->movePointerToTheTop();
  RefCountedKinematicVertex bDecayVertexMC = vertexFitTree->currentDecayVertex();
  if ( !bDecayVertexMC->vertexIsValid()) return false; 
 
  if ( bDecayVertexMC->chiSquared()<0
       || bDecayVertexMC->chiSquared()>1000 ) return false; 

  return true; 
}

bool
BToKstarMuMu::hasGoodDimuonKshortMass(const reco::TrackRef muTrackm, 
				      const reco::TrackRef muTrackp, 
				      const reco::VertexCompositeCandidate iKshort, 
				      double &dimu_ks_mass)
{
  // book the mumu mass histogram 
  TLorentzVector mu1, mu2, dimu; 
  mu1.SetXYZM(muTrackm->px(), muTrackm->py(), muTrackm->pz(), MuonMass_); 
  mu2.SetXYZM(muTrackp->px(), muTrackp->py(), muTrackp->pz(), MuonMass_); 
  dimu = mu1 + mu2; 
  BToKstarMuMuFigures[h_mumumass]->Fill(dimu.M()); 
  
  BToKstarMuMuFigures[h_kshortmass]->Fill(iKshort.mass()); 
  
  TLorentzVector kshort, dimukshort; 
  kshort.SetXYZM(iKshort.px(), iKshort.py(), iKshort.pz(), KshortMass_); 
  dimukshort = dimu + kshort; 

  dimu_ks_mass = dimukshort.M(); 

  if ( dimukshort.M() > DimukshortMinMass_ 
       && dimukshort.M() < DimukshortMaxMass_ ) return true; 
  return false; 
}


void 
BToKstarMuMu::saveBuToKstarMuMu(RefCountedKinematicTree vertexFitTree){

  vertexFitTree->movePointerToTheTop(); // B+ or B- 
  RefCountedKinematicParticle b_KP = vertexFitTree->currentParticle();

  bpx->push_back(b_KP->currentState().globalMomentum().x());
  bpxerr->push_back( sqrt( b_KP->currentState().kinematicParametersError().matrix()(3,3) ) );
  bpy->push_back(b_KP->currentState().globalMomentum().y());
  bpyerr->push_back( sqrt( b_KP->currentState().kinematicParametersError().matrix()(4,4) ) );
  bpz->push_back(b_KP->currentState().globalMomentum().z());
  bpzerr->push_back( sqrt( b_KP->currentState().kinematicParametersError().matrix()(5,5) ) );
  bmass->push_back(b_KP->currentState().mass()); 
  bmasserr->push_back( sqrt( b_KP->currentState().kinematicParametersError().matrix()(6,6) ) );

  vertexFitTree->movePointerToTheFirstChild(); // mu1 
  RefCountedKinematicParticle mu1_KP = vertexFitTree->currentParticle();
  vertexFitTree->movePointerToTheNextChild();  // mu2 
  RefCountedKinematicParticle mu2_KP = vertexFitTree->currentParticle();

  RefCountedKinematicParticle mup_KP, mum_KP ;

  if ( mu1_KP->currentState().particleCharge() > 0 ) mup_KP = mu1_KP;
  if ( mu1_KP->currentState().particleCharge() < 0 ) mum_KP = mu1_KP;
  if ( mu2_KP->currentState().particleCharge() > 0 ) mup_KP = mu2_KP;
  if ( mu2_KP->currentState().particleCharge() < 0 ) mum_KP = mu2_KP;

  
  muppx->push_back(mup_KP->currentState().globalMomentum().x());
  muppy->push_back(mup_KP->currentState().globalMomentum().y());
  muppz->push_back(mup_KP->currentState().globalMomentum().z());
  mupmass->push_back(mup_KP->currentState().mass());
  mupmasserr->push_back(mup_KP->currentState().kinematicParametersError().matrix()(6,6));
  
  mumpx->push_back(mum_KP->currentState().globalMomentum().x());
  mumpy->push_back(mum_KP->currentState().globalMomentum().y());
  mumpz->push_back(mum_KP->currentState().globalMomentum().z());
  mummass->push_back(mum_KP->currentState().mass());
  mummasserr->push_back(mum_KP->currentState().kinematicParametersError().matrix()(6,6));


  vertexFitTree->movePointerToTheNextChild();  // pion track 
  RefCountedKinematicParticle pi1_KP = vertexFitTree->currentParticle();
  trkchg->push_back(pi1_KP->currentState().particleCharge()); 
  trkpx->push_back(pi1_KP->currentState().globalMomentum().x());
  trkpy->push_back(pi1_KP->currentState().globalMomentum().y());
  trkpz->push_back(pi1_KP->currentState().globalMomentum().z());
  trkmass->push_back(pi1_KP->currentState().mass());
  trkmasserr->push_back(pi1_KP->currentState().kinematicParametersError().matrix()(6,6));

  vertexFitTree->movePointerToTheNextChild();  // ks 
  RefCountedKinematicParticle ks_KP = vertexFitTree->currentParticle();
  kspx->push_back(ks_KP->currentState().globalMomentum().x());
  kspy->push_back(ks_KP->currentState().globalMomentum().y());
  kspz->push_back(ks_KP->currentState().globalMomentum().z());
  ksmass->push_back(ks_KP->currentState().mass());
  ksmasserr->push_back(ks_KP->currentState().kinematicParametersError().matrix()(6,6));

}


void 
BToKstarMuMu::saveBuVertex(RefCountedKinematicTree vertexFitTree){
  vertexFitTree->movePointerToTheTop();
  RefCountedKinematicVertex b_KV = vertexFitTree->currentDecayVertex();
  bvtxcl->push_back( ChiSquaredProbability((double)(b_KV->chiSquared()),
					   (double)(b_KV->degreesOfFreedom())) );

  bvtxx->push_back((*b_KV).position().x());  
  bvtxxerr->push_back(sqrt( abs(b_KV->error().cxx()) ));
  bvtxy->push_back((*b_KV).position().y());  
  bvtxyerr->push_back(sqrt( abs(b_KV->error().cyy()) ));
  bvtxz->push_back((*b_KV).position().z()); 
  bvtxzerr->push_back(sqrt( abs(b_KV->error().czz()) ));
  
}

void 
BToKstarMuMu::saveBuCosAlpha(RefCountedKinematicTree vertexFitTree)
{
  // alpha is the angle in the transverse plane between the B0 momentum
  // and the seperation between the B0 vertex and the beamspot

  vertexFitTree->movePointerToTheTop();
  RefCountedKinematicParticle b_KP = vertexFitTree->currentParticle();
  RefCountedKinematicVertex b_KV = vertexFitTree->currentDecayVertex();
  
  
  double cosAlphaBS, cosAlphaBSErr;
  computeCosAlpha(b_KP->currentState().globalMomentum().x(),
		  b_KP->currentState().globalMomentum().y(), 
		  b_KP->currentState().globalMomentum().z(),
		  b_KV->position().x() - beamSpot_.position().x(),
		  b_KV->position().y() - beamSpot_.position().y(),
		  b_KV->position().z() - beamSpot_.position().z(),
		  b_KP->currentState().kinematicParametersError().matrix()(3,3),
		  b_KP->currentState().kinematicParametersError().matrix()(4,4),
		  b_KP->currentState().kinematicParametersError().matrix()(5,5),
		  b_KP->currentState().kinematicParametersError().matrix()(3,4),
		  b_KP->currentState().kinematicParametersError().matrix()(3,5),
		  b_KP->currentState().kinematicParametersError().matrix()(4,5),
		  b_KV->error().cxx() + beamSpot_.covariance()(0,0),
		  b_KV->error().cyy() + beamSpot_.covariance()(1,1),
		  b_KV->error().czz() + beamSpot_.covariance()(2,2),
		  b_KV->error().matrix()(0,1) + beamSpot_.covariance()(0,1),
		  b_KV->error().matrix()(0,2) + beamSpot_.covariance()(0,2),
		  b_KV->error().matrix()(1,2) + beamSpot_.covariance()(1,2),
		  &cosAlphaBS,&cosAlphaBSErr);	  

  bcosalphabs->push_back(cosAlphaBS);
  bcosalphabserr->push_back(cosAlphaBSErr); 
  
}



bool 
BToKstarMuMu::matchPrimaryVertexTracks () 
{
  vector<reco::TransientTrack> vertexTracks;
  for (vector<reco::TrackBaseRef>::const_iterator iTrack = primaryVertex_.tracks_begin(); 
       iTrack != primaryVertex_.tracks_end(); iTrack++){
    reco::TrackRef trackRef = iTrack->castTo<reco::TrackRef>();
  }
  return false; 
}

void 
BToKstarMuMu::saveBuLsig(RefCountedKinematicTree vertexFitTree)
{
  vertexFitTree->movePointerToTheTop();
  RefCountedKinematicVertex b_KV = vertexFitTree->currentDecayVertex();
  double LSBS, LSBSErr; 

  computeLS (b_KV->position().x(), b_KV->position().y(), 0.0,
	     beamSpot_.position().x(), beamSpot_.position().y(), 0.0,
	     b_KV->error().cxx(), b_KV->error().cyy(), 0.0,
	     b_KV->error().matrix()(0,1), 0.0, 0.0, 
	     beamSpot_.covariance()(0,0), beamSpot_.covariance()(1,1), 0.0,
	     beamSpot_.covariance()(0,1), 0.0, 0.0,
	     &LSBS,&LSBSErr);
 
  blsbs->push_back(LSBS);
  blsbserr->push_back(LSBSErr); 
 
}

void 
BToKstarMuMu::computeCtau(RefCountedKinematicTree vertexFitTree, 
			  double &bctau, double &bctauerr)
{
  //calculate ctau = (mB*(Bvtx-Pvtx)*pB)/(|pB|**2)

  vertexFitTree->movePointerToTheTop();
  RefCountedKinematicParticle b_KP = vertexFitTree->currentParticle();
  RefCountedKinematicVertex b_KV = vertexFitTree->currentDecayVertex();

  double betagamma = (b_KP->currentState().globalMomentum().mag()/BuMass_);

  // calculate ctau error. Momentum error is negligible compared to
  // the vertex errors, so don't worry about it

  GlobalPoint BVP = GlobalPoint( b_KV->position() );
  GlobalPoint PVP = GlobalPoint( primaryVertex_.position().x(),
				 primaryVertex_.position().y(), 
				 primaryVertex_.position().z() );
  GlobalVector sep3D = BVP-PVP;
  GlobalVector pBV = b_KP->currentState().globalMomentum();	      
  bctau = (BuMass_* (sep3D.dot(pBV)))/(pBV.dot(pBV));
  
  GlobalError BVE = b_KV->error();
  GlobalError PVE = GlobalError( primaryVertex_.error() );
  VertexDistance3D theVertexDistance3D; 
  Measurement1D TheMeasurement = theVertexDistance3D.distance( VertexState(BVP, BVE), 
							       VertexState(PVP, PVE) );
  double myError = TheMeasurement.error();	 
  
  //  ctau is defined by the portion of the flight distance along
  //  the compoenent of the B momementum, so only consider the error
  //  of that component, too, which is accomplished by scaling by
  //  ((VB-VP)(dot)PB)/|VB-VP|*|PB|	 
  
  double scale = abs( (sep3D.dot(pBV))/(sep3D.mag()*pBV.mag()) );    	       
  bctauerr =  (myError*scale)/betagamma;

}

double 
BToKstarMuMu::computeEta (double Px,
			  double Py,
			  double Pz)
{
  double P = sqrt(Px*Px + Py*Py + Pz*Pz);
  return 0.5*log((P + Pz) / (P - Pz));
}

double 
BToKstarMuMu::computePhi (double Px,
			  double Py,
			  double Pz)
{
  double phi = atan(Py / Px);
  if (Px < 0 && Py < 0) phi = phi - PI;
  if (Px < 0 && Py > 0) phi = phi + PI;
  return phi;
}

double
BToKstarMuMu::computeEtaPhiDistance (double Px1,
				     double Py1,
				     double Pz1,
				     double Px2,
				     double Py2,
				     double Pz2)
{
  double phi1 = computePhi (Px1,Py1,Pz1);
  double eta1 = computeEta (Px1,Py1,Pz1);
  double phi2 = computePhi (Px2,Py2,Pz2);
  double eta2 = computeEta (Px2,Py2,Pz2);
  return sqrt((eta1-eta2) * (eta1-eta2) + (phi1-phi2) * (phi1-phi2));
}


void 
BToKstarMuMu::saveBuCtau(RefCountedKinematicTree vertexFitTree)
{
  double bctau_temp, bctauerr_temp; 
  computeCtau(vertexFitTree, bctau_temp, bctauerr_temp); 
  bctau->push_back(bctau_temp);
  bctauerr->push_back(bctauerr_temp); 
}

void
BToKstarMuMu::saveGenInfo(const edm::Event& iEvent){
  edm::Handle<reco::GenParticleCollection> genparticles;
  iEvent.getByLabel(GenParticlesLabel_, genparticles );

  // loop over all gen particles
  for( size_t i = 0; i < genparticles->size(); ++i ) {
    const reco::GenParticle & b = (*genparticles)[i];

    // only select B+ or B- candidate 
    if ( abs(b.pdgId()) != BPLUS_PDG_ID ) continue; 
    
    int imum(-1), imup(-1), ikst(-1), iks(-1), ipi(-1), ipip(-1), ipim(-1); 

    // loop over all B+/- daughters 
    for ( size_t j = 0; j < b.numberOfDaughters(); ++j){
      const reco::Candidate  &dau = *(b.daughter(j));
      if (dau.pdgId() == MUONMINUS_PDG_ID) imum = j; 
      if (dau.pdgId() == -MUONMINUS_PDG_ID) imup = j; 
      if (abs(dau.pdgId()) == KSTARPLUS_PDG_ID) ikst = j;  
    }

    if (ikst == -1 || imum == -1 || imup == -1) continue; 

    const reco::Candidate & kst = *(b.daughter(ikst));
    
    for ( size_t j = 0; j < kst.numberOfDaughters(); ++j){
      const reco::Candidate  &dau = *(kst.daughter(j));
      if (abs(dau.pdgId()) == KSHORTZERO_PDG_ID) iks = j; 
      if (abs(dau.pdgId()) == PIONPLUS_PDG_ID) ipi = j; 
    }

    if (iks == -1 || ipi == -1) continue; 
   
    const reco::Candidate & ks = *(kst.daughter(iks));
    for ( size_t j = 0; j < ks.numberOfDaughters(); ++j){
      const reco::Candidate  &dau = *(ks.daughter(j));
      if ( dau.pdgId() == PIONPLUS_PDG_ID) ipip = j; 
      if ( dau.pdgId() == -PIONPLUS_PDG_ID) ipim = j; 
    }
   
    if (ipip == -1 || ipim == -1) continue; 
   
    const reco::Candidate & pip = *(ks.daughter(ipip));
    const reco::Candidate & pim = *(ks.daughter(ipim));
    const reco::Candidate & pi = *(kst.daughter(ipi));
    const reco::Candidate & mum = *(b.daughter(imum));
    const reco::Candidate & mup = *(b.daughter(imup));

    genbchg = b.charge(); 
    genbpx = b.px();
    genbpy = b.py();
    genbpz = b.pz();
 
    genmumpx = mum.px();
    genmumpy = mum.py();
    genmumpz = mum.pz();

    genmuppx = mup.px();
    genmuppy = mup.py();
    genmuppz = mup.pz();

    genkstpx = kst.px();
    genkstpy = kst.py();
    genkstpz = kst.pz();
    
    genkspx = ks.px();
    genkspy = ks.py();
    genkspz = ks.pz();

    genksvtxx = ks.vx();
    genksvtxy = ks.vy();
    genksvtxz = ks.vz();

    gentrkchg = pi.charge(); 
    gentrkpx = pi.px();
    gentrkpy = pi.py();
    gentrkpz = pi.pz();

    // save kshort pions 
    genpippx = pip.px();
    genpippy = pip.py();
    genpippz = pip.pz();
    genpimpx = pim.px();
    genpimpy = pim.py();
    genpimpz = pim.pz();
  }
}

bool
BToKstarMuMu::isGenKstarCharged(const reco::Candidate *p){
  if ( abs(p->pdgId()) != KSTARPLUS_PDG_ID ) 
    return false; 
  
  int nGenKshort = 0; 
  for ( size_t j = 0; j < p->numberOfDaughters(); ++j){
    const reco::Candidate * dau = p->daughter(j);
    
    if ( isGenKshort(dau) )  nGenKshort ++;  
  }
  return true; 
}


bool
BToKstarMuMu::isGenKshort(const reco::Candidate *p){
  if ( abs(p->pdgId()) != KSHORTZERO_PDG_ID )  return false; 
  return true; 
}


bool
BToKstarMuMu::isGenMuonP(const reco::Candidate *p){
  if ( abs(p->pdgId()) != MUONMINUS_PDG_ID )
    return false; 
  
  return true; 
}


void 
BToKstarMuMu::saveKshortVariables(RefCountedKinematicTree ksVertexFitTree)
{
  ksVertexFitTree->movePointerToTheTop();
  RefCountedKinematicVertex ks_vertex = ksVertexFitTree->currentDecayVertex();

  ksvtxcl->push_back( ChiSquaredProbability((double)(ks_vertex->chiSquared()),
					    (double)(ks_vertex->degreesOfFreedom())) );

  ksvtxx->push_back(ks_vertex->position().x());
  ksvtxy->push_back(ks_vertex->position().y());
  ksvtxz->push_back(ks_vertex->position().z());
  
  // Kshort vertex distance to the beam spot 

  double KshortLSBS, KshortLSBSErr; 

  computeLS (ks_vertex->position().x(),ks_vertex->position().y(),0.0,
	     beamSpot_.position().x(),beamSpot_.position().y(),0.0,
	     ks_vertex->error().cxx(),ks_vertex->error().cyy(),0.0,
	     ks_vertex->error().matrix()(0,1),0.0,0.0,
	     beamSpot_.covariance()(0,0),beamSpot_.covariance()(1,1),0.0,
	     beamSpot_.covariance()(0,1),0.0,0.0,
	     &KshortLSBS,&KshortLSBSErr);
  
  kslsbs->push_back(KshortLSBS); 
  kslsbserr->push_back(KshortLSBSErr); 


  ksVertexFitTree->movePointerToTheFirstChild(); // Kshort pion1
  RefCountedKinematicParticle ksPi1 = ksVertexFitTree->currentParticle();

  ksVertexFitTree->movePointerToTheNextChild(); // Kshort pion2 
  RefCountedKinematicParticle ksPi2 = ksVertexFitTree->currentParticle();


  KinematicParameters ksPi1KP = ksPi1->currentState().kinematicParameters();
  KinematicParameters ksPi2KP = ksPi2->currentState().kinematicParameters();
  KinematicParameters ksPipKP;
  KinematicParameters ksPimKP;

  if ( ksPi1->currentState().particleCharge() > 0 ) ksPipKP = ksPi1KP;
  if ( ksPi1->currentState().particleCharge() < 0 ) ksPimKP = ksPi1KP;
  if ( ksPi2->currentState().particleCharge() > 0 ) ksPipKP = ksPi2KP;
  if ( ksPi2->currentState().particleCharge() < 0 ) ksPimKP = ksPi2KP;

  pippx->push_back(ksPipKP.momentum().x());
  pippy->push_back(ksPipKP.momentum().y());
  pippz->push_back(ksPipKP.momentum().z());

  pimpx->push_back(ksPimKP.momentum().x());
  pimpy->push_back(ksPimKP.momentum().y());
  pimpz->push_back(ksPimKP.momentum().z());
}


void 
BToKstarMuMu::saveSoftMuonVariables(pat::Muon iMuonM, pat::Muon iMuonP, 
				    reco::TrackRef muTrackm, reco::TrackRef muTrackp)
{
  mumisgoodmuon->push_back(muon::isGoodMuon(iMuonM, muon::TMOneStationTight)); 
  mupisgoodmuon->push_back(muon::isGoodMuon(iMuonP, muon::TMOneStationTight)); 
  mumnpixhits->push_back(muTrackm->hitPattern().numberOfValidPixelHits()); 
  mupnpixhits->push_back(muTrackp->hitPattern().numberOfValidPixelHits()); 
  mumnpixlayers->push_back(muTrackm->hitPattern().pixelLayersWithMeasurement()); 
  mupnpixlayers->push_back(muTrackp->hitPattern().pixelLayersWithMeasurement()); 

  mumntrkhits->push_back(muTrackm->hitPattern().numberOfValidTrackerHits()); 
  mupntrkhits->push_back(muTrackp->hitPattern().numberOfValidTrackerHits()); 
  mumntrklayers->push_back(muTrackm->hitPattern().trackerLayersWithMeasurement()); 
  mupntrklayers->push_back(muTrackp->hitPattern().trackerLayersWithMeasurement()); 

  mumnormchi2->push_back(muTrackm->normalizedChi2()); 
  mupnormchi2->push_back(muTrackp->normalizedChi2()); 

  mumdxyvtx->push_back(muTrackm->dxy(primaryVertex_.position())); 
  mupdxyvtx->push_back(muTrackp->dxy(primaryVertex_.position())); 
	  
  mumdzvtx->push_back(muTrackm->dz(primaryVertex_.position())); 
  mupdzvtx->push_back(muTrackp->dz(primaryVertex_.position())); 

  mumpt->push_back(muTrackm->pt()); 
  muppt->push_back(muTrackp->pt()); 
  
  mumeta->push_back(muTrackm->eta()); 
  mupeta->push_back(muTrackp->eta()); 
  
}

void 
BToKstarMuMu::saveDimuVariables(double DCAmumBS, double DCAmumBSErr, 
				double DCAmupBS, double DCAmupBSErr,
				double DCAmumu,  double mu_mu_vtx_cl, 
				double MuMuLSBS, double MuMuLSBSErr, 
				double MuMuCosAlphaBS, double MuMuCosAlphaBSErr, 
				double dimu_ks_mass)
{
  mumdcabs->push_back(DCAmumBS);
  mumdcabserr->push_back(DCAmumBSErr);

  mupdcabs->push_back(DCAmupBS);
  mupdcabserr->push_back(DCAmupBSErr);

  mumudca->push_back(DCAmumu);
  mumuvtxcl->push_back(mu_mu_vtx_cl); 
  mumulsbs->push_back(MuMuLSBS);
  mumulsbserr->push_back(MuMuLSBSErr);
  mumucosalphabs->push_back(MuMuCosAlphaBS);
  mumucosalphabserr->push_back(MuMuCosAlphaBSErr); 

  dimuksmass->push_back(dimu_ks_mass); 
}


void 
BToKstarMuMu::saveMuonTriggerMatches(const pat::Muon iMuonM, const pat::Muon iMuonP)
{	
  string mum_matched_lastfilter_name = ""; 
  string mup_matched_lastfilter_name = ""; 

  for(vector<string>::iterator it = triggernames->begin(); 
      it != triggernames->end(); ++it) {
  
    string hltLastFilterName = mapTriggerToLastFilter_[*it] ; 

    const pat::TriggerObjectStandAloneCollection mumHLTMatches 
      = iMuonM.triggerObjectMatchesByFilter( hltLastFilterName );
    const pat::TriggerObjectStandAloneCollection mupHLTMatches 
      = iMuonP.triggerObjectMatchesByFilter( hltLastFilterName );
    
    if ( mumHLTMatches.size() > 0 ) 
      mum_matched_lastfilter_name.append(hltLastFilterName+" ") ; 
    
    if ( mupHLTMatches.size() > 0 ) 
      mup_matched_lastfilter_name.append(hltLastFilterName+" ") ; 
  }
  
  mumtriglastfilter->push_back(mum_matched_lastfilter_name); 
  muptriglastfilter->push_back(mup_matched_lastfilter_name); 
}

void
BToKstarMuMu::saveTruthMatch(const edm::Event& iEvent){
  double deltaEtaPhi; 

  for (vector<int>::size_type i = 0; i < bmass->size(); i++) {
    // truth match with mu-
    deltaEtaPhi = computeEtaPhiDistance(genmumpx, genmumpy, genmumpz, 
					mumpx->at(i), mumpy->at(i), mumpz->at(i)); 
    if (deltaEtaPhi < TruthMatchMuonMaxR_) {
      istruemum->push_back(true); 
    } else {
      istruemum->push_back(false); 
    }
    
    // truth match with mu+
    deltaEtaPhi = computeEtaPhiDistance(genmuppx, genmuppy, genmuppz, 
					muppx->at(i), muppy->at(i), muppz->at(i)); 
  
    if (deltaEtaPhi < TruthMatchMuonMaxR_) {
      istruemup->push_back(true); 
    }
    else {
      istruemup->push_back(false); 
    }

    // truth match with Ks pi+
    bool istruepip = false; 
    
    deltaEtaPhi = computeEtaPhiDistance(genpippx, genpippy, genpippz, 
					pippx->at(i), pippy->at(i), pippz->at(i)); 
    
    if (deltaEtaPhi < TruthMatchPionMaxR_) istruepip = true; 
    
    // truth match with Ks pi-

    bool istruepim = false; 
    deltaEtaPhi = computeEtaPhiDistance(genpimpx, genpimpy, genpimpz, 
					pimpx->at(i), pimpy->at(i), pimpz->at(i)); 
    
    if (deltaEtaPhi < TruthMatchPionMaxR_) istruepim = true; 
    
    // truth match Ks vertex 
    float deltaRksvtx = sqrt( (genksvtxx - ksvtxx->at(i))*
			      (genksvtxx - ksvtxx->at(i)) +
			      (genksvtxy - ksvtxy->at(i))*
			      (genksvtxy - ksvtxy->at(i)) +
			      (genksvtxz - ksvtxz->at(i))*
			      (genksvtxz - ksvtxz->at(i)) );	     

    if ( istruepip & istruepim && (deltaRksvtx<TruthMatchKsMaxVtx_) ){
      istrueks->push_back(true);
    } else {
      istrueks->push_back(false); 
    }
    
    // truth match with pion track
    deltaEtaPhi = computeEtaPhiDistance(gentrkpx, gentrkpy, gentrkpz, 
					trkpx->at(i), trkpy->at(i), trkpz->at(i)); 
    if (deltaEtaPhi < TruthMatchPionMaxR_){
      istruetrk->push_back(true);
    } else {
      istruetrk->push_back(false); 
    }
 
    // truth match with B+/B-
    if ( istruemum->back() && istruemup->back() 
	 && istrueks->back() && istruetrk->back()) {
      istruebu->push_back(true); 
    } else {
      istruebu->push_back(false); 
    }
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(BToKstarMuMu);
