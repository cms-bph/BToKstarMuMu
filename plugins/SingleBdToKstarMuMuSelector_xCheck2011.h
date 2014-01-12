//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Dec 28 19:32:01 2013 by ROOT version 5.32/03
// from TTree B0KstMuMuSingleCandNTuple/B0KstMuMuSingleCandNTuple
// found on file: singleCand_B0ToKstMuMu_DataRRPRv4v5v6v1B_NTuples_Merged.root
//////////////////////////////////////////////////////////

#ifndef SingleBdToKstarMuMuSelector_xCheck2011_h
#define SingleBdToKstarMuMuSelector_xCheck2011_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

using namespace std;

// Fixed size dimensions of array or collections stored in the TTree if any.

class SingleBdToKstarMuMuSelector_xCheck2011 : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   UInt_t          runN;
   UInt_t          eventN;
   UInt_t          recoVtxN;
   Double_t        evWeight;
   vector<string>  *TrigTable;
   vector<int>     *TrigPrescales;
   UInt_t          nB;
   Double_t        priVtxCL;
   Double_t        priVtxX;
   Double_t        priVtxY;
   Double_t        priVtxZ;
   Double_t        bsX;
   Double_t        bsY;
   vector<double>  *bMass;
   vector<double>  *bMassE;
   vector<double>  *bBarMass;
   vector<double>  *bBarMassE;
   vector<double>  *bPx;
   vector<double>  *bPy;
   vector<double>  *bPz;
   vector<double>  *bunchXingMC;
   vector<double>  *numInteractionsMC;
   vector<double>  *trueNumInteractionsMC;
   vector<double>  *bVtxCL;
   vector<double>  *bVtxX;
   vector<double>  *bVtxY;
   vector<double>  *bVtxZ;
   vector<double>  *bCosAlphaVtx;
   vector<double>  *bCosAlphaVtxE;
   vector<double>  *bCosAlphaBS;
   vector<double>  *bCosAlphaBSE;
   vector<double>  *bLVtx;
   vector<double>  *bLVtxE;
   vector<double>  *bLBS;
   vector<double>  *bLBSE;
   vector<double>  *bDCAVtx;
   vector<double>  *bDCAVtxE;
   vector<double>  *bDCABS;
   vector<double>  *bDCABSE;
   vector<double>  *bctauPVBS;
   vector<double>  *bctauPVBSE;
   vector<double>  *kstMass;
   vector<double>  *kstMassE;
   vector<double>  *kstBarMass;
   vector<double>  *kstBarMassE;
   vector<double>  *kstPx;
   vector<double>  *kstPy;
   vector<double>  *kstPz;
   vector<double>  *kstVtxCL;
   vector<double>  *kstVtxX;
   vector<double>  *kstVtxY;
   vector<double>  *kstVtxZ;
   vector<double>  *mumuMass;
   vector<double>  *mumuMassE;
   vector<double>  *mumuPx;
   vector<double>  *mumuPy;
   vector<double>  *mumuPz;
   vector<double>  *mumuVtxCL;
   vector<double>  *mumuVtxX;
   vector<double>  *mumuVtxY;
   vector<double>  *mumuVtxZ;
   vector<double>  *mumuCosAlphaBS;
   vector<double>  *mumuCosAlphaBSE;
   vector<double>  *mumuLBS;
   vector<double>  *mumuLBSE;
   vector<double>  *mumuDCA;
   vector<bool>    *mumHighPurity;
   vector<double>  *mumCL;
   vector<double>  *mumNormChi2;
   vector<double>  *mumPx;
   vector<double>  *mumPy;
   vector<double>  *mumPz;
   vector<double>  *mumDCAVtx;
   vector<double>  *mumDCAVtxE;
   vector<double>  *mumDCABS;
   vector<double>  *mumDCABSE;
   vector<double>  *mumKinkChi2;
   vector<double>  *mumFracHits;
   vector<double>  *mumdxyVtx;
   vector<double>  *mumdzVtx;
   vector<string>  *mumCat;
   vector<int>     *mumNPixHits;
   vector<int>     *mumNPixLayers;
   vector<int>     *mumNTrkHits;
   vector<int>     *mumNTrkLayers;
   vector<int>     *mumNMuonHits;
   vector<int>     *mumNMatchStation;
   vector<string>  *mumTrig;
   vector<bool>    *mum2LastTrigFilter;
   vector<bool>    *mupHighPurity;
   vector<double>  *mupCL;
   vector<double>  *mupNormChi2;
   vector<double>  *mupPx;
   vector<double>  *mupPy;
   vector<double>  *mupPz;
   vector<double>  *mupDCAVtx;
   vector<double>  *mupDCAVtxE;
   vector<double>  *mupDCABS;
   vector<double>  *mupDCABSE;
   vector<double>  *mupKinkChi2;
   vector<double>  *mupFracHits;
   vector<double>  *mupdxyVtx;
   vector<double>  *mupdzVtx;
   vector<string>  *mupCat;
   vector<int>     *mupNPixHits;
   vector<int>     *mupNPixLayers;
   vector<int>     *mupNTrkHits;
   vector<int>     *mupNTrkLayers;
   vector<int>     *mupNMuonHits;
   vector<int>     *mupNMatchStation;
   vector<string>  *mupTrig;
   vector<bool>    *mup2LastTrigFilter;
   vector<bool>    *kstTrkmHighPurity;
   vector<double>  *kstTrkmCL;
   vector<double>  *kstTrkmNormChi2;
   vector<double>  *kstTrkmPx;
   vector<double>  *kstTrkmPy;
   vector<double>  *kstTrkmPz;
   vector<double>  *kstTrkmDCAVtx;
   vector<double>  *kstTrkmDCAVtxE;
   vector<double>  *kstTrkmDCABS;
   vector<double>  *kstTrkmDCABSE;
   vector<double>  *kstTrkmFracHits;
   vector<double>  *kstTrkmdxyVtx;
   vector<double>  *kstTrkmdzVtx;
   vector<int>     *kstTrkmNPixHits;
   vector<int>     *kstTrkmNPixLayers;
   vector<int>     *kstTrkmNTrkHits;
   vector<int>     *kstTrkmNTrkLayers;
   vector<string>  *kstTrkmMuMatch;
   vector<bool>    *kstTrkpHighPurity;
   vector<double>  *kstTrkpCL;
   vector<double>  *kstTrkpNormChi2;
   vector<double>  *kstTrkpPx;
   vector<double>  *kstTrkpPy;
   vector<double>  *kstTrkpPz;
   vector<double>  *kstTrkpDCAVtx;
   vector<double>  *kstTrkpDCAVtxE;
   vector<double>  *kstTrkpDCABS;
   vector<double>  *kstTrkpDCABSE;
   vector<double>  *kstTrkpFracHits;
   vector<double>  *kstTrkpdxyVtx;
   vector<double>  *kstTrkpdzVtx;
   vector<int>     *kstTrkpNPixHits;
   vector<int>     *kstTrkpNPixLayers;
   vector<int>     *kstTrkpNTrkHits;
   vector<int>     *kstTrkpNTrkLayers;
   vector<string>  *kstTrkpMuMatch;
   Int_t           genSignal;
   Int_t           genMuMuBG;
   Int_t           genMuMuBGnTrksm;
   Int_t           genMuMuBGnTrksp;
   Bool_t          genPsiPrompt;
   Bool_t          genSignHasFSR;
   Bool_t          genSignKstHasFSR;
   Bool_t          genSignPsiHasFSR;
   Double_t        genPriVtxX;
   Double_t        genPriVtxY;
   Double_t        genPriVtxZ;
   Double_t        genB0Mass;
   Double_t        genB0Px;
   Double_t        genB0Py;
   Double_t        genB0Pz;
   Double_t        genB0VtxX;
   Double_t        genB0VtxY;
   Double_t        genB0VtxZ;
   Double_t        genKstMass;
   Double_t        genKstPx;
   Double_t        genKstPy;
   Double_t        genKstPz;
   Double_t        genKstVtxX;
   Double_t        genKstVtxY;
   Double_t        genKstVtxZ;
   Double_t        genPsiMass;
   Double_t        genPsiVtxX;
   Double_t        genPsiVtxY;
   Double_t        genPsiVtxZ;
   Double_t        genMumPx;
   Double_t        genMumPy;
   Double_t        genMumPz;
   Bool_t          trueMumTriggered;
   Bool_t          trueMumInAcceptance;
   Double_t        genMupPx;
   Double_t        genMupPy;
   Double_t        genMupPz;
   Bool_t          trueMupTriggered;
   Bool_t          trueMupInAcceptance;
   Double_t        genKstTrkmPx;
   Double_t        genKstTrkmPy;
   Double_t        genKstTrkmPz;
   Double_t        genKstTrkpPx;
   Double_t        genKstTrkpPy;
   Double_t        genKstTrkpPz;
   vector<bool>    *truthMatchSignal;
   vector<bool>    *truthMatchMum;
   vector<bool>    *truthMatchMup;
   vector<bool>    *truthMatchTrkm;
   vector<bool>    *truthMatchTrkp;
   Double_t        evWeightE2;
   Double_t        B0MassArb;
   Double_t        B0pT;
   Double_t        B0Eta;
   Double_t        B0Phi;
   Bool_t          B0notB0bar;
   Int_t           TrigCat;
   Double_t        CosThetaKArb;
   Double_t        CosThetaMuArb;
   Double_t        PhiKstMuMuPlaneArb;

   // List of branches
   TBranch        *b_runN;   //!
   TBranch        *b_eventN;   //!
   TBranch        *b_recoVtxN;   //!
   TBranch        *b_evWeight;   //!
   TBranch        *b_TrigTable;   //!
   TBranch        *b_TrigPrescales;   //!
   TBranch        *b_nB;   //!
   TBranch        *b_priVtxCL;   //!
   TBranch        *b_priVtxX;   //!
   TBranch        *b_priVtxY;   //!
   TBranch        *b_priVtxZ;   //!
   TBranch        *b_bsX;   //!
   TBranch        *b_bsY;   //!
   TBranch        *b_bMass;   //!
   TBranch        *b_bMassE;   //!
   TBranch        *b_bBarMass;   //!
   TBranch        *b_bBarMassE;   //!
   TBranch        *b_bPx;   //!
   TBranch        *b_bPy;   //!
   TBranch        *b_bPz;   //!
   TBranch        *b_bunchXingMC;   //!
   TBranch        *b_numInteractionsMC;   //!
   TBranch        *b_trueNumInteractionsMC;   //!
   TBranch        *b_bVtxCL;   //!
   TBranch        *b_bVtxX;   //!
   TBranch        *b_bVtxY;   //!
   TBranch        *b_bVtxZ;   //!
   TBranch        *b_bCosAlphaVtx;   //!
   TBranch        *b_bCosAlphaVtxE;   //!
   TBranch        *b_bCosAlphaBS;   //!
   TBranch        *b_bCosAlphaBSE;   //!
   TBranch        *b_bLVtx;   //!
   TBranch        *b_bLVtxE;   //!
   TBranch        *b_bLBS;   //!
   TBranch        *b_bLBSE;   //!
   TBranch        *b_bDCAVtx;   //!
   TBranch        *b_bDCAVtxE;   //!
   TBranch        *b_bDCABS;   //!
   TBranch        *b_bDCABSE;   //!
   TBranch        *b_bctauPVBS;   //!
   TBranch        *b_bctauPVBSE;   //!
   TBranch        *b_kstMass;   //!
   TBranch        *b_kstMassE;   //!
   TBranch        *b_kstBarMass;   //!
   TBranch        *b_kstBarMassE;   //!
   TBranch        *b_kstPx;   //!
   TBranch        *b_kstPy;   //!
   TBranch        *b_kstPz;   //!
   TBranch        *b_kstVtxCL;   //!
   TBranch        *b_kstVtxX;   //!
   TBranch        *b_kstVtxY;   //!
   TBranch        *b_kstVtxZ;   //!
   TBranch        *b_mumuMass;   //!
   TBranch        *b_mumuMassE;   //!
   TBranch        *b_mumuPx;   //!
   TBranch        *b_mumuPy;   //!
   TBranch        *b_mumuPz;   //!
   TBranch        *b_mumuVtxCL;   //!
   TBranch        *b_mumuVtxX;   //!
   TBranch        *b_mumuVtxY;   //!
   TBranch        *b_mumuVtxZ;   //!
   TBranch        *b_mumuCosAlphaBS;   //!
   TBranch        *b_mumuCosAlphaBSE;   //!
   TBranch        *b_mumuLBS;   //!
   TBranch        *b_mumuLBSE;   //!
   TBranch        *b_mumuDCA;   //!
   TBranch        *b_mumHighPurity;   //!
   TBranch        *b_mumCL;   //!
   TBranch        *b_mumNormChi2;   //!
   TBranch        *b_mumPx;   //!
   TBranch        *b_mumPy;   //!
   TBranch        *b_mumPz;   //!
   TBranch        *b_mumDCAVtx;   //!
   TBranch        *b_mumDCAVtxE;   //!
   TBranch        *b_mumDCABS;   //!
   TBranch        *b_mumDCABSE;   //!
   TBranch        *b_mumKinkChi2;   //!
   TBranch        *b_mumFracHits;   //!
   TBranch        *b_mumdxyVtx;   //!
   TBranch        *b_mumdzVtx;   //!
   TBranch        *b_mumCat;   //!
   TBranch        *b_mumNPixHits;   //!
   TBranch        *b_mumNPixLayers;   //!
   TBranch        *b_mumNTrkHits;   //!
   TBranch        *b_mumNTrkLayers;   //!
   TBranch        *b_mumNMuonHits;   //!
   TBranch        *b_mumNMatchStation;   //!
   TBranch        *b_mumTrig;   //!
   TBranch        *b_mum2LastTrigFilter;   //!
   TBranch        *b_mupHighPurity;   //!
   TBranch        *b_mupCL;   //!
   TBranch        *b_mupNormChi2;   //!
   TBranch        *b_mupPx;   //!
   TBranch        *b_mupPy;   //!
   TBranch        *b_mupPz;   //!
   TBranch        *b_mupDCAVtx;   //!
   TBranch        *b_mupDCAVtxE;   //!
   TBranch        *b_mupDCABS;   //!
   TBranch        *b_mupDCABSE;   //!
   TBranch        *b_mupKinkChi2;   //!
   TBranch        *b_mupFracHits;   //!
   TBranch        *b_mupdxyVtx;   //!
   TBranch        *b_mupdzVtx;   //!
   TBranch        *b_mupCat;   //!
   TBranch        *b_mupNPixHits;   //!
   TBranch        *b_mupNPixLayers;   //!
   TBranch        *b_mupNTrkHits;   //!
   TBranch        *b_mupNTrkLayers;   //!
   TBranch        *b_mupNMuonHits;   //!
   TBranch        *b_mupNMatchStation;   //!
   TBranch        *b_mupTrig;   //!
   TBranch        *b_mup2LastTrigFilter;   //!
   TBranch        *b_kstTrkmHighPurity;   //!
   TBranch        *b_kstTrkmCL;   //!
   TBranch        *b_kstTrkmNormChi2;   //!
   TBranch        *b_kstTrkmPx;   //!
   TBranch        *b_kstTrkmPy;   //!
   TBranch        *b_kstTrkmPz;   //!
   TBranch        *b_kstTrkmDCAVtx;   //!
   TBranch        *b_kstTrkmDCAVtxE;   //!
   TBranch        *b_kstTrkmDCABS;   //!
   TBranch        *b_kstTrkmDCABSE;   //!
   TBranch        *b_kstTrkmFracHits;   //!
   TBranch        *b_kstTrkmdxyVtx;   //!
   TBranch        *b_kstTrkmdzVtx;   //!
   TBranch        *b_kstTrkmNPixHits;   //!
   TBranch        *b_kstTrkmNPixLayers;   //!
   TBranch        *b_kstTrkmNTrkHits;   //!
   TBranch        *b_kstTrkmNTrkLayers;   //!
   TBranch        *b_kstTrkmMuMatch;   //!
   TBranch        *b_kstTrkpHighPurity;   //!
   TBranch        *b_kstTrkpCL;   //!
   TBranch        *b_kstTrkpNormChi2;   //!
   TBranch        *b_kstTrkpPx;   //!
   TBranch        *b_kstTrkpPy;   //!
   TBranch        *b_kstTrkpPz;   //!
   TBranch        *b_kstTrkpDCAVtx;   //!
   TBranch        *b_kstTrkpDCAVtxE;   //!
   TBranch        *b_kstTrkpDCABS;   //!
   TBranch        *b_kstTrkpDCABSE;   //!
   TBranch        *b_kstTrkpFracHits;   //!
   TBranch        *b_kstTrkpdxyVtx;   //!
   TBranch        *b_kstTrkpdzVtx;   //!
   TBranch        *b_kstTrkpNPixHits;   //!
   TBranch        *b_kstTrkpNPixLayers;   //!
   TBranch        *b_kstTrkpNTrkHits;   //!
   TBranch        *b_kstTrkpNTrkLayers;   //!
   TBranch        *b_kstTrkpMuMatch;   //!
   TBranch        *b_genSignal;   //!
   TBranch        *b_genMuMuBG;   //!
   TBranch        *b_genMuMuBGnTrksm;   //!
   TBranch        *b_genMuMuBGnTrksp;   //!
   TBranch        *b_genPsiPrompt;   //!
   TBranch        *b_genSignHasFSR;   //!
   TBranch        *b_genSignKstHasFSR;   //!
   TBranch        *b_genSignPsiHasFSR;   //!
   TBranch        *b_genPriVtxX;   //!
   TBranch        *b_genPriVtxY;   //!
   TBranch        *b_genPriVtxZ;   //!
   TBranch        *b_genB0Mass;   //!
   TBranch        *b_genB0Px;   //!
   TBranch        *b_genB0Py;   //!
   TBranch        *b_genB0Pz;   //!
   TBranch        *b_genB0VtxX;   //!
   TBranch        *b_genB0VtxY;   //!
   TBranch        *b_genB0VtxZ;   //!
   TBranch        *b_genKstMass;   //!
   TBranch        *b_genKstPx;   //!
   TBranch        *b_genKstPy;   //!
   TBranch        *b_genKstPz;   //!
   TBranch        *b_genKstVtxX;   //!
   TBranch        *b_genKstVtxY;   //!
   TBranch        *b_genKstVtxZ;   //!
   TBranch        *b_genPsiMass;   //!
   TBranch        *b_genPsiVtxX;   //!
   TBranch        *b_genPsiVtxY;   //!
   TBranch        *b_genPsiVtxZ;   //!
   TBranch        *b_genMumPx;   //!
   TBranch        *b_genMumPy;   //!
   TBranch        *b_genMumPz;   //!
   TBranch        *b_trueMumTriggered;   //!
   TBranch        *b_trueMumInAcceptance;   //!
   TBranch        *b_genMupPx;   //!
   TBranch        *b_genMupPy;   //!
   TBranch        *b_genMupPz;   //!
   TBranch        *b_trueMupTriggered;   //!
   TBranch        *b_trueMupInAcceptance;   //!
   TBranch        *b_genKstTrkmPx;   //!
   TBranch        *b_genKstTrkmPy;   //!
   TBranch        *b_genKstTrkmPz;   //!
   TBranch        *b_genKstTrkpPx;   //!
   TBranch        *b_genKstTrkpPy;   //!
   TBranch        *b_genKstTrkpPz;   //!
   TBranch        *b_truthMatchSignal;   //!
   TBranch        *b_truthMatchMum;   //!
   TBranch        *b_truthMatchMup;   //!
   TBranch        *b_truthMatchTrkm;   //!
   TBranch        *b_truthMatchTrkp;   //!
   TBranch        *b_evWeightE2;   //!
   TBranch        *b_B0MassArb;   //!
   TBranch        *b_B0pT;   //!
   TBranch        *b_B0Eta;   //!
   TBranch        *b_B0Phi;   //!
   TBranch        *b_B0notB0bar;   //!
   TBranch        *b_TrigCat;   //!
   TBranch        *b_CosThetaKArb;   //!
   TBranch        *b_CosThetaMuArb;   //!
   TBranch        *b_PhiKstMuMuPlaneArb;   //!

   SingleBdToKstarMuMuSelector_xCheck2011(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~SingleBdToKstarMuMuSelector_xCheck2011() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   // ClassDef(SingleBdToKstarMuMuSelector_xCheck2011,0);
   int SelectB(string);
   bool HasGoodDimuon(int); 
   void SaveEvent(int); 
   void SaveGen();

};

#endif

#ifdef SingleBdToKstarMuMuSelector_xCheck2011_cxx
void SingleBdToKstarMuMuSelector_xCheck2011::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   TrigTable = 0;
   TrigPrescales = 0;
   bMass = 0;
   bMassE = 0;
   bBarMass = 0;
   bBarMassE = 0;
   bPx = 0;
   bPy = 0;
   bPz = 0;
   bunchXingMC = 0;
   numInteractionsMC = 0;
   trueNumInteractionsMC = 0;
   bVtxCL = 0;
   bVtxX = 0;
   bVtxY = 0;
   bVtxZ = 0;
   bCosAlphaVtx = 0;
   bCosAlphaVtxE = 0;
   bCosAlphaBS = 0;
   bCosAlphaBSE = 0;
   bLVtx = 0;
   bLVtxE = 0;
   bLBS = 0;
   bLBSE = 0;
   bDCAVtx = 0;
   bDCAVtxE = 0;
   bDCABS = 0;
   bDCABSE = 0;
   bctauPVBS = 0;
   bctauPVBSE = 0;
   kstMass = 0;
   kstMassE = 0;
   kstBarMass = 0;
   kstBarMassE = 0;
   kstPx = 0;
   kstPy = 0;
   kstPz = 0;
   kstVtxCL = 0;
   kstVtxX = 0;
   kstVtxY = 0;
   kstVtxZ = 0;
   mumuMass = 0;
   mumuMassE = 0;
   mumuPx = 0;
   mumuPy = 0;
   mumuPz = 0;
   mumuVtxCL = 0;
   mumuVtxX = 0;
   mumuVtxY = 0;
   mumuVtxZ = 0;
   mumuCosAlphaBS = 0;
   mumuCosAlphaBSE = 0;
   mumuLBS = 0;
   mumuLBSE = 0;
   mumuDCA = 0;
   mumHighPurity = 0;
   mumCL = 0;
   mumNormChi2 = 0;
   mumPx = 0;
   mumPy = 0;
   mumPz = 0;
   mumDCAVtx = 0;
   mumDCAVtxE = 0;
   mumDCABS = 0;
   mumDCABSE = 0;
   mumKinkChi2 = 0;
   mumFracHits = 0;
   mumdxyVtx = 0;
   mumdzVtx = 0;
   mumCat = 0;
   mumNPixHits = 0;
   mumNPixLayers = 0;
   mumNTrkHits = 0;
   mumNTrkLayers = 0;
   mumNMuonHits = 0;
   mumNMatchStation = 0;
   mumTrig = 0;
   mum2LastTrigFilter = 0;
   mupHighPurity = 0;
   mupCL = 0;
   mupNormChi2 = 0;
   mupPx = 0;
   mupPy = 0;
   mupPz = 0;
   mupDCAVtx = 0;
   mupDCAVtxE = 0;
   mupDCABS = 0;
   mupDCABSE = 0;
   mupKinkChi2 = 0;
   mupFracHits = 0;
   mupdxyVtx = 0;
   mupdzVtx = 0;
   mupCat = 0;
   mupNPixHits = 0;
   mupNPixLayers = 0;
   mupNTrkHits = 0;
   mupNTrkLayers = 0;
   mupNMuonHits = 0;
   mupNMatchStation = 0;
   mupTrig = 0;
   mup2LastTrigFilter = 0;
   kstTrkmHighPurity = 0;
   kstTrkmCL = 0;
   kstTrkmNormChi2 = 0;
   kstTrkmPx = 0;
   kstTrkmPy = 0;
   kstTrkmPz = 0;
   kstTrkmDCAVtx = 0;
   kstTrkmDCAVtxE = 0;
   kstTrkmDCABS = 0;
   kstTrkmDCABSE = 0;
   kstTrkmFracHits = 0;
   kstTrkmdxyVtx = 0;
   kstTrkmdzVtx = 0;
   kstTrkmNPixHits = 0;
   kstTrkmNPixLayers = 0;
   kstTrkmNTrkHits = 0;
   kstTrkmNTrkLayers = 0;
   kstTrkmMuMatch = 0;
   kstTrkpHighPurity = 0;
   kstTrkpCL = 0;
   kstTrkpNormChi2 = 0;
   kstTrkpPx = 0;
   kstTrkpPy = 0;
   kstTrkpPz = 0;
   kstTrkpDCAVtx = 0;
   kstTrkpDCAVtxE = 0;
   kstTrkpDCABS = 0;
   kstTrkpDCABSE = 0;
   kstTrkpFracHits = 0;
   kstTrkpdxyVtx = 0;
   kstTrkpdzVtx = 0;
   kstTrkpNPixHits = 0;
   kstTrkpNPixLayers = 0;
   kstTrkpNTrkHits = 0;
   kstTrkpNTrkLayers = 0;
   kstTrkpMuMatch = 0;
   truthMatchSignal = 0;
   truthMatchMum = 0;
   truthMatchMup = 0;
   truthMatchTrkm = 0;
   truthMatchTrkp = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runN", &runN, &b_runN);
   fChain->SetBranchAddress("eventN", &eventN, &b_eventN);
   fChain->SetBranchAddress("recoVtxN", &recoVtxN, &b_recoVtxN);
   fChain->SetBranchAddress("evWeight", &evWeight, &b_evWeight);
   fChain->SetBranchAddress("TrigTable", &TrigTable, &b_TrigTable);
   fChain->SetBranchAddress("TrigPrescales", &TrigPrescales, &b_TrigPrescales);
   fChain->SetBranchAddress("nB", &nB, &b_nB);
   fChain->SetBranchAddress("priVtxCL", &priVtxCL, &b_priVtxCL);
   fChain->SetBranchAddress("priVtxX", &priVtxX, &b_priVtxX);
   fChain->SetBranchAddress("priVtxY", &priVtxY, &b_priVtxY);
   fChain->SetBranchAddress("priVtxZ", &priVtxZ, &b_priVtxZ);
   fChain->SetBranchAddress("bsX", &bsX, &b_bsX);
   fChain->SetBranchAddress("bsY", &bsY, &b_bsY);
   fChain->SetBranchAddress("bMass", &bMass, &b_bMass);
   fChain->SetBranchAddress("bMassE", &bMassE, &b_bMassE);
   fChain->SetBranchAddress("bBarMass", &bBarMass, &b_bBarMass);
   fChain->SetBranchAddress("bBarMassE", &bBarMassE, &b_bBarMassE);
   fChain->SetBranchAddress("bPx", &bPx, &b_bPx);
   fChain->SetBranchAddress("bPy", &bPy, &b_bPy);
   fChain->SetBranchAddress("bPz", &bPz, &b_bPz);
   fChain->SetBranchAddress("bunchXingMC", &bunchXingMC, &b_bunchXingMC);
   fChain->SetBranchAddress("numInteractionsMC", &numInteractionsMC, &b_numInteractionsMC);
   fChain->SetBranchAddress("trueNumInteractionsMC", &trueNumInteractionsMC, &b_trueNumInteractionsMC);
   fChain->SetBranchAddress("bVtxCL", &bVtxCL, &b_bVtxCL);
   fChain->SetBranchAddress("bVtxX", &bVtxX, &b_bVtxX);
   fChain->SetBranchAddress("bVtxY", &bVtxY, &b_bVtxY);
   fChain->SetBranchAddress("bVtxZ", &bVtxZ, &b_bVtxZ);
   fChain->SetBranchAddress("bCosAlphaVtx", &bCosAlphaVtx, &b_bCosAlphaVtx);
   fChain->SetBranchAddress("bCosAlphaVtxE", &bCosAlphaVtxE, &b_bCosAlphaVtxE);
   fChain->SetBranchAddress("bCosAlphaBS", &bCosAlphaBS, &b_bCosAlphaBS);
   fChain->SetBranchAddress("bCosAlphaBSE", &bCosAlphaBSE, &b_bCosAlphaBSE);
   fChain->SetBranchAddress("bLVtx", &bLVtx, &b_bLVtx);
   fChain->SetBranchAddress("bLVtxE", &bLVtxE, &b_bLVtxE);
   fChain->SetBranchAddress("bLBS", &bLBS, &b_bLBS);
   fChain->SetBranchAddress("bLBSE", &bLBSE, &b_bLBSE);
   fChain->SetBranchAddress("bDCAVtx", &bDCAVtx, &b_bDCAVtx);
   fChain->SetBranchAddress("bDCAVtxE", &bDCAVtxE, &b_bDCAVtxE);
   fChain->SetBranchAddress("bDCABS", &bDCABS, &b_bDCABS);
   fChain->SetBranchAddress("bDCABSE", &bDCABSE, &b_bDCABSE);
   fChain->SetBranchAddress("bctauPVBS", &bctauPVBS, &b_bctauPVBS);
   fChain->SetBranchAddress("bctauPVBSE", &bctauPVBSE, &b_bctauPVBSE);
   fChain->SetBranchAddress("kstMass", &kstMass, &b_kstMass);
   fChain->SetBranchAddress("kstMassE", &kstMassE, &b_kstMassE);
   fChain->SetBranchAddress("kstBarMass", &kstBarMass, &b_kstBarMass);
   fChain->SetBranchAddress("kstBarMassE", &kstBarMassE, &b_kstBarMassE);
   fChain->SetBranchAddress("kstPx", &kstPx, &b_kstPx);
   fChain->SetBranchAddress("kstPy", &kstPy, &b_kstPy);
   fChain->SetBranchAddress("kstPz", &kstPz, &b_kstPz);
   fChain->SetBranchAddress("kstVtxCL", &kstVtxCL, &b_kstVtxCL);
   fChain->SetBranchAddress("kstVtxX", &kstVtxX, &b_kstVtxX);
   fChain->SetBranchAddress("kstVtxY", &kstVtxY, &b_kstVtxY);
   fChain->SetBranchAddress("kstVtxZ", &kstVtxZ, &b_kstVtxZ);
   fChain->SetBranchAddress("mumuMass", &mumuMass, &b_mumuMass);
   fChain->SetBranchAddress("mumuMassE", &mumuMassE, &b_mumuMassE);
   fChain->SetBranchAddress("mumuPx", &mumuPx, &b_mumuPx);
   fChain->SetBranchAddress("mumuPy", &mumuPy, &b_mumuPy);
   fChain->SetBranchAddress("mumuPz", &mumuPz, &b_mumuPz);
   fChain->SetBranchAddress("mumuVtxCL", &mumuVtxCL, &b_mumuVtxCL);
   fChain->SetBranchAddress("mumuVtxX", &mumuVtxX, &b_mumuVtxX);
   fChain->SetBranchAddress("mumuVtxY", &mumuVtxY, &b_mumuVtxY);
   fChain->SetBranchAddress("mumuVtxZ", &mumuVtxZ, &b_mumuVtxZ);
   fChain->SetBranchAddress("mumuCosAlphaBS", &mumuCosAlphaBS, &b_mumuCosAlphaBS);
   fChain->SetBranchAddress("mumuCosAlphaBSE", &mumuCosAlphaBSE, &b_mumuCosAlphaBSE);
   fChain->SetBranchAddress("mumuLBS", &mumuLBS, &b_mumuLBS);
   fChain->SetBranchAddress("mumuLBSE", &mumuLBSE, &b_mumuLBSE);
   fChain->SetBranchAddress("mumuDCA", &mumuDCA, &b_mumuDCA);
   fChain->SetBranchAddress("mumHighPurity", &mumHighPurity, &b_mumHighPurity);
   fChain->SetBranchAddress("mumCL", &mumCL, &b_mumCL);
   fChain->SetBranchAddress("mumNormChi2", &mumNormChi2, &b_mumNormChi2);
   fChain->SetBranchAddress("mumPx", &mumPx, &b_mumPx);
   fChain->SetBranchAddress("mumPy", &mumPy, &b_mumPy);
   fChain->SetBranchAddress("mumPz", &mumPz, &b_mumPz);
   fChain->SetBranchAddress("mumDCAVtx", &mumDCAVtx, &b_mumDCAVtx);
   fChain->SetBranchAddress("mumDCAVtxE", &mumDCAVtxE, &b_mumDCAVtxE);
   fChain->SetBranchAddress("mumDCABS", &mumDCABS, &b_mumDCABS);
   fChain->SetBranchAddress("mumDCABSE", &mumDCABSE, &b_mumDCABSE);
   fChain->SetBranchAddress("mumKinkChi2", &mumKinkChi2, &b_mumKinkChi2);
   fChain->SetBranchAddress("mumFracHits", &mumFracHits, &b_mumFracHits);
   fChain->SetBranchAddress("mumdxyVtx", &mumdxyVtx, &b_mumdxyVtx);
   fChain->SetBranchAddress("mumdzVtx", &mumdzVtx, &b_mumdzVtx);
   fChain->SetBranchAddress("mumCat", &mumCat, &b_mumCat);
   fChain->SetBranchAddress("mumNPixHits", &mumNPixHits, &b_mumNPixHits);
   fChain->SetBranchAddress("mumNPixLayers", &mumNPixLayers, &b_mumNPixLayers);
   fChain->SetBranchAddress("mumNTrkHits", &mumNTrkHits, &b_mumNTrkHits);
   fChain->SetBranchAddress("mumNTrkLayers", &mumNTrkLayers, &b_mumNTrkLayers);
   fChain->SetBranchAddress("mumNMuonHits", &mumNMuonHits, &b_mumNMuonHits);
   fChain->SetBranchAddress("mumNMatchStation", &mumNMatchStation, &b_mumNMatchStation);
   fChain->SetBranchAddress("mumTrig", &mumTrig, &b_mumTrig);
   fChain->SetBranchAddress("mum2LastTrigFilter", &mum2LastTrigFilter, &b_mum2LastTrigFilter);
   fChain->SetBranchAddress("mupHighPurity", &mupHighPurity, &b_mupHighPurity);
   fChain->SetBranchAddress("mupCL", &mupCL, &b_mupCL);
   fChain->SetBranchAddress("mupNormChi2", &mupNormChi2, &b_mupNormChi2);
   fChain->SetBranchAddress("mupPx", &mupPx, &b_mupPx);
   fChain->SetBranchAddress("mupPy", &mupPy, &b_mupPy);
   fChain->SetBranchAddress("mupPz", &mupPz, &b_mupPz);
   fChain->SetBranchAddress("mupDCAVtx", &mupDCAVtx, &b_mupDCAVtx);
   fChain->SetBranchAddress("mupDCAVtxE", &mupDCAVtxE, &b_mupDCAVtxE);
   fChain->SetBranchAddress("mupDCABS", &mupDCABS, &b_mupDCABS);
   fChain->SetBranchAddress("mupDCABSE", &mupDCABSE, &b_mupDCABSE);
   fChain->SetBranchAddress("mupKinkChi2", &mupKinkChi2, &b_mupKinkChi2);
   fChain->SetBranchAddress("mupFracHits", &mupFracHits, &b_mupFracHits);
   fChain->SetBranchAddress("mupdxyVtx", &mupdxyVtx, &b_mupdxyVtx);
   fChain->SetBranchAddress("mupdzVtx", &mupdzVtx, &b_mupdzVtx);
   fChain->SetBranchAddress("mupCat", &mupCat, &b_mupCat);
   fChain->SetBranchAddress("mupNPixHits", &mupNPixHits, &b_mupNPixHits);
   fChain->SetBranchAddress("mupNPixLayers", &mupNPixLayers, &b_mupNPixLayers);
   fChain->SetBranchAddress("mupNTrkHits", &mupNTrkHits, &b_mupNTrkHits);
   fChain->SetBranchAddress("mupNTrkLayers", &mupNTrkLayers, &b_mupNTrkLayers);
   fChain->SetBranchAddress("mupNMuonHits", &mupNMuonHits, &b_mupNMuonHits);
   fChain->SetBranchAddress("mupNMatchStation", &mupNMatchStation, &b_mupNMatchStation);
   fChain->SetBranchAddress("mupTrig", &mupTrig, &b_mupTrig);
   fChain->SetBranchAddress("mup2LastTrigFilter", &mup2LastTrigFilter, &b_mup2LastTrigFilter);
   fChain->SetBranchAddress("kstTrkmHighPurity", &kstTrkmHighPurity, &b_kstTrkmHighPurity);
   fChain->SetBranchAddress("kstTrkmCL", &kstTrkmCL, &b_kstTrkmCL);
   fChain->SetBranchAddress("kstTrkmNormChi2", &kstTrkmNormChi2, &b_kstTrkmNormChi2);
   fChain->SetBranchAddress("kstTrkmPx", &kstTrkmPx, &b_kstTrkmPx);
   fChain->SetBranchAddress("kstTrkmPy", &kstTrkmPy, &b_kstTrkmPy);
   fChain->SetBranchAddress("kstTrkmPz", &kstTrkmPz, &b_kstTrkmPz);
   fChain->SetBranchAddress("kstTrkmDCAVtx", &kstTrkmDCAVtx, &b_kstTrkmDCAVtx);
   fChain->SetBranchAddress("kstTrkmDCAVtxE", &kstTrkmDCAVtxE, &b_kstTrkmDCAVtxE);
   fChain->SetBranchAddress("kstTrkmDCABS", &kstTrkmDCABS, &b_kstTrkmDCABS);
   fChain->SetBranchAddress("kstTrkmDCABSE", &kstTrkmDCABSE, &b_kstTrkmDCABSE);
   fChain->SetBranchAddress("kstTrkmFracHits", &kstTrkmFracHits, &b_kstTrkmFracHits);
   fChain->SetBranchAddress("kstTrkmdxyVtx", &kstTrkmdxyVtx, &b_kstTrkmdxyVtx);
   fChain->SetBranchAddress("kstTrkmdzVtx", &kstTrkmdzVtx, &b_kstTrkmdzVtx);
   fChain->SetBranchAddress("kstTrkmNPixHits", &kstTrkmNPixHits, &b_kstTrkmNPixHits);
   fChain->SetBranchAddress("kstTrkmNPixLayers", &kstTrkmNPixLayers, &b_kstTrkmNPixLayers);
   fChain->SetBranchAddress("kstTrkmNTrkHits", &kstTrkmNTrkHits, &b_kstTrkmNTrkHits);
   fChain->SetBranchAddress("kstTrkmNTrkLayers", &kstTrkmNTrkLayers, &b_kstTrkmNTrkLayers);
   fChain->SetBranchAddress("kstTrkmMuMatch", &kstTrkmMuMatch, &b_kstTrkmMuMatch);
   fChain->SetBranchAddress("kstTrkpHighPurity", &kstTrkpHighPurity, &b_kstTrkpHighPurity);
   fChain->SetBranchAddress("kstTrkpCL", &kstTrkpCL, &b_kstTrkpCL);
   fChain->SetBranchAddress("kstTrkpNormChi2", &kstTrkpNormChi2, &b_kstTrkpNormChi2);
   fChain->SetBranchAddress("kstTrkpPx", &kstTrkpPx, &b_kstTrkpPx);
   fChain->SetBranchAddress("kstTrkpPy", &kstTrkpPy, &b_kstTrkpPy);
   fChain->SetBranchAddress("kstTrkpPz", &kstTrkpPz, &b_kstTrkpPz);
   fChain->SetBranchAddress("kstTrkpDCAVtx", &kstTrkpDCAVtx, &b_kstTrkpDCAVtx);
   fChain->SetBranchAddress("kstTrkpDCAVtxE", &kstTrkpDCAVtxE, &b_kstTrkpDCAVtxE);
   fChain->SetBranchAddress("kstTrkpDCABS", &kstTrkpDCABS, &b_kstTrkpDCABS);
   fChain->SetBranchAddress("kstTrkpDCABSE", &kstTrkpDCABSE, &b_kstTrkpDCABSE);
   fChain->SetBranchAddress("kstTrkpFracHits", &kstTrkpFracHits, &b_kstTrkpFracHits);
   fChain->SetBranchAddress("kstTrkpdxyVtx", &kstTrkpdxyVtx, &b_kstTrkpdxyVtx);
   fChain->SetBranchAddress("kstTrkpdzVtx", &kstTrkpdzVtx, &b_kstTrkpdzVtx);
   fChain->SetBranchAddress("kstTrkpNPixHits", &kstTrkpNPixHits, &b_kstTrkpNPixHits);
   fChain->SetBranchAddress("kstTrkpNPixLayers", &kstTrkpNPixLayers, &b_kstTrkpNPixLayers);
   fChain->SetBranchAddress("kstTrkpNTrkHits", &kstTrkpNTrkHits, &b_kstTrkpNTrkHits);
   fChain->SetBranchAddress("kstTrkpNTrkLayers", &kstTrkpNTrkLayers, &b_kstTrkpNTrkLayers);
   fChain->SetBranchAddress("kstTrkpMuMatch", &kstTrkpMuMatch, &b_kstTrkpMuMatch);
   fChain->SetBranchAddress("genSignal", &genSignal, &b_genSignal);
   fChain->SetBranchAddress("genMuMuBG", &genMuMuBG, &b_genMuMuBG);
   fChain->SetBranchAddress("genMuMuBGnTrksm", &genMuMuBGnTrksm, &b_genMuMuBGnTrksm);
   fChain->SetBranchAddress("genMuMuBGnTrksp", &genMuMuBGnTrksp, &b_genMuMuBGnTrksp);
   fChain->SetBranchAddress("genPsiPrompt", &genPsiPrompt, &b_genPsiPrompt);
   fChain->SetBranchAddress("genSignHasFSR", &genSignHasFSR, &b_genSignHasFSR);
   fChain->SetBranchAddress("genSignKstHasFSR", &genSignKstHasFSR, &b_genSignKstHasFSR);
   fChain->SetBranchAddress("genSignPsiHasFSR", &genSignPsiHasFSR, &b_genSignPsiHasFSR);
   fChain->SetBranchAddress("genPriVtxX", &genPriVtxX, &b_genPriVtxX);
   fChain->SetBranchAddress("genPriVtxY", &genPriVtxY, &b_genPriVtxY);
   fChain->SetBranchAddress("genPriVtxZ", &genPriVtxZ, &b_genPriVtxZ);
   fChain->SetBranchAddress("genB0Mass", &genB0Mass, &b_genB0Mass);
   fChain->SetBranchAddress("genB0Px", &genB0Px, &b_genB0Px);
   fChain->SetBranchAddress("genB0Py", &genB0Py, &b_genB0Py);
   fChain->SetBranchAddress("genB0Pz", &genB0Pz, &b_genB0Pz);
   fChain->SetBranchAddress("genB0VtxX", &genB0VtxX, &b_genB0VtxX);
   fChain->SetBranchAddress("genB0VtxY", &genB0VtxY, &b_genB0VtxY);
   fChain->SetBranchAddress("genB0VtxZ", &genB0VtxZ, &b_genB0VtxZ);
   fChain->SetBranchAddress("genKstMass", &genKstMass, &b_genKstMass);
   fChain->SetBranchAddress("genKstPx", &genKstPx, &b_genKstPx);
   fChain->SetBranchAddress("genKstPy", &genKstPy, &b_genKstPy);
   fChain->SetBranchAddress("genKstPz", &genKstPz, &b_genKstPz);
   fChain->SetBranchAddress("genKstVtxX", &genKstVtxX, &b_genKstVtxX);
   fChain->SetBranchAddress("genKstVtxY", &genKstVtxY, &b_genKstVtxY);
   fChain->SetBranchAddress("genKstVtxZ", &genKstVtxZ, &b_genKstVtxZ);
   fChain->SetBranchAddress("genPsiMass", &genPsiMass, &b_genPsiMass);
   fChain->SetBranchAddress("genPsiVtxX", &genPsiVtxX, &b_genPsiVtxX);
   fChain->SetBranchAddress("genPsiVtxY", &genPsiVtxY, &b_genPsiVtxY);
   fChain->SetBranchAddress("genPsiVtxZ", &genPsiVtxZ, &b_genPsiVtxZ);
   fChain->SetBranchAddress("genMumPx", &genMumPx, &b_genMumPx);
   fChain->SetBranchAddress("genMumPy", &genMumPy, &b_genMumPy);
   fChain->SetBranchAddress("genMumPz", &genMumPz, &b_genMumPz);
   fChain->SetBranchAddress("trueMumTriggered", &trueMumTriggered, &b_trueMumTriggered);
   fChain->SetBranchAddress("trueMumInAcceptance", &trueMumInAcceptance, &b_trueMumInAcceptance);
   fChain->SetBranchAddress("genMupPx", &genMupPx, &b_genMupPx);
   fChain->SetBranchAddress("genMupPy", &genMupPy, &b_genMupPy);
   fChain->SetBranchAddress("genMupPz", &genMupPz, &b_genMupPz);
   fChain->SetBranchAddress("trueMupTriggered", &trueMupTriggered, &b_trueMupTriggered);
   fChain->SetBranchAddress("trueMupInAcceptance", &trueMupInAcceptance, &b_trueMupInAcceptance);
   fChain->SetBranchAddress("genKstTrkmPx", &genKstTrkmPx, &b_genKstTrkmPx);
   fChain->SetBranchAddress("genKstTrkmPy", &genKstTrkmPy, &b_genKstTrkmPy);
   fChain->SetBranchAddress("genKstTrkmPz", &genKstTrkmPz, &b_genKstTrkmPz);
   fChain->SetBranchAddress("genKstTrkpPx", &genKstTrkpPx, &b_genKstTrkpPx);
   fChain->SetBranchAddress("genKstTrkpPy", &genKstTrkpPy, &b_genKstTrkpPy);
   fChain->SetBranchAddress("genKstTrkpPz", &genKstTrkpPz, &b_genKstTrkpPz);
   fChain->SetBranchAddress("truthMatchSignal", &truthMatchSignal, &b_truthMatchSignal);
   fChain->SetBranchAddress("truthMatchMum", &truthMatchMum, &b_truthMatchMum);
   fChain->SetBranchAddress("truthMatchMup", &truthMatchMup, &b_truthMatchMup);
   fChain->SetBranchAddress("truthMatchTrkm", &truthMatchTrkm, &b_truthMatchTrkm);
   fChain->SetBranchAddress("truthMatchTrkp", &truthMatchTrkp, &b_truthMatchTrkp);
   fChain->SetBranchAddress("evWeightE2", &evWeightE2, &b_evWeightE2);
   fChain->SetBranchAddress("B0MassArb", &B0MassArb, &b_B0MassArb);
   fChain->SetBranchAddress("B0pT", &B0pT, &b_B0pT);
   fChain->SetBranchAddress("B0Eta", &B0Eta, &b_B0Eta);
   fChain->SetBranchAddress("B0Phi", &B0Phi, &b_B0Phi);
   fChain->SetBranchAddress("B0notB0bar", &B0notB0bar, &b_B0notB0bar);
   fChain->SetBranchAddress("TrigCat", &TrigCat, &b_TrigCat);
   fChain->SetBranchAddress("CosThetaKArb", &CosThetaKArb, &b_CosThetaKArb);
   fChain->SetBranchAddress("CosThetaMuArb", &CosThetaMuArb, &b_CosThetaMuArb);
   fChain->SetBranchAddress("PhiKstMuMuPlaneArb", &PhiKstMuMuPlaneArb, &b_PhiKstMuMuPlaneArb);
}

Bool_t SingleBdToKstarMuMuSelector_xCheck2011::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef SingleBdToKstarMuMuSelector_xCheck2011_cxx
