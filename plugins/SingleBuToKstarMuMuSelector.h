//////////////////////////////////////////////////////////
// [2013-07-18 Thu 10:52] Update based on BuToKstarMuMu_v2 
// 
// This class was automatically generated on
// Tue Feb 26 11:42:08 2013 by ROOT version 5.27/06b
// from TTree tree/BToKstarMuMu
// found on file: /afs/cern.ch/user/x/xshi/work/cms/afb/dat/ntp/data/MuOnia/Run2011A_May10ReReco_v1_10_6/BToKstarMuMu_692_1_H0J.root
//////////////////////////////////////////////////////////

#ifndef SingleBuToKstarMuMuSelector_h
#define SingleBuToKstarMuMuSelector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

using namespace std;  

class SingleBuToKstarMuMuSelector : public TSelector {
public :

  TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   UInt_t          run;
   UInt_t          event;
   UInt_t          lumiblock;
   UInt_t          nprivtx;
   vector<string>  *triggernames;
   vector<int>     *triggerprescales;
   vector<double>  *mumdcabs;
   vector<double>  *mumdcabserr;
   vector<double>  *mumpx;
   vector<double>  *mumpy;
   vector<double>  *mumpz;
   vector<double>  *mummass;
   vector<double>  *mummasserr;
   vector<double>  *mupdcabs;
   vector<double>  *mupdcabserr;
   vector<double>  *muppx;
   vector<double>  *muppy;
   vector<double>  *muppz;
   vector<double>  *mupmass;
   vector<double>  *mupmasserr;
   vector<double>  *mumudca;
   vector<double>  *mumuvtxcl;
   vector<double>  *mumulsbs;
   vector<double>  *mumulsbserr;
   vector<double>  *mumucosalphabs;
   vector<double>  *mumucosalphabserr;
   vector<bool>    *mumisgoodmuon;
   vector<bool>    *mupisgoodmuon;
   vector<int>     *mumnpixhits;
   vector<int>     *mupnpixhits;
   vector<int>     *mumnpixlayers;
   vector<int>     *mupnpixlayers;
   vector<int>     *mumntrkhits;
   vector<int>     *mupntrkhits;
   vector<int>     *mumntrklayers;
   vector<int>     *mupntrklayers;
   vector<double>  *mumnormchi2;
   vector<double>  *mupnormchi2;
   vector<double>  *mumdxyvtx;
   vector<double>  *mupdxyvtx;
   vector<double>  *mumdzvtx;
   vector<double>  *mupdzvtx;
   vector<string>  *mumtriglastfilter;
   vector<string>  *muptriglastfilter;
   vector<int>     *trkchg;
   vector<double>  *trkpx;
   vector<double>  *trkpy;
   vector<double>  *trkpz;
   vector<double>  *trkmass;
   vector<double>  *trkmasserr;
   vector<double>  *pimpx;
   vector<double>  *pimpy;
   vector<double>  *pimpz;
   vector<double>  *pimmass;
   vector<double>  *pimd0;
   vector<double>  *pimd0err;
   vector<double>  *pippx;
   vector<double>  *pippy;
   vector<double>  *pippz;
   vector<double>  *pipmass;
   vector<double>  *pipd0;
   vector<double>  *pipd0err;
   vector<double>  *kspx;
   vector<double>  *kspy;
   vector<double>  *kspz;
   vector<double>  *ksmass;
   vector<double>  *ksmasserr;
   vector<double>  *ksvtxx;
   vector<double>  *ksvtxy;
   vector<double>  *ksvtxz;
   vector<double>  *ksvtxcl;
   vector<double>  *kslsbs;
   vector<double>  *kslsbserr;
   vector<int>     *bchg;
   vector<double>  *bpx;
   vector<double>  *bpxerr;
   vector<double>  *bpy;
   vector<double>  *bpyerr;
   vector<double>  *bpz;
   vector<double>  *bpzerr;
   vector<double>  *bmass;
   vector<double>  *bmasserr;
   vector<double>  *bvtxcl;
   vector<double>  *bvtxx;
   vector<double>  *bvtxxerr;
   vector<double>  *bvtxy;
   vector<double>  *bvtxyerr;
   vector<double>  *bvtxz;
   vector<double>  *bvtxzerr;
   vector<double>  *bcosalphabs;
   vector<double>  *bcosalphabserr;
   vector<double>  *blsbs;
   vector<double>  *blsbserr;
   vector<double>  *bctau;
   vector<double>  *bctauerr;
   Int_t           genbchg;
   Double_t        genbpx;
   Double_t        genbpy;
   Double_t        genbpz;
   Double_t        genkstpx;
   Double_t        genkstpy;
   Double_t        genkstpz;
   Double_t        genkspx;
   Double_t        genkspy;
   Double_t        genkspz;
   Double_t        genksvtxx;
   Double_t        genksvtxy;
   Double_t        genksvtxz;
   Int_t           gentrkchg;
   Double_t        gentrkpx;
   Double_t        gentrkpy;
   Double_t        gentrkpz;
   Double_t        genmumpx;
   Double_t        genmumpy;
   Double_t        genmumpz;
   Double_t        genmuppx;
   Double_t        genmuppy;
   Double_t        genmuppz;
   Double_t        genpippx;
   Double_t        genpippy;
   Double_t        genpippz;
   Double_t        genpimpx;
   Double_t        genpimpy;
   Double_t        genpimpz;
   vector<bool>    *istruemum;
   vector<bool>    *istruemup;
   vector<bool>    *istrueks;
   vector<bool>    *istruetrk;
   vector<bool>    *istruebu;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumiblock;   //!
   TBranch        *b_nprivtx;   //!
   TBranch        *b_triggernames;   //!
   TBranch        *b_triggerprescales;   //!
   TBranch        *b_mumdcabs;   //!
   TBranch        *b_mumdcabserr;   //!
   TBranch        *b_mumpx;   //!
   TBranch        *b_mumpy;   //!
   TBranch        *b_mumpz;   //!
   TBranch        *b_mummass;   //!
   TBranch        *b_mummasserr;   //!
   TBranch        *b_mupdcabs;   //!
   TBranch        *b_mupdcabserr;   //!
   TBranch        *b_muppx;   //!
   TBranch        *b_muppy;   //!
   TBranch        *b_muppz;   //!
   TBranch        *b_mupmass;   //!
   TBranch        *b_mupmasserr;   //!
   TBranch        *b_mumudca;   //!
   TBranch        *b_mumuvtxcl;   //!
   TBranch        *b_mumulsbs;   //!
   TBranch        *b_mumulsbserr;   //!
   TBranch        *b_mumucosalphabs;   //!
   TBranch        *b_mumucosalphabserr;   //!
   TBranch        *b_mumisgoodmuon;   //!
   TBranch        *b_mupisgoodmuon;   //!
   TBranch        *b_mumnpixhits;   //!
   TBranch        *b_mupnpixhits;   //!
   TBranch        *b_mumnpixlayers;   //!
   TBranch        *b_mupnpixlayers;   //!
   TBranch        *b_mumntrkhits;   //!
   TBranch        *b_mupntrkhits;   //!
   TBranch        *b_mumntrklayers;   //!
   TBranch        *b_mupntrklayers;   //!
   TBranch        *b_mumnormchi2;   //!
   TBranch        *b_mupnormchi2;   //!
   TBranch        *b_mumdxyvtx;   //!
   TBranch        *b_mupdxyvtx;   //!
   TBranch        *b_mumdzvtx;   //!
   TBranch        *b_mupdzvtx;   //!
   TBranch        *b_mumtriglastfilter;   //!
   TBranch        *b_muptriglastfilter;   //!
   TBranch        *b_trkchg;   //!
   TBranch        *b_trkpx;   //!
   TBranch        *b_trkpy;   //!
   TBranch        *b_trkpz;   //!
   TBranch        *b_trkmass;   //!
   TBranch        *b_trkmasserr;   //!
   TBranch        *b_pimpx;   //!
   TBranch        *b_pimpy;   //!
   TBranch        *b_pimpz;   //!
   TBranch        *b_pimmass;   //!
   TBranch        *b_pimd0;   //!
   TBranch        *b_pimd0err;   //!
   TBranch        *b_pippx;   //!
   TBranch        *b_pippy;   //!
   TBranch        *b_pippz;   //!
   TBranch        *b_pipmass;   //!
   TBranch        *b_pipd0;   //!
   TBranch        *b_pipd0err;   //!
   TBranch        *b_kspx;   //!
   TBranch        *b_kspy;   //!
   TBranch        *b_kspz;   //!
   TBranch        *b_ksmass;   //!
   TBranch        *b_ksmasserr;   //!
   TBranch        *b_ksvtxx;   //!
   TBranch        *b_ksvtxy;   //!
   TBranch        *b_ksvtxz;   //!
   TBranch        *b_ksvtxcl;   //!
   TBranch        *b_kslsbs;   //!
   TBranch        *b_kslsbserr;   //!
   TBranch        *b_bchg;   //!
   TBranch        *b_bpx;   //!
   TBranch        *b_bpxerr;   //!
   TBranch        *b_bpy;   //!
   TBranch        *b_bpyerr;   //!
   TBranch        *b_bpz;   //!
   TBranch        *b_bpzerr;   //!
   TBranch        *b_bmass;   //!
   TBranch        *b_bmasserr;   //!
   TBranch        *b_bvtxcl;   //!
   TBranch        *b_bvtxx;   //!
   TBranch        *b_bvtxxerr;   //!
   TBranch        *b_bvtxy;   //!
   TBranch        *b_bvtxyerr;   //!
   TBranch        *b_bvtxz;   //!
   TBranch        *b_bvtxzerr;   //!
   TBranch        *b_bcosalphabs;   //!
   TBranch        *b_bcosalphabserr;   //!
   TBranch        *b_blsbs;   //!
   TBranch        *b_blsbserr;   //!
   TBranch        *b_bctau;   //!
   TBranch        *b_bctauerr;   //!
   TBranch        *b_genbchg;   //!
   TBranch        *b_genbpx;   //!
   TBranch        *b_genbpy;   //!
   TBranch        *b_genbpz;   //!
   TBranch        *b_genkstpx;   //!
   TBranch        *b_genkstpy;   //!
   TBranch        *b_genkstpz;   //!
   TBranch        *b_genkspx;   //!
   TBranch        *b_genkspy;   //!
   TBranch        *b_genkspz;   //!
   TBranch        *b_genksvtxx;   //!
   TBranch        *b_genksvtxy;   //!
   TBranch        *b_genksvtxz;   //!
   TBranch        *b_gentrkchg;   //!
   TBranch        *b_gentrkpx;   //!
   TBranch        *b_gentrkpy;   //!
   TBranch        *b_gentrkpz;   //!
   TBranch        *b_genmumpx;   //!
   TBranch        *b_genmumpy;   //!
   TBranch        *b_genmumpz;   //!
   TBranch        *b_genmuppx;   //!
   TBranch        *b_genmuppy;   //!
   TBranch        *b_genmuppz;   //!
   TBranch        *b_genpippx;   //!
   TBranch        *b_genpippy;   //!
   TBranch        *b_genpippz;   //!
   TBranch        *b_genpimpx;   //!
   TBranch        *b_genpimpy;   //!
   TBranch        *b_genpimpz;   //!
   TBranch        *b_istruemum;   //!
   TBranch        *b_istruemup;   //!
   TBranch        *b_istrueks;   //!
   TBranch        *b_istruetrk;   //!
   TBranch        *b_istruebu;   //!

   SingleBuToKstarMuMuSelector(TTree * /*tree*/ =0) { }
   virtual ~SingleBuToKstarMuMuSelector() { }
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

   // ClassDef(SingleBuToKstarMuMuSelector,0);
   int SelectB(string); 
   bool HasGoodDimuon(int); 
   void SaveMuMu(int); 
   // void SaveKstar(int); 
   void SaveB(int); 
   //  double GetKstarMass(int); 
};

#endif

#ifdef SingleBuToKstarMuMuSelector_cxx
void SingleBuToKstarMuMuSelector::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   triggernames = 0;
   triggerprescales = 0;
   mumdcabs = 0;
   mumdcabserr = 0;
   mumpx = 0;
   mumpy = 0;
   mumpz = 0;
   mummass = 0;
   mummasserr = 0;
   mupdcabs = 0;
   mupdcabserr = 0;
   muppx = 0;
   muppy = 0;
   muppz = 0;
   mupmass = 0;
   mupmasserr = 0;
   mumudca = 0;
   mumuvtxcl = 0;
   mumulsbs = 0;
   mumulsbserr = 0;
   mumucosalphabs = 0;
   mumucosalphabserr = 0;
   mumisgoodmuon = 0;
   mupisgoodmuon = 0;
   mumnpixhits = 0;
   mupnpixhits = 0;
   mumnpixlayers = 0;
   mupnpixlayers = 0;
   mumntrkhits = 0;
   mupntrkhits = 0;
   mumntrklayers = 0;
   mupntrklayers = 0;
   mumnormchi2 = 0;
   mupnormchi2 = 0;
   mumdxyvtx = 0;
   mupdxyvtx = 0;
   mumdzvtx = 0;
   mupdzvtx = 0;
   mumtriglastfilter = 0;
   muptriglastfilter = 0;
   trkchg = 0;
   trkpx = 0;
   trkpy = 0;
   trkpz = 0;
   trkmass = 0;
   trkmasserr = 0;
   pimpx = 0;
   pimpy = 0;
   pimpz = 0;
   pimmass = 0;
   pimd0 = 0;
   pimd0err = 0;
   pippx = 0;
   pippy = 0;
   pippz = 0;
   pipmass = 0;
   pipd0 = 0;
   pipd0err = 0;
   kspx = 0;
   kspy = 0;
   kspz = 0;
   ksmass = 0;
   ksmasserr = 0;
   ksvtxx = 0;
   ksvtxy = 0;
   ksvtxz = 0;
   ksvtxcl = 0;
   kslsbs = 0;
   kslsbserr = 0;
   bchg = 0;
   bpx = 0;
   bpxerr = 0;
   bpy = 0;
   bpyerr = 0;
   bpz = 0;
   bpzerr = 0;
   bmass = 0;
   bmasserr = 0;
   bvtxcl = 0;
   bvtxx = 0;
   bvtxxerr = 0;
   bvtxy = 0;
   bvtxyerr = 0;
   bvtxz = 0;
   bvtxzerr = 0;
   bcosalphabs = 0;
   bcosalphabserr = 0;
   blsbs = 0;
   blsbserr = 0;
   bctau = 0;
   bctauerr = 0;
   istruemum = 0;
   istruemup = 0;
   istrueks = 0;
   istruetrk = 0;
   istruebu = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumiblock", &lumiblock, &b_lumiblock);
   fChain->SetBranchAddress("nprivtx", &nprivtx, &b_nprivtx);
   fChain->SetBranchAddress("triggernames", &triggernames, &b_triggernames);
   fChain->SetBranchAddress("triggerprescales", &triggerprescales, &b_triggerprescales);
   fChain->SetBranchAddress("mumdcabs", &mumdcabs, &b_mumdcabs);
   fChain->SetBranchAddress("mumdcabserr", &mumdcabserr, &b_mumdcabserr);
   fChain->SetBranchAddress("mumpx", &mumpx, &b_mumpx);
   fChain->SetBranchAddress("mumpy", &mumpy, &b_mumpy);
   fChain->SetBranchAddress("mumpz", &mumpz, &b_mumpz);
   fChain->SetBranchAddress("mummass", &mummass, &b_mummass);
   fChain->SetBranchAddress("mummasserr", &mummasserr, &b_mummasserr);
   fChain->SetBranchAddress("mupdcabs", &mupdcabs, &b_mupdcabs);
   fChain->SetBranchAddress("mupdcabserr", &mupdcabserr, &b_mupdcabserr);
   fChain->SetBranchAddress("muppx", &muppx, &b_muppx);
   fChain->SetBranchAddress("muppy", &muppy, &b_muppy);
   fChain->SetBranchAddress("muppz", &muppz, &b_muppz);
   fChain->SetBranchAddress("mupmass", &mupmass, &b_mupmass);
   fChain->SetBranchAddress("mupmasserr", &mupmasserr, &b_mupmasserr);
   fChain->SetBranchAddress("mumudca", &mumudca, &b_mumudca);
   fChain->SetBranchAddress("mumuvtxcl", &mumuvtxcl, &b_mumuvtxcl);
   fChain->SetBranchAddress("mumulsbs", &mumulsbs, &b_mumulsbs);
   fChain->SetBranchAddress("mumulsbserr", &mumulsbserr, &b_mumulsbserr);
   fChain->SetBranchAddress("mumucosalphabs", &mumucosalphabs, &b_mumucosalphabs);
   fChain->SetBranchAddress("mumucosalphabserr", &mumucosalphabserr, &b_mumucosalphabserr);
   fChain->SetBranchAddress("mumisgoodmuon", &mumisgoodmuon, &b_mumisgoodmuon);
   fChain->SetBranchAddress("mupisgoodmuon", &mupisgoodmuon, &b_mupisgoodmuon);
   fChain->SetBranchAddress("mumnpixhits", &mumnpixhits, &b_mumnpixhits);
   fChain->SetBranchAddress("mupnpixhits", &mupnpixhits, &b_mupnpixhits);
   fChain->SetBranchAddress("mumnpixlayers", &mumnpixlayers, &b_mumnpixlayers);
   fChain->SetBranchAddress("mupnpixlayers", &mupnpixlayers, &b_mupnpixlayers);
   fChain->SetBranchAddress("mumntrkhits", &mumntrkhits, &b_mumntrkhits);
   fChain->SetBranchAddress("mupntrkhits", &mupntrkhits, &b_mupntrkhits);
   fChain->SetBranchAddress("mumntrklayers", &mumntrklayers, &b_mumntrklayers);
   fChain->SetBranchAddress("mupntrklayers", &mupntrklayers, &b_mupntrklayers);
   fChain->SetBranchAddress("mumnormchi2", &mumnormchi2, &b_mumnormchi2);
   fChain->SetBranchAddress("mupnormchi2", &mupnormchi2, &b_mupnormchi2);
   fChain->SetBranchAddress("mumdxyvtx", &mumdxyvtx, &b_mumdxyvtx);
   fChain->SetBranchAddress("mupdxyvtx", &mupdxyvtx, &b_mupdxyvtx);
   fChain->SetBranchAddress("mumdzvtx", &mumdzvtx, &b_mumdzvtx);
   fChain->SetBranchAddress("mupdzvtx", &mupdzvtx, &b_mupdzvtx);
   fChain->SetBranchAddress("mumtriglastfilter", &mumtriglastfilter, &b_mumtriglastfilter);
   fChain->SetBranchAddress("muptriglastfilter", &muptriglastfilter, &b_muptriglastfilter);
   fChain->SetBranchAddress("trkchg", &trkchg, &b_trkchg);
   fChain->SetBranchAddress("trkpx", &trkpx, &b_trkpx);
   fChain->SetBranchAddress("trkpy", &trkpy, &b_trkpy);
   fChain->SetBranchAddress("trkpz", &trkpz, &b_trkpz);
   fChain->SetBranchAddress("trkmass", &trkmass, &b_trkmass);
   fChain->SetBranchAddress("trkmasserr", &trkmasserr, &b_trkmasserr);
   fChain->SetBranchAddress("pimpx", &pimpx, &b_pimpx);
   fChain->SetBranchAddress("pimpy", &pimpy, &b_pimpy);
   fChain->SetBranchAddress("pimpz", &pimpz, &b_pimpz);
   fChain->SetBranchAddress("pimmass", &pimmass, &b_pimmass);
   fChain->SetBranchAddress("pimd0", &pimd0, &b_pimd0);
   fChain->SetBranchAddress("pimd0err", &pimd0err, &b_pimd0err);
   fChain->SetBranchAddress("pippx", &pippx, &b_pippx);
   fChain->SetBranchAddress("pippy", &pippy, &b_pippy);
   fChain->SetBranchAddress("pippz", &pippz, &b_pippz);
   fChain->SetBranchAddress("pipmass", &pipmass, &b_pipmass);
   fChain->SetBranchAddress("pipd0", &pipd0, &b_pipd0);
   fChain->SetBranchAddress("pipd0err", &pipd0err, &b_pipd0err);
   fChain->SetBranchAddress("kspx", &kspx, &b_kspx);
   fChain->SetBranchAddress("kspy", &kspy, &b_kspy);
   fChain->SetBranchAddress("kspz", &kspz, &b_kspz);
   fChain->SetBranchAddress("ksmass", &ksmass, &b_ksmass);
   fChain->SetBranchAddress("ksmasserr", &ksmasserr, &b_ksmasserr);
   fChain->SetBranchAddress("ksvtxx", &ksvtxx, &b_ksvtxx);
   fChain->SetBranchAddress("ksvtxy", &ksvtxy, &b_ksvtxy);
   fChain->SetBranchAddress("ksvtxz", &ksvtxz, &b_ksvtxz);
   fChain->SetBranchAddress("ksvtxcl", &ksvtxcl, &b_ksvtxcl);
   fChain->SetBranchAddress("kslsbs", &kslsbs, &b_kslsbs);
   fChain->SetBranchAddress("kslsbserr", &kslsbserr, &b_kslsbserr);
   fChain->SetBranchAddress("bchg", &bchg, &b_bchg);
   fChain->SetBranchAddress("bpx", &bpx, &b_bpx);
   fChain->SetBranchAddress("bpxerr", &bpxerr, &b_bpxerr);
   fChain->SetBranchAddress("bpy", &bpy, &b_bpy);
   fChain->SetBranchAddress("bpyerr", &bpyerr, &b_bpyerr);
   fChain->SetBranchAddress("bpz", &bpz, &b_bpz);
   fChain->SetBranchAddress("bpzerr", &bpzerr, &b_bpzerr);
   fChain->SetBranchAddress("bmass", &bmass, &b_bmass);
   fChain->SetBranchAddress("bmasserr", &bmasserr, &b_bmasserr);
   fChain->SetBranchAddress("bvtxcl", &bvtxcl, &b_bvtxcl);
   fChain->SetBranchAddress("bvtxx", &bvtxx, &b_bvtxx);
   fChain->SetBranchAddress("bvtxxerr", &bvtxxerr, &b_bvtxxerr);
   fChain->SetBranchAddress("bvtxy", &bvtxy, &b_bvtxy);
   fChain->SetBranchAddress("bvtxyerr", &bvtxyerr, &b_bvtxyerr);
   fChain->SetBranchAddress("bvtxz", &bvtxz, &b_bvtxz);
   fChain->SetBranchAddress("bvtxzerr", &bvtxzerr, &b_bvtxzerr);
   fChain->SetBranchAddress("bcosalphabs", &bcosalphabs, &b_bcosalphabs);
   fChain->SetBranchAddress("bcosalphabserr", &bcosalphabserr, &b_bcosalphabserr);
   fChain->SetBranchAddress("blsbs", &blsbs, &b_blsbs);
   fChain->SetBranchAddress("blsbserr", &blsbserr, &b_blsbserr);
   fChain->SetBranchAddress("bctau", &bctau, &b_bctau);
   fChain->SetBranchAddress("bctauerr", &bctauerr, &b_bctauerr);
   fChain->SetBranchAddress("genbchg", &genbchg, &b_genbchg);
   fChain->SetBranchAddress("genbpx", &genbpx, &b_genbpx);
   fChain->SetBranchAddress("genbpy", &genbpy, &b_genbpy);
   fChain->SetBranchAddress("genbpz", &genbpz, &b_genbpz);
   fChain->SetBranchAddress("genkstpx", &genkstpx, &b_genkstpx);
   fChain->SetBranchAddress("genkstpy", &genkstpy, &b_genkstpy);
   fChain->SetBranchAddress("genkstpz", &genkstpz, &b_genkstpz);
   fChain->SetBranchAddress("genkspx", &genkspx, &b_genkspx);
   fChain->SetBranchAddress("genkspy", &genkspy, &b_genkspy);
   fChain->SetBranchAddress("genkspz", &genkspz, &b_genkspz);
   fChain->SetBranchAddress("genksvtxx", &genksvtxx, &b_genksvtxx);
   fChain->SetBranchAddress("genksvtxy", &genksvtxy, &b_genksvtxy);
   fChain->SetBranchAddress("genksvtxz", &genksvtxz, &b_genksvtxz);
   fChain->SetBranchAddress("gentrkchg", &gentrkchg, &b_gentrkchg);
   fChain->SetBranchAddress("gentrkpx", &gentrkpx, &b_gentrkpx);
   fChain->SetBranchAddress("gentrkpy", &gentrkpy, &b_gentrkpy);
   fChain->SetBranchAddress("gentrkpz", &gentrkpz, &b_gentrkpz);
   fChain->SetBranchAddress("genmumpx", &genmumpx, &b_genmumpx);
   fChain->SetBranchAddress("genmumpy", &genmumpy, &b_genmumpy);
   fChain->SetBranchAddress("genmumpz", &genmumpz, &b_genmumpz);
   fChain->SetBranchAddress("genmuppx", &genmuppx, &b_genmuppx);
   fChain->SetBranchAddress("genmuppy", &genmuppy, &b_genmuppy);
   fChain->SetBranchAddress("genmuppz", &genmuppz, &b_genmuppz);
   fChain->SetBranchAddress("genpippx", &genpippx, &b_genpippx);
   fChain->SetBranchAddress("genpippy", &genpippy, &b_genpippy);
   fChain->SetBranchAddress("genpippz", &genpippz, &b_genpippz);
   fChain->SetBranchAddress("genpimpx", &genpimpx, &b_genpimpx);
   fChain->SetBranchAddress("genpimpy", &genpimpy, &b_genpimpy);
   fChain->SetBranchAddress("genpimpz", &genpimpz, &b_genpimpz);
   fChain->SetBranchAddress("istruemum", &istruemum, &b_istruemum);
   fChain->SetBranchAddress("istruemup", &istruemup, &b_istruemup);
   fChain->SetBranchAddress("istrueks", &istrueks, &b_istrueks);
   fChain->SetBranchAddress("istruetrk", &istruetrk, &b_istruetrk);
   fChain->SetBranchAddress("istruebu", &istruebu, &b_istruebu);
}

Bool_t SingleBuToKstarMuMuSelector::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef SingleBuToKstarMuMuSelector_cxx
