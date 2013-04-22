//////////////////////////////////////////////////////////
// This class has been automatically generated on
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
   vector<double>  *mupdcabs;
   vector<double>  *mupdcabserr;
   vector<double>  *muppx;
   vector<double>  *muppy;
   vector<double>  *muppz;
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
   vector<double>  *ksvtxx;
   vector<double>  *ksvtxy;
   vector<double>  *ksvtxz;
   vector<double>  *ksvtxcl;
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
   vector<int>     *bmu1chg;
   vector<int>     *bmu2chg;
   vector<int>     *bpi1chg;
   vector<double>  *bmu1px;
   vector<double>  *bmu1py;
   vector<double>  *bmu1pz;
   vector<double>  *bmu2px;
   vector<double>  *bmu2py;
   vector<double>  *bmu2pz;
   vector<double>  *bpi1px;
   vector<double>  *bpi1py;
   vector<double>  *bpi1pz;
   vector<double>  *bkspx;
   vector<double>  *bkspy;
   vector<double>  *bkspz;
   vector<int>     *b3mu1chg;
   vector<int>     *b3mu2chg;
   vector<int>     *b3pi1chg;
   vector<double>  *b3mu1px;
   vector<double>  *b3mu1py;
   vector<double>  *b3mu1pz;
   vector<double>  *b3mu2px;
   vector<double>  *b3mu2py;
   vector<double>  *b3mu2pz;
   vector<double>  *b3pi1px;
   vector<double>  *b3pi1py;
   vector<double>  *b3pi1pz;
   vector<int>     *b3chg;
   vector<double>  *b3px;
   vector<double>  *b3pxerr;
   vector<double>  *b3py;
   vector<double>  *b3pyerr;
   vector<double>  *b3pz;
   vector<double>  *b3pzerr;
   vector<double>  *b3mass;
   vector<double>  *b3masserr;
   vector<double>  *b3vtxcl;
   vector<double>  *b3vtxx;
   vector<double>  *b3vtxxerr;
   vector<double>  *b3vtxy;
   vector<double>  *b3vtxyerr;
   vector<double>  *b3vtxz;
   vector<double>  *b3vtxzerr;
   vector<double>  *b3cosalphabs;
   vector<double>  *b3cosalphabserr;
   vector<double>  *b3lsbs;
   vector<double>  *b3lsbserr;
   vector<double>  *b3ctau;
   vector<double>  *b3ctauerr;

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
   TBranch        *b_mupdcabs;   //!
   TBranch        *b_mupdcabserr;   //!
   TBranch        *b_muppx;   //!
   TBranch        *b_muppy;   //!
   TBranch        *b_muppz;   //!
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
   TBranch        *b_ksvtxx;   //!
   TBranch        *b_ksvtxy;   //!
   TBranch        *b_ksvtxz;   //!
   TBranch        *b_ksvtxcl;   //!
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
   TBranch        *b_bmu1chg;   //!
   TBranch        *b_bmu2chg;   //!
   TBranch        *b_bpi1chg;   //!
   TBranch        *b_bmu1px;   //!
   TBranch        *b_bmu1py;   //!
   TBranch        *b_bmu1pz;   //!
   TBranch        *b_bmu2px;   //!
   TBranch        *b_bmu2py;   //!
   TBranch        *b_bmu2pz;   //!
   TBranch        *b_bpi1px;   //!
   TBranch        *b_bpi1py;   //!
   TBranch        *b_bpi1pz;   //!
   TBranch        *b_bkspx;   //!
   TBranch        *b_bkspy;   //!
   TBranch        *b_bkspz;   //!
   TBranch        *b_b3mu1chg;   //!
   TBranch        *b_b3mu2chg;   //!
   TBranch        *b_b3pi1chg;   //!
   TBranch        *b_b3mu1px;   //!
   TBranch        *b_b3mu1py;   //!
   TBranch        *b_b3mu1pz;   //!
   TBranch        *b_b3mu2px;   //!
   TBranch        *b_b3mu2py;   //!
   TBranch        *b_b3mu2pz;   //!
   TBranch        *b_b3pi1px;   //!
   TBranch        *b_b3pi1py;   //!
   TBranch        *b_b3pi1pz;   //!
   TBranch        *b_b3chg;   //!
   TBranch        *b_b3px;   //!
   TBranch        *b_b3pxerr;   //!
   TBranch        *b_b3py;   //!
   TBranch        *b_b3pyerr;   //!
   TBranch        *b_b3pz;   //!
   TBranch        *b_b3pzerr;   //!
   TBranch        *b_b3mass;   //!
   TBranch        *b_b3masserr;   //!
   TBranch        *b_b3vtxcl;   //!
   TBranch        *b_b3vtxx;   //!
   TBranch        *b_b3vtxxerr;   //!
   TBranch        *b_b3vtxy;   //!
   TBranch        *b_b3vtxyerr;   //!
   TBranch        *b_b3vtxz;   //!
   TBranch        *b_b3vtxzerr;   //!
   TBranch        *b_b3cosalphabs;   //!
   TBranch        *b_b3cosalphabserr;   //!
   TBranch        *b_b3lsbs;   //!
   TBranch        *b_b3lsbserr;   //!
   TBranch        *b_b3ctau;   //!
   TBranch        *b_b3ctauerr;   //!

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
   int SelectB(); 
   bool HasGoodDimuon(); 
   void SaveMuMu(int); 
   void SaveKstar(int); 
   void SaveB(int); 
   string TString_to_string(TString);
   TString get_option_value(string, string);    
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
   mupdcabs = 0;
   mupdcabserr = 0;
   muppx = 0;
   muppy = 0;
   muppz = 0;
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
   ksvtxx = 0;
   ksvtxy = 0;
   ksvtxz = 0;
   ksvtxcl = 0;
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
   bmu1chg = 0;
   bmu2chg = 0;
   bpi1chg = 0;
   bmu1px = 0;
   bmu1py = 0;
   bmu1pz = 0;
   bmu2px = 0;
   bmu2py = 0;
   bmu2pz = 0;
   bpi1px = 0;
   bpi1py = 0;
   bpi1pz = 0;
   bkspx = 0;
   bkspy = 0;
   bkspz = 0;
   b3mu1chg = 0;
   b3mu2chg = 0;
   b3pi1chg = 0;
   b3mu1px = 0;
   b3mu1py = 0;
   b3mu1pz = 0;
   b3mu2px = 0;
   b3mu2py = 0;
   b3mu2pz = 0;
   b3pi1px = 0;
   b3pi1py = 0;
   b3pi1pz = 0;
   b3chg = 0;
   b3px = 0;
   b3pxerr = 0;
   b3py = 0;
   b3pyerr = 0;
   b3pz = 0;
   b3pzerr = 0;
   b3mass = 0;
   b3masserr = 0;
   b3vtxcl = 0;
   b3vtxx = 0;
   b3vtxxerr = 0;
   b3vtxy = 0;
   b3vtxyerr = 0;
   b3vtxz = 0;
   b3vtxzerr = 0;
   b3cosalphabs = 0;
   b3cosalphabserr = 0;
   b3lsbs = 0;
   b3lsbserr = 0;
   b3ctau = 0;
   b3ctauerr = 0;
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
   fChain->SetBranchAddress("mupdcabs", &mupdcabs, &b_mupdcabs);
   fChain->SetBranchAddress("mupdcabserr", &mupdcabserr, &b_mupdcabserr);
   fChain->SetBranchAddress("muppx", &muppx, &b_muppx);
   fChain->SetBranchAddress("muppy", &muppy, &b_muppy);
   fChain->SetBranchAddress("muppz", &muppz, &b_muppz);
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
   fChain->SetBranchAddress("ksvtxx", &ksvtxx, &b_ksvtxx);
   fChain->SetBranchAddress("ksvtxy", &ksvtxy, &b_ksvtxy);
   fChain->SetBranchAddress("ksvtxz", &ksvtxz, &b_ksvtxz);
   fChain->SetBranchAddress("ksvtxcl", &ksvtxcl, &b_ksvtxcl);
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
   fChain->SetBranchAddress("bmu1chg", &bmu1chg, &b_bmu1chg);
   fChain->SetBranchAddress("bmu2chg", &bmu2chg, &b_bmu2chg);
   fChain->SetBranchAddress("bpi1chg", &bpi1chg, &b_bpi1chg);
   fChain->SetBranchAddress("bmu1px", &bmu1px, &b_bmu1px);
   fChain->SetBranchAddress("bmu1py", &bmu1py, &b_bmu1py);
   fChain->SetBranchAddress("bmu1pz", &bmu1pz, &b_bmu1pz);
   fChain->SetBranchAddress("bmu2px", &bmu2px, &b_bmu2px);
   fChain->SetBranchAddress("bmu2py", &bmu2py, &b_bmu2py);
   fChain->SetBranchAddress("bmu2pz", &bmu2pz, &b_bmu2pz);
   fChain->SetBranchAddress("bpi1px", &bpi1px, &b_bpi1px);
   fChain->SetBranchAddress("bpi1py", &bpi1py, &b_bpi1py);
   fChain->SetBranchAddress("bpi1pz", &bpi1pz, &b_bpi1pz);
   fChain->SetBranchAddress("bkspx", &bkspx, &b_bkspx);
   fChain->SetBranchAddress("bkspy", &bkspy, &b_bkspy);
   fChain->SetBranchAddress("bkspz", &bkspz, &b_bkspz);
   fChain->SetBranchAddress("b3mu1chg", &b3mu1chg, &b_b3mu1chg);
   fChain->SetBranchAddress("b3mu2chg", &b3mu2chg, &b_b3mu2chg);
   fChain->SetBranchAddress("b3pi1chg", &b3pi1chg, &b_b3pi1chg);
   fChain->SetBranchAddress("b3mu1px", &b3mu1px, &b_b3mu1px);
   fChain->SetBranchAddress("b3mu1py", &b3mu1py, &b_b3mu1py);
   fChain->SetBranchAddress("b3mu1pz", &b3mu1pz, &b_b3mu1pz);
   fChain->SetBranchAddress("b3mu2px", &b3mu2px, &b_b3mu2px);
   fChain->SetBranchAddress("b3mu2py", &b3mu2py, &b_b3mu2py);
   fChain->SetBranchAddress("b3mu2pz", &b3mu2pz, &b_b3mu2pz);
   fChain->SetBranchAddress("b3pi1px", &b3pi1px, &b_b3pi1px);
   fChain->SetBranchAddress("b3pi1py", &b3pi1py, &b_b3pi1py);
   fChain->SetBranchAddress("b3pi1pz", &b3pi1pz, &b_b3pi1pz);
   fChain->SetBranchAddress("b3chg", &b3chg, &b_b3chg);
   fChain->SetBranchAddress("b3px", &b3px, &b_b3px);
   fChain->SetBranchAddress("b3pxerr", &b3pxerr, &b_b3pxerr);
   fChain->SetBranchAddress("b3py", &b3py, &b_b3py);
   fChain->SetBranchAddress("b3pyerr", &b3pyerr, &b_b3pyerr);
   fChain->SetBranchAddress("b3pz", &b3pz, &b_b3pz);
   fChain->SetBranchAddress("b3pzerr", &b3pzerr, &b_b3pzerr);
   fChain->SetBranchAddress("b3mass", &b3mass, &b_b3mass);
   fChain->SetBranchAddress("b3masserr", &b3masserr, &b_b3masserr);
   fChain->SetBranchAddress("b3vtxcl", &b3vtxcl, &b_b3vtxcl);
   fChain->SetBranchAddress("b3vtxx", &b3vtxx, &b_b3vtxx);
   fChain->SetBranchAddress("b3vtxxerr", &b3vtxxerr, &b_b3vtxxerr);
   fChain->SetBranchAddress("b3vtxy", &b3vtxy, &b_b3vtxy);
   fChain->SetBranchAddress("b3vtxyerr", &b3vtxyerr, &b_b3vtxyerr);
   fChain->SetBranchAddress("b3vtxz", &b3vtxz, &b_b3vtxz);
   fChain->SetBranchAddress("b3vtxzerr", &b3vtxzerr, &b_b3vtxzerr);
   fChain->SetBranchAddress("b3cosalphabs", &b3cosalphabs, &b_b3cosalphabs);
   fChain->SetBranchAddress("b3cosalphabserr", &b3cosalphabserr, &b_b3cosalphabserr);
   fChain->SetBranchAddress("b3lsbs", &b3lsbs, &b_b3lsbs);
   fChain->SetBranchAddress("b3lsbserr", &b3lsbserr, &b_b3lsbserr);
   fChain->SetBranchAddress("b3ctau", &b3ctau, &b_b3ctau);
   fChain->SetBranchAddress("b3ctauerr", &b3ctauerr, &b_b3ctauerr);
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
