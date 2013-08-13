// -----------------------------------------------
//       Author: Xin Shi <Xin.Shi@cern.ch> 
//       Created:   <2013-02-26 Tue 12:08> 
// -----------------------------------------------
#define SingleBuToKstarMuMuSelector_cxx

#include <iostream>
#include <sstream>
#include "SingleBuToKstarMuMuSelector.h"
#include <TH2.h>
#include <TStyle.h>
#include <TProof.h>
#include <TLorentzVector.h>

// Global Constants
const double KSTAR_MASS = 0.89166; // GeV
const double KSTAR_WIDTH = 0.0508; // GeV 
const double MUON_MASS = 0.10565837;
const double KAON_MASS = 0.493677;
const double PION_MASS = 0.13957018;
const double KSHORT_MASS = 0.497614;

// user defined variables
TDatime t_begin_ , t_now_ ;
int n_processed_, n_selected_; 

TTree *tree_; 

// Branch variables for new tree
int    Nb = 0;
double Mumumass = 0; 
double Kstarmass = 0; 
double Trkpt = 0; 
double Trkdcasigbs = 0; 

double Bmass = 0; 
double Bpt = 0; 
int    Bchg = 0; 
double Bvtxcl = 0; 
double Blxysig = 0; 
double Bcosalphabs = 0; 
double Bctau = 0; 


void str_replace(std::string& str, const std::string& oldStr,
		 const std::string& newStr)
{
  size_t pos = 0;
  while((pos = str.find(oldStr, pos)) != std::string::npos)
  {
    str.replace(pos, oldStr.length(), newStr);
    pos += newStr.length();
  }
}


string get_option_value(string option, string name){
  vector<string> args;
  istringstream f(option);
  string s;    
  while (getline(f, s, ';')) {
    args.push_back(s);
  }
  
  string value; 
  for(vector<string>::iterator it = args.begin(); it != args.end(); ++it) {
    value = *it; 
    unsigned found = value.find(name);
    if (found == 0) {
      str_replace(value, name+"=", ""); 
      break; 
    }
  }
  return value; 
}


void SingleBuToKstarMuMuSelector::Begin(TTree * /*tree*/){
  t_begin_.Set(); 
  printf("\n ---------- Begin Job ---------- \n");
  t_begin_.Print();

  n_processed_ = 0;
  n_selected_ = 0;
}

void SingleBuToKstarMuMuSelector::SlaveBegin(TTree * /*tree*/){
  // printf("\n ---------- Begin of Slave Job ---------- \n");
  TString option = GetOption();
  tree_ = new TTree("tree", "tree"); 
  
  tree_->Branch("Mumumass", &Mumumass, "Mumumass/D");
  tree_->Branch("Kstarmass", &Kstarmass, "Kstarmass/D");
  tree_->Branch("Trkpt", &Trkpt, "Trkpt/D");
  tree_->Branch("Trkdcasigbs", &Trkdcasigbs, "Trkdcasigbs/D");

  tree_->Branch("Bmass", &Bmass, "Bmass/D");
  tree_->Branch("Bpt", &Bpt, "Bpt/D");
  tree_->Branch("Bchg", &Bchg, "Bchg/I");

  tree_->Branch("Bvtxcl", &Bvtxcl, "Bvtxcl/D");
  tree_->Branch("Blxysig", &Blxysig, "Blxysig/D");
  tree_->Branch("Bcosalphabs", &Bcosalphabs, "Bcosalphabs/D");
  tree_->Branch("Bctau", &Bctau, "Bctau/D");
  
  fOutput->AddAll(gDirectory->GetList()); 
}


Bool_t SingleBuToKstarMuMuSelector::Process(Long64_t entry){

  string option = GetOption();
  string datatype = get_option_value(option, "datatype"); 
  string cut = get_option_value(option, "cut"); 
  
  GetEntry(entry); 
  n_processed_ += 1; 
  Nb = nb; 

  int i = SelectB(cut); 
  if ( i != -1) {
    
    if ( datatype == "data" || 
	 ( datatype == *decname && istruebu->at(i) ) ) {
      // cout << "true Bu? :" << istruebu->at(i) 
      // 	   << "  decname = " << *decname << endl; 
      n_selected_ += 1; 
      
      SaveEvent(i);     
      // SaveMuMu(i);
      tree_->Fill();	   
    }
  }

  return kTRUE;
}


void SingleBuToKstarMuMuSelector::SlaveTerminate(){
  // printf ( "\n ---------- End of Slave Job ---------- \n" ) ;
  // t_now_.Set() ; 
  // t_now_.Print() ;
  // printf(" processed: %i \n selected: %i \n duration: %i sec \n rate: %g evts/sec\n",
  // 	 n_processed_, n_selected_, 
  // 	 t_now_.Convert() - t_begin_.Convert(), 
  // 	 float(n_processed_)/(t_now_.Convert()-t_begin_.Convert()) );
}

void SingleBuToKstarMuMuSelector::Terminate(){
  string option = GetOption();
  TString outfile = get_option_value(option, "outfile"); 
  
  TFile file(outfile.Data(), "recreate"); 
  fOutput->Write();
  
  t_now_.Set(); 
  printf(" \n ---------- End Job ---------- \n" ) ;
  t_now_.Print();  
  printf(" processed: %i \n selected: %i \n duration: %i sec \n rate: %g evts/sec\n",
	 n_processed_, n_selected_, 
	 t_now_.Convert() - t_begin_.Convert(), 
	 float(n_processed_)/(t_now_.Convert()-t_begin_.Convert()) );
}


int SingleBuToKstarMuMuSelector::SelectB(string cut){

  int best_idx = -1; 
  double best_bvtxcl = 0.0; 

  if (cut == "cut0") {
    for (int i = 0; i < nb; i++) {

      if (bvtxcl->at(i) > best_bvtxcl) {
	best_bvtxcl = bvtxcl->at(i); 
	best_idx = i; 
      }
    }
  }




  // if ( ! HasGoodDimuon() ) return -1; 

  // for (vector<int>::size_type i = 0; i < bmass->size(); i++) {
  // for (int i = 0; i < Nb; i++) {
    // for (int i = 0; i < nb; i++) {

    
    // if ( ! HasGoodDimuon(i) ) continue;  
    
    // cout << "passed dimuon cut" << endl; 

    // if (bvtxcl->at(i) < 0.09) continue; 
    // double blxysig = blsbs->at(i)/blsbserr->at(i); 
    // if (blxysig < 12 ) continue; 
    // if (bcosalphabs->at(i) < 0.99) continue; 
    
    // if ( label == "Run2011v10.1" ) 
    //   if (bctau->at(i) < 0.03) continue; 

    // // Kstarmass = GetKstarMass(i);
    // Kstarmass = kstarmass->at(i);
    // float kstar_mass_delta; 
    // // if ( label == "Run2011v11.1" ) 
    // kstar_mass_delta = 0.06; 
      
    // if (Kstarmass < (KSTAR_MASS - kstar_mass_delta) 
    // 	|| Kstarmass > (KSTAR_MASS + kstar_mass_delta))
    //   continue; 
    // cout << "n processed: " << n_processed_ << ", bvtxcl[" << i << "] ="
    // 	   << bvtxcl->at(i)   << endl; 
    

    // if (bvtxcl->at(i) > best_bvtxcl) {
    //   best_bvtxcl = bvtxcl->at(i); 
    //   best_idx = i; 
    // }

  
  
  return best_idx;
}

bool SingleBuToKstarMuMuSelector::HasGoodDimuon(int i){
  //  for (vector<int>::size_type i = 0; i < mumpx->size(); i++) {
  if  ( 
       // soft muon 
       mumisgoodmuon->at(i)
       && mupisgoodmuon->at(i) 
       && mumntrkhits->at(i) > 10 
       && mupntrkhits->at(i) > 10 
       && mumnpixlayers->at(i) > 1
       && mupnpixlayers->at(i) > 1
       && mumnormchi2->at(i) < 1.8 
       && mupnormchi2->at(i) < 1.8 
       && mumdxyvtx->at(i) < 3
       && mupdxyvtx->at(i) < 3
       && mumdzvtx->at(i) < 30 
       && mupdzvtx->at(i) < 30 
       
	) return true; 
  // }
  return false; 
}


// void SingleBuToKstarMuMuSelector::SaveMuMu(int i){
//   TLorentzVector mup, mum, dimu; 
//   mup.SetXYZM(muppx->at(i), muppy->at(i), muppz->at(i), MUON_MASS); 
//   mum.SetXYZM(mumpx->at(i), mumpy->at(i), mumpz->at(i), MUON_MASS); 
//   dimu = mup + mum; 
//   Mumumass = dimu.M(); 
// }


// void SingleBuToKstarMuMuSelector::SaveKstar(int i){
//   // TLorentzVector ks, pi, kstar; 
//   // ks.SetXYZM(bkspx->at(i), bkspy->at(i), bkspz->at(i), KSHORT_MASS); 
//   // pi.SetXYZM(bpi1px->at(i), bpi1py->at(i), bpi1pz->at(i), PION_MASS); 
//   // kstar = ks + pi; 
//   Kstarmass = GetKstarMass(i);
// }


void SingleBuToKstarMuMuSelector::SaveEvent(int i){
  Bmass = bmass->at(i); 
  Bchg = bchg->at(i); 
  Bvtxcl = bvtxcl->at(i); 
  Blxysig = (blsbs->at(i)/blsbserr->at(i)); 
  Bcosalphabs = bcosalphabs->at(i); 
  Bctau = bctau->at(i); 

  TLorentzVector b; 
  b.SetXYZM(bpx->at(i), bpy->at(i), bpz->at(i), bmass->at(i)); 
  Bpt = b.Pt(); 

  Mumumass = mumumass->at(i); 
  Kstarmass = kstarmass->at(i); 
  Trkpt = trkpt->at(i); 
  Trkdcasigbs = fabs( trkdcabs->at(i)/trkdcabserr->at(i) ); 
  
}

// double SingleBuToKstarMuMuSelector::GetKstarMass(int i){
//   TLorentzVector ks, pi, kstar; 
//   ks.SetXYZM(bkspx->at(i), bkspy->at(i), bkspz->at(i), KSHORT_MASS); 
//   pi.SetXYZM(bpi1px->at(i), bpi1py->at(i), bpi1pz->at(i), PION_MASS); 
//   kstar = ks + pi; 
//   return kstar.M();
// }


#ifndef __CINT__ 
#include <algorithm>

char* get_option(char ** begin, char ** end, const std::string & option){
  char ** itr = std::find(begin, end, option);
  if (itr != end && ++itr != end)  return *itr;
  return 0;
}


bool option_exists(char** begin, char** end, const std::string& option){
  return std::find(begin, end, option) != end;
}

void print_usage(){
  cerr << "Usage: SingleBuToKstarMuMuSelector infile outfile [-n] [-j] [-h]\n"
       << "Options: \n" 
       << " -h \t\tPrint this info\n"
       << " -n \t\tNumber of entries\n" 
       << " -j \t\tNumber of workers\n" 
       << endl; 
}



int main(int argc, char** argv) {
 if ( (argc < 2) or option_exists(argv, argv+argc, "-h") ){
    print_usage() ;  
    return -1; 
  }

  TString datatype = argv[1]; 
  TString cut = argv[2]; 
  TString infile = argv[3]; 
  TString outfile = argv[4]; 
  
  Printf("datatype: '%s'", datatype.Data());
  Printf("cut: '%s'", cut.Data());
  Printf("input file: '%s'", infile.Data());
  Printf("output file: '%s'", outfile.Data());
 
  TString option; 
  option.Form("datatype=%s;cut=%s;outfile=%s", datatype.Data(), 
	      cut.Data(), outfile.Data()); 

  TChain *ch = new TChain("tree"); 
  ch->Add(infile.Data()); 
   
  char * j = get_option(argv, argv + argc, "-j");
  if (j) {
    TProof::Open(Form("workers=%s", j));
    ch->SetProof(); 
  }

  Long64_t nentries = 1000000000; 
  char * n = get_option(argv, argv+argc, "-n");  
  if (n) nentries = atoi(n);
  
  ch->Process("SingleBuToKstarMuMuSelector.cc+", option, nentries); 

  gSystem->Exit(0);

  return 0 ;
}

#endif

