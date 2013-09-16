// -----------------------------------------------
//       Author: Xin Shi <Xin.Shi@cern.ch> 
//       Created:   [2013-04-04 Thu 08:42] 
// -----------------------------------------------
#include <iostream>
#include <TSystem.h>
#include <TH1.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TMath.h>
#include "tools.h" 

using namespace std; 

const double KSTARPLUS_MASS = 0.89166; // GeV
const double BPLUS_MASS = 5.279 ; // GeV 
const double JPSI_MASS = 3.097 ; // GeV  ±0.011

void bmass(TString datatype, TString label, TString cut, TString outfile){
  // plot the BuToKstarJPsi part
  TH1F* h_bmass = new TH1F("h_bmass", "B^{+/-} mass; M(B^{+/-}) [GeV/c^{2}]",
			   100, 5, 5.6); 
  TChain* ch = add_chain(datatype, label, cut); 
  if (ch == NULL) gSystem->Exit(0);

  double Bmass = 0; 
  ch->SetBranchAddress("Bmass", &Bmass);

  double Mumumass = 0; 
  ch->SetBranchAddress("Mumumass", &Mumumass); 
  double Mumumasserr = 0; 
  ch->SetBranchAddress("Mumumasserr", &Mumumasserr); 

  Int_t nentries = (Int_t)ch->GetEntries();
  if (nentries == 0) {
    cerr << "No entries found!" << endl; 
    gSystem->Exit(0);
  }
  
  for (Int_t i=0;i<nentries;i++) {
    ch->GetEntry(i);
    if ( cut.Contains("/jpsi")) {
      if ( Mumumass > ( JPSI_MASS - 5*Mumumasserr) && 
	   Mumumass < ( JPSI_MASS + 3*Mumumasserr) ) {
	h_bmass->Fill(Bmass); 
      }
      
    } else {
      h_bmass->Fill(Bmass); 
    }
  }

  TCanvas* c = new TCanvas("c","c", 400, 400); 
  set_root_style(); 
  c->UseCurrentStyle() ;

  h_bmass->SetMinimum(0); 
  h_bmass->Draw();
  c->Print(outfile.Data());
  h_bmass->Delete(); 
  delete c; 
}

void bmass_snb(TString datatype, TString label, TString cut, TString outfile){
  float sigma_mb = 0.05 ; // GeV 


  // Signal part with BuToKstarJPsi
  datatype = "BuToKstarJPsi"; 

  TH1F* h_bmass = new TH1F("h_bmass", "B^{+/-} mass; M(B^{+/-}) [GeV/c^{2}]",
			   100, 5, 5.6); 

  TChain* ch = add_chain(datatype, label, cut); 
  if (ch == NULL) gSystem->Exit(0);

  double Bmass = 0; 
  ch->SetBranchAddress("Bmass", &Bmass);

  Int_t nentries = (Int_t)ch->GetEntries();
  if (nentries == 0) {
    cerr << "No entries found!" << endl; 
    gSystem->Exit(0);
  }
  
  for (Int_t i=0;i<nentries;i++) {
    ch->GetEntry(i);
    if (fabs( Bmass - BPLUS_MASS) < 2.5*sigma_mb )
      h_bmass->Fill(Bmass); 
  }

  // double scaleFactor = 10; 

  double scaleFactor = calc_scale_factor(datatype, "7TeV"); 
  cout << "scalefactor = " << scaleFactor << endl; 
  h_bmass->Scale(scaleFactor); 

  // Backgroud from data sideband
  datatype = "data" ;
 
  TH1F* h_bmass_bkg = new TH1F("h_bmass_bkg", "B^{+/-} mass; M(B^{+/-}) [GeV/c^{2}]",
			       100, 5, 5.6); 
  
  TChain* ch_bkg = add_chain(datatype, label, cut); 
  if (ch_bkg == NULL) gSystem->Exit(0);

  Bmass = 0; 
  ch_bkg->SetBranchAddress("Bmass", &Bmass);

  nentries = (Int_t)ch_bkg->GetEntries();
  if (nentries == 0) {
    cerr << "No entries found!" << endl; 
    gSystem->Exit(0);
  }
  
  for (Int_t i=0;i<nentries;i++) {
    ch_bkg->GetEntry(i);
    
    if (fabs( Bmass - BPLUS_MASS) > 2.5*sigma_mb &&  
	fabs( Bmass - BPLUS_MASS) < 5.0*sigma_mb)
      h_bmass_bkg->Fill(Bmass); 
  }

  TCanvas* c = new TCanvas("c","c", 400, 400); 
  set_root_style(); 
  c->UseCurrentStyle() ;

  h_bmass->SetMinimum(0); 

  h_bmass->Draw();
  h_bmass_bkg->Draw("same"); 
  c->Print(outfile.Data());
  h_bmass->Delete(); 
  delete c; 
}



int main(int argc, char** argv) {
  TString func     = argv[1]; 
  TString datatype = argv[2]; 
  TString label    = argv[3]; 
  TString cut      = argv[4]; 
  TString outfile  = argv[5]; 
  
  if (func == "bmass") 
    bmass(datatype, label, cut, outfile); 
  else if (func == "bmass_snb") 
    bmass_snb(datatype, label, cut, outfile); 
  else 
    cerr << "No function available for: " << func.Data() << endl; 

  gSystem->Exit(0);
  return 0 ;
}


