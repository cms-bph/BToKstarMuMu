// -----------------------------------------------
//       Draw pre-selection variables 
//       Author: Xin Shi <Xin.Shi@cern.ch> 
//       Created:   [2013-07-30 Tue 10:52] 
// -----------------------------------------------
#include <iostream>
#include <TSystem.h>
#include <TH1.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TStyle.h>
#include <TROOT.h>

using namespace std; 

// Global Constants
const double JPSI_MASS = 3.09692; // GeV
const double PSI2S_MASS = 3.68609; // GeV
const double KSTAR_MASS = 0.89166; // GeV
const double KSTAR_WIDTH = 0.0508; // GeV 
const double MUON_MASS = 0.10565837;
const double KAON_MASS = 0.493677;
const double PION_MASS = 0.13957018;
const double KSHORT_MASS = 0.497614;

// Structures 
struct HistArgs{
    char name[128];
    char title[128];
    int n_bins;
    double x_min;
    double x_max;
};

enum HistName{
  h_truemupt,
  kHistNameSize
};

// Global hist args

HistArgs hist_args[kHistNameSize] = {
  // name, title, n_bins, x_min, x_max  
  {"h_truemupt", "True muon pT; pT [GeV]", 100, 0, 30}
};

// Define histograms 
TH1F *histos[kHistNameSize];

TChain *ch; 
TCanvas *c;

void set_root_style(int stat=1110, int grid=0){
  
  gROOT->Reset();

  gStyle->SetTitleFillColor(0) ; 
  gStyle->SetTitleBorderSize(0); 
    
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasDefX(0); 
  gStyle->SetCanvasDefY(0); 
  gStyle->SetFrameBorderMode(0); 
  gStyle->SetFrameBorderSize(1); 
  gStyle->SetFrameFillColor(0); 
  gStyle->SetFrameFillStyle(0); 
  gStyle->SetFrameLineColor(1); 
  gStyle->SetFrameLineStyle(1); 
  gStyle->SetFrameLineWidth(1); 

  gStyle->SetPadBorderMode(0);  
  gStyle->SetPadColor(0);  
  gStyle->SetPadTickX(1); 
  gStyle->SetPadTickY(1); 
  gStyle->SetPadGridX(grid); 
  gStyle->SetPadGridY(grid); 

  gStyle->SetOptStat(stat); 
  gStyle->SetStatColor(0); 
  gStyle->SetStatBorderSize(1); 
}




void summary(TString infile, TString outfile){
  // create histograms 
  for(int i=0; i<kHistNameSize; i++) {
    histos[i] = new TH1F(hist_args[i].name, hist_args[i].title,
			 hist_args[i].n_bins, hist_args[i].x_min, 
			 hist_args[i].x_max);
  }

  TFile *fi = TFile::Open(infile); 
  TH1F *h = (TH1F*) fi->Get("h_mupt"); 

  TTree *t = (TTree*) fi->Get("tree"); 
  vector<double>  *bmass = 0;  // must init with 0! 
  vector<bool>    *istruebu = 0;
  vector<double>  *muppt = 0; 
  vector<double>  *mumpt = 0; 
  t->SetBranchAddress("bmass", &bmass); 
  t->SetBranchAddress("istruebu", &istruebu); 
  t->SetBranchAddress("muppt", &muppt); 
  t->SetBranchAddress("mumpt", &mumpt); 

  Int_t nentries = (Int_t)t->GetEntries();
  // Int_t nentries = 1000;  

  for (Int_t i=0;i<nentries;i++) {
    t->GetEntry(i);
    if (bmass->size() < 1) continue; 
    
    for (vector<int>::size_type i = 0; i < bmass->size(); i++) {
      if (istruebu->at(i)) histos[h_truemupt]->Fill(muppt->at(i)); 
      if (istruebu->at(i)) histos[h_truemupt]->Fill(mumpt->at(i)); 
    }
  }

  // save figures
  set_root_style(); 

  c = new TCanvas("c","c", 800, 400); 
  c->UseCurrentStyle() ;
  TFile *fo = new TFile(outfile, "recreate");

  TString pdffile = outfile; 
  pdffile.ReplaceAll(".root", ".pdf"); 

  c->Print(Form("%s[", pdffile.Data()));

  c->Divide(2, 1); 
  // p1 
  c->cd(1); 
  h->Draw(); 

  c->cd(2); 
  histos[h_truemupt]->Draw(); 
  c->Print(pdffile.Data());

  h->Write(); 
  h->Delete(); 
  


  c->Print(Form("%s]", pdffile.Data()));
  delete c; 
  fo->Close(); 
 
}

  

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
  cerr << "Usage: Singlehistos label infile outfile \n"
       << "Options: \n" 
       << " -h \t\tPrint this info\n"
       << endl; 
}

int main(int argc, char** argv) {
  if ( (argc < 3) or option_exists(argv, argv+argc, "-h") ){
    print_usage() ;  
    return -1; 
  }

  TString infile = argv[1]; 
  TString outfile = argv[2]; 
  
  Printf("input file: '%s'", infile.Data());
  Printf("output file: '%s'", outfile.Data());
  
  summary(infile, outfile); 
  
  gSystem->Exit(0);

  return 0 ;
}

#endif

