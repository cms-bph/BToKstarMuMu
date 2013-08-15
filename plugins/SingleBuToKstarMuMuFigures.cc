// -----------------------------------------------
//       Author: Xin Shi <Xin.Shi@cern.ch> 
//       Created:   [2013-04-04 Thu 08:42] 
// -----------------------------------------------
#include <iostream>
#include <fstream>
#include <TSystem.h>
#include <TH1.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TStyle.h>
#include <TROOT.h>

using namespace std; 

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

TChain* add_chain(TString datatype, TString label, int verbose=0){
  TChain *globChain = new TChain("tree");
  TString base = Form("%s/sel", getenv("dat")); 
  TString fNameList = Form("%s/db/%s/%s/rootfiles.list", base.Data(), 
			   datatype.Data(), label.Data());
  if (verbose >0 ) 
    cout << ">>> Load Chain from file: " << fNameList << endl;

  ifstream fList(fNameList.Data());
  if (!fList)
    { cerr << "!!! Can't open file " << fNameList << endl;
      return NULL;
    }

  char lineFromFile[255];
  while(fList.getline(lineFromFile, 250))
    {
      TString fileName = lineFromFile;
      fileName = Form("%s/%s/%s", base.Data(), datatype.Data(), fileName.Data());
      
      if(globChain->Add(fileName)){
	if (verbose >0 ) 
	  cout << ">> File '" << fileName << "' has been loaded" << endl;
      }
      else
        cerr << ">> Can't load file '" << fileName << "'" << endl;
    }
  
  if (verbose >0 ) 
    cout << ">> Total number of entries: " << 
      globChain->GetEntries() << endl;
  
  fList.close();
  return globChain; 
}


void bpmass(TString label, TString outfile){
  // plot the BuToKstarJPsi part
  TH1F* h_bpmass = new TH1F("h_bpmass", "B^{+} mass; M(B^{+}) [GeV]", 100, 5, 5.6); 

  TChain* ch = add_chain("data", label); 

  double Bmass = 0; 
  int    Bchg = 0; 

  ch->SetBranchAddress("Bchg", &Bchg);
  ch->SetBranchAddress("Bmass", &Bmass);

  // fill histograms
  Int_t nentries = (Int_t)ch->GetEntries();
  
  if (nentries == 0) {
    cerr << "No entries found!" << endl; 
    gSystem->Exit(0);
  }
  
  for (Int_t i=0;i<nentries;i++) {
    ch->GetEntry(i);
    if (Bchg > 0) 
      h_bpmass->Fill(Bmass); 
  }
  
  // save figures

  TCanvas* c = new TCanvas("c","c", 400, 400); 
  set_root_style(); 
  c->UseCurrentStyle() ;

  h_bpmass->Draw();
  c->Print(outfile.Data());
  h_bpmass->Delete(); 
  delete c; 
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
  cerr << "Usage: SingleBuToKstarMuMuFigures label infile outfile \n"
       << "Options: \n" 
       << " -h \t\tPrint this info\n"
       << endl; 
}

int main(int argc, char** argv) {
  if ( (argc < 4) or option_exists(argv, argv+argc, "-h") ){
    print_usage() ;  
    return -1; 
  }

  TString func = argv[1]; 
  TString label = argv[2]; 
  TString outfile = argv[3]; 
  
  if (func == "bpmass") 
    bpmass(label, outfile); 
  else 
    cerr << "No function available for: " << func.Data() << endl; 

  gSystem->Exit(0);

  return 0 ;
}

#endif

