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

TChain* add_chain(TString datatype, TString label, TString cut, int verbose=0){
  TChain *globChain = new TChain("tree");
  // TString base = Form("%s/sel", getenv("dat")); 
  TString fNameList = Form("../data/sel_%s_%s_%s_rootfiles.txt",
			   datatype.Data(), label.Data(), cut.Data());
  // if (verbose >0 ) 
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
      fileName = Form("%s/sel/%s/%s", getenv("dat"), 
		      datatype.Data(), fileName.Data());
      
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


void bmass(TString datatype, TString label, TString cut, TString outfile){
  // plot the BuToKstarJPsi part
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
    h_bmass->Fill(Bmass); 
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

  TString func     = argv[1]; 
  TString datatype = argv[2]; 
  TString label    = argv[3]; 
  TString cut      = argv[4]; 
  TString outfile  = argv[5]; 
  
  if (func == "bmass") 
    bmass(datatype, label, cut, outfile); 
  else 
    cerr << "No function available for: " << func.Data() << endl; 

  gSystem->Exit(0);

  return 0 ;
}

#endif

