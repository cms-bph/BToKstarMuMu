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
  h_mumumass,
  h_kstarpmass,
  h_bpmass,
  h_bpmass_res,
  h_bpmass_res_zoom_5_6,
  h_bpmass_res_zoom_5p3_5p8,
  h_bpmass_nonres,
  h_bpmass_nonres_zoom_5_6,
  kHistNameSize
};

// Global hist args

HistArgs hist_args[kHistNameSize] = {
  // name, title, n_bins, x_min, x_max  
  {"h_mumumass", "#mu^{+}#mu^{-} invariant mass; M(#mu^{+}#mu^{-}) [GeV]", 100, 2, 4},
  {"h_kstarpmass", "K^{*+} mass; M(K^{*+}) [GeV]", 100, 0.6, 1.2},
  {"h_bpmass", "B^{+} mass; M(B^{+}) [GeV]", 100, 1, 10},
  {"h_bpmass_res", "B^{+} mass with J/#psi and #psi(2S) resonant; M(B^{+}) [GeV]", 100, 1, 10},
  {"h_bpmass_res_zoom_5_6", "B^{+} mass with J/#psi and #psi(2S) resonant; M(B^{+}) [GeV]", 100, 5, 6},
  {"h_bpmass_res_zoom_5p3_5p8", "B^{+} mass with J/#psi and #psi(2S) resonant; M(B^{+}) [GeV]", 100, 5.3, 5.8},
  {"h_bpmass_nonres", "B^{+} mass non-resonant; M(B^{+}) [GeV]", 100, 1, 10},
  {"h_bpmass_nonres_zoom_5_6", "B^{+} mass non-resonant; M(B^{+}) [GeV]", 100, 5, 6},
};

// Define histograms 
TH1F *BuToKstarMuMuFigures[kHistNameSize];

TChain *ch; 
TCanvas *c;

double Mumumass = 0; 
double Kstarmass = 0; 

double Bmass = 0; 
double Bpt = 0; 
int Bchg = 0; 
double Bvtxcl = 0; 
double Blxysig = 0; 
double Bcosalphabs = 0; 
double Bctau = 0; 

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



bool sel_bmass_res(TString label, double mumumass){
  if (label == "Run2011v10.3" or 
      label == "Run2011v10.4" or
      label == "Run2011v10.5" or
      label == "Run2011v11.1"
      ){
    double jpsimass_min = JPSI_MASS - 5*0.03 ; // 2.94692
    double jpsimass_max = JPSI_MASS + 5*0.03 ; // 3.24692
    
    double psi2smass_min = PSI2S_MASS - 5*0.03 ; // 3.53609 
    double psi2smass_max = PSI2S_MASS + 5*0.03 ; // 3.83609
    
    if ( (mumumass > jpsimass_min and mumumass < jpsimass_max)
	 or ( mumumass > psi2smass_min and mumumass < psi2smass_max)) 
      return true ;
  }
  
  else{
    cerr << "\n !!! Unkonwn label: " << label << endl; 
    gSystem->Exit(0);
  }

  return false; 
}
  


void summary(TString label, TString infile, TString outfile){
  ch = new TChain("tree"); 
  ch->Add(infile.Data()); 
  
  ch->SetBranchAddress("Mumumass", &Mumumass);
  ch->SetBranchAddress("Bchg", &Bchg);
  ch->SetBranchAddress("Kstarmass", &Kstarmass);
  ch->SetBranchAddress("Bmass", &Bmass);

  // create histograms 
  for(int i=0; i<kHistNameSize; i++) {
    BuToKstarMuMuFigures[i] = new TH1F(hist_args[i].name, hist_args[i].title,
				       hist_args[i].n_bins,
				       hist_args[i].x_min, hist_args[i].x_max);
  }

  // fill histograms
  Int_t nentries = (Int_t)ch->GetEntries();
  
  if (nentries == 0) {
    cerr << "No entries found!" << endl; 
    gSystem->Exit(0);
  }

  for (Int_t i=0;i<nentries;i++) {
    ch->GetEntry(i);
    BuToKstarMuMuFigures[h_mumumass]->Fill(Mumumass); 
  
    if (Bchg > 0) {
      BuToKstarMuMuFigures[h_kstarpmass]->Fill(Kstarmass); 
      BuToKstarMuMuFigures[h_bpmass]->Fill(Bmass); 
      if (sel_bmass_res(label, Mumumass)) {
	BuToKstarMuMuFigures[h_bpmass_res]->Fill(Bmass); 

	if (Bmass > 5 and Bmass < 6) {
	  BuToKstarMuMuFigures[h_bpmass_res_zoom_5_6]->Fill(Bmass); 
	}

	if (Bmass > 5.3 and Bmass < 5.8) {
	  BuToKstarMuMuFigures[h_bpmass_res_zoom_5p3_5p8]->Fill(Bmass); 
	}
	
      } else {
	BuToKstarMuMuFigures[h_bpmass_nonres]->Fill(Bmass); 
	if (Bmass > 5 and Bmass < 6) {
	  BuToKstarMuMuFigures[h_bpmass_nonres_zoom_5_6]->Fill(Bmass); 
	}
	
      }
      
    }
  }
  
  // save figures
  set_root_style(); 

  c = new TCanvas("c","c", 640, 640); 
  c->UseCurrentStyle() ;
  TFile *f = new TFile(outfile, "recreate");

  TString pdffile = outfile; 
  pdffile.ReplaceAll(".root", ".pdf"); 
  c->Print(Form("%s[", pdffile.Data()));

  for(int i = 0; i < kHistNameSize; i++) {
    BuToKstarMuMuFigures[i]->Draw();
    c->Print(pdffile.Data());
    BuToKstarMuMuFigures[i]->Write();
    BuToKstarMuMuFigures[i]->Delete();
  }
  c->Print(Form("%s]", pdffile.Data()));
  delete c; 
  f->Close(); 
  
}



//----------------------------------------------------------
TChain* add_chain(TString datatype, TString label){
  TChain *globChain = new TChain("tree");
  
  TString base = "/Users/xshi/work/cms/afb/dat/sel"; 

  TString fNameList = Form("%s/db/%s/%s/rootfiles.list", base.Data(), 
			   datatype.Data(), label.Data());

  cout << ">>> Load Chain from file: " << fNameList << endl;

  // ifstream fList((char*)fNameList);
  ifstream fList(fNameList.Data());
  // ifstream fList; 
  // fList.open(fNameList.Data());
  if (!fList)
    {
      cerr << "!!! Can't open file " << fNameList << endl;
      return NULL;
    }


  char lineFromFile[255];
  while(fList.getline(lineFromFile, 250))
    {
      TString fileName = lineFromFile;

      fileName = Form("%s/%s/%s", base.Data(), datatype.Data(), fileName.Data());

      if(globChain->Add(fileName))
        cout << ">> File '" << fileName << "' has been loaded" << endl;
      else
        cout << ">> Can't load file '" << fileName << "'" << endl;
    }
  
  cout << ">> Total number of entries: " << 
    globChain->GetEntries() << endl;
  fList.close();
  return globChain; 

}
//----------------------------------------------------------



vector<TString> add_datafiles(TString label){
  vector<TString> datafiles(0); 
  
  if (label == "run2011v0/cut0") {
    TString base = "/Users/xshi/work/cms/afb/dat/sel/data/"; 
    datafiles.push_back(Form("%s%s/cut0/SingleBuToKstarMuMu_%s.root", base.Data(),
			     "Run2011A_May10ReReco_v1_run2011v0_2", "1")); 
    
    datafiles.push_back(Form("%s%s/cut0/SingleBuToKstarMuMu_%s.root", base.Data(),
			     "Run2011A_May10ReReco_v1_run2011v0_2", "2")); 

    datafiles.push_back(Form("%s%s/cut0/SingleBuToKstarMuMu_%s.root", base.Data(),
			     "Run2011A_PromptReco_v4_run2011v0_1", "1")); 

    datafiles.push_back(Form("%s%s/cut0/SingleBuToKstarMuMu_%s.root", base.Data(),
			     "Run2011A_PromptReco_v4_run2011v0_1", "2")); 

    datafiles.push_back(Form("%s%s/cut0/SingleBuToKstarMuMu_%s.root", base.Data(),
			     "Run2011A_PromptReco_v4_run2011v0_1", "3")); 
    
    datafiles.push_back(Form("%s%s/cut0/SingleBuToKstarMuMu_%s.root", base.Data(),
			     "Run2011A_PromptReco_v5_run2011v0", "1")); 
  
    datafiles.push_back(Form("%s%s/cut0/SingleBuToKstarMuMu_%s.root",  base.Data(),
			     "Run2011A_PromptReco_v5_run2011v0", "1")); 

    datafiles.push_back(Form("%s%s/cut0/SingleBuToKstarMuMu_%s.root",  base.Data(),
			     "Run2011A_PromptReco_v6_run2011v0", "1")); 

    datafiles.push_back(Form("%s%s/cut0/SingleBuToKstarMuMu_%s.root", base.Data(),
			     "Run2011B_PromptReco_v1_run2011v0", "1")); 

    datafiles.push_back(Form("%s%s/cut0/SingleBuToKstarMuMu_%s.root", base.Data(),
			     "Run2011B_PromptReco_v1_run2011v0", "2")); 

    datafiles.push_back(Form("%s%s/cut0/SingleBuToKstarMuMu_%s.root", base.Data(),
			     "Run2011B_PromptReco_v1_run2011v0", "3")); 
    
  }
    
  return datafiles; 
  
}


void bpmass(TString label, TString outfile){
  // cout << "Making bpmass figure ..." << endl; 

  // plot the BuToKstarJPsi part

  // ch = new TChain("tree"); 

  // TString infile = "/Users/xshi/work/cms/afb/dat/sel/mc/BuToKstarJPsi_7TeV_5E5_v1_run2011v0_2/cut0/SingleBuToKstarMuMu.root"; 
  // ch->Add(infile.Data()); 

  // vector<TString> datafiles = get_datafiles(label); 
  // for (vector<TString>::iterator it = datafiles.begin(); 
  //      it != datafiles.end(); ++it) {
  //   ch->Add(*it); 
  // }

  
  ch = add_chain("data", label); 

  
  ch->SetBranchAddress("Bchg", &Bchg);
  ch->SetBranchAddress("Bmass", &Bmass);

  // fill histograms
  Int_t nentries = (Int_t)ch->GetEntries();
  
  if (nentries == 0) {
    cerr << "No entries found!" << endl; 
    gSystem->Exit(0);
  }
  
  TH1F* h_bpmass = new TH1F("h_bpmass", "B^{+} mass; M(B^{+}) [GeV]", 100, 5, 5.6); 
  
  for (Int_t i=0;i<nentries;i++) {
    ch->GetEntry(i);
    
    if (Bchg > 0) 
      h_bpmass->Fill(Bmass); 
  }
  
  // save figures
  set_root_style(); 

  c = new TCanvas("c","c", 400, 400); 
  c->UseCurrentStyle() ;
  // TString pdffile = "test.pdf"; 
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

