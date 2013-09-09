// -----------------------------------------------
//       Author: Xin Shi <Xin.Shi@cern.ch> 
//       Created:   [2013-08-15 Thu 15:56] 
// -----------------------------------------------
#include <algorithm>
#include <iostream>
#include <fstream>
#include <TStyle.h>
#include <TROOT.h>
#include <TFile.h>
#include <TChain.h>
#include "tools.h"

using namespace std; 

void set_root_style(int stat, int grid){
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

TChain* add_chain(TString datatype, TString label, TString cut, int verbose){

  TString fNameList = Form("../data/sel_%s_%s_%s_rootfiles.txt",
			   datatype.Data(), label.Data(), cut.Data());
  // if (verbose >0 ) 
  cout << ">>> Load Chain from file: " << fNameList << endl;

  TString fileName; 
  int num_lines = get_number_of_lines(fNameList, datatype, fileName); 

  // cout << "Number of lines: " << num_lines << endl; 
  if (num_lines == 1) {
    // cout << fileName.Data() << endl;
    TFile *f = new TFile(fileName.Data());
    TChain *ch = (TChain*)f->Get("tree");
    return ch; 
  }

  
  TChain *globChain = new TChain("tree");
  ifstream fList(fNameList.Data());
  if (!fList)
    { cerr << "!!! Can't open file " << fNameList << endl;
      return NULL;
    }

  string lineFromFile;
  while( getline(fList, lineFromFile) )
    {
      if (lineFromFile.empty()) continue; 

      // TString fileName = lineFromFile;
      fileName = lineFromFile;
      
      if (datatype == "data") {
	fileName = Form("%s/sel/data/%s", getenv("dat"), fileName.Data());
      } else {
	fileName = Form("%s/sel/mc/%s", getenv("dat"), fileName.Data());
      }

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

char* get_option(char ** begin, char ** end, const std::string & option){
  char ** itr = std::find(begin, end, option);
  if (itr != end && ++itr != end)  return *itr;
  return 0;
}


bool option_exists(char** begin, char** end, const std::string& option){
  return std::find(begin, end, option) != end;
}


int get_number_of_lines(TString fNameList, TString datatype, 
			TString &fileName){
  string lineFromFile;
  
  ifstream fList(fNameList.Data());
  
  int n = 0 ; 
  while( getline(fList, lineFromFile) )
    {
      if (lineFromFile.empty()) continue; 

      // TString fileName = lineFromFile;
      fileName = lineFromFile;
      if (datatype == "data") {
	fileName = Form("%s/sel/data/%s", getenv("dat"), fileName.Data());
      } else {
	fileName = Form("%s/sel/mc/%s", getenv("dat"), fileName.Data());
      }
      n++; 
    }
  return n; 
}
