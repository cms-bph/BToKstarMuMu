// -----------------------------------------------
//       Author: Xin Shi <Xin.Shi@cern.ch> 
//       Created:   [2013-08-15 Thu 14:54] 
// -----------------------------------------------
#include <sstream>

#include <TSystem.h>
#include <TH1.h>
#include <TFile.h>
#include <TPad.h> 
#include <TCanvas.h> 
#include <TChain.h> 
#include <TPaveText.h>

#include <RooRealVar.h>
#include <RooGaussian.h>
#include <RooChebychev.h> 
#include <RooAddPdf.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooGenericPdf.h> 

#include "tools.h" 

using namespace std; 
using namespace RooFit ;

void bmass(TString datatype, TString label, TString cut, TString outfile)
{
  // bool test = true; 
  bool test = false; 

  // Importing a  TTree into a RooDataSet with cuts 
  // --------------------------------------------------------------------------
  TChain* ch = add_chain(datatype, label, cut); 
  if (ch == NULL) gSystem->Exit(0);

  RooRealVar x("Bmass", "B^{+/-} mass(GeV/c^{2})", 5.27953-0.28, 5.27953+0.28) ;
  RooDataSet data("data", "data", RooArgSet(x), Import(*ch)) ;
  data.Print();

  // Create model and dataset
  // -----------------------------------------------
  // Gaussian signal 
  RooRealVar mean("mean","mean of gaussians", 5.27, 5.23, 5.32) ;
  RooRealVar sigma("sigma","width of gaussians", 0.0285, 0, 1) ;
  RooGaussian sig("sig","Signal component", x, mean, sigma) ;  

  // Build Chebychev polynomial p.d.f.  
  RooRealVar a0("a0", "constant", 0.5, -1, 1.) ;
  RooRealVar a1("a1", "linear", 0.6, -1, 1) ;
  RooRealVar a2("a2", "quadratic", 0.1, -1, 1) ;
  RooChebychev bkg("bkg", "Background", x, RooArgSet(a0, a1, a2)) ;

  // Construct signal+background PDF
  RooRealVar nsig("nsig", "number of signal events", 4648, 0, 1000000); 
  RooRealVar nbkg("nbkg", "number of background events", 21472, 0, 1000000);
  RooAddPdf  model("model", "g+c", RooArgList(bkg, sig), RooArgList(nbkg, nsig)) ;
  
  // Print structure of composite p.d.f.
  model.Print("t") ;

  // Fit model to data, save fitresult 
  // -----------------------------------------------------------------
  RooFitResult* fitres; 
  if (! test) {
    fitres = model.fitTo(data, Extended(true), Save(true)) ;
    fitres->Print("v"); 
  }
  
  // Plot model 
  // ---------------------------------------------------------
  TString title = "B^{+/-} mass";
  int nbins = 50; 
  RooPlot* xframe = x.frame(Title(title), Bins(nbins));
  data.plotOn(xframe) ;
  model.plotOn(xframe) ;

  // Overlay the background component of model with a dashed line
  model.plotOn(xframe,Components("bkg"), LineStyle(kDashed)) ;

  // Draw the frame on the canvas
  TCanvas* c = new TCanvas("c", "c", 400, 400); 
  set_root_style(); 
  c->UseCurrentStyle() ;

  gPad->SetLeftMargin(0.15) ;
  xframe->GetYaxis()->SetTitleOffset(1.7) ; 
  xframe->Draw();

  TPaveText* paveText = new TPaveText(0.17, 0.80, 0.41, 0.88, "NDC"); 
  paveText->SetBorderSize(0.0);
  paveText->SetFillColor(kWhite);
  paveText->AddText(Form("nsig = %.0f #pm %.0f ", nsig.getVal(), nsig.getError())); 
  paveText->AddText(Form("nbkg = %.0f #pm %.0f ", nbkg.getVal(), nbkg.getError())); 
  paveText->AddText(Form("mean = %.3f #pm %.3f ", mean.getVal(), mean.getError())); 
  paveText->AddText(Form("sigma = %.3f #pm %.3f ", sigma.getVal(), sigma.getError())); 
  paveText->Draw(); 

  TString pdffile = outfile + ".pdf"; 
  c->Print(pdffile); 
  
  // Persist fit result in root file 
  // -------------------------------------------------------------
  TString resfile = outfile + ".root"; 
  TFile resf(resfile, "RECREATE") ;
  gPad->Write("plot"); 
  if (! test) fitres->Write("fitres") ;
  resf.Close() ;

  // In a clean ROOT session retrieve the persisted fit result as follows:
  // RooFitResult* r = gDirectory->Get("fitres") ;
  
  delete c;
  delete paveText; 

}

void fl(TString outfile)
{
  // From fomula (2) in LHCb 2012 PRL108, 181806(2012)
  // integrated over theta_l and phi: 
  // 
  // 1/Gamma * d^2 Gamma/d cos(theta_K) dq^2 = 3/2 * F_L cos^2(theta_K)
  // + 3/4(1-F_L)(1-cos^2theta_K)
  // 

  RooRealVar cosk("cosk", "cos#theta_{K}", 0, 1); 
  RooRealVar fl("fl", "F_{L}l", 0.5);

  RooGenericPdf f("f", "0.5*3*fl*cosk*cosk+0.25*3*(1-fl)*(1-cosk*cosk)", 
		  RooArgSet(cosk,fl));
  RooDataSet* data = f.generate(RooArgSet(cosk), 10000); 
  
  f.fitTo(*data); 

  RooPlot* framecosk = cosk.frame(); 
  data->plotOn(framecosk); 
  f.plotOn(framecosk); 

  // Draw the frame on the canvas
  TCanvas* c = new TCanvas("c", "c", 400, 400); 
  framecosk->Draw(); 
  TString pdffile = outfile + ".pdf"; 
  c->Print(pdffile); 
}


int main(int argc, char** argv) {
  TString func     = argv[1]; 
  if (func == "bmass") {
    TString datatype = argv[2]; 
    TString label    = argv[3]; 
    TString cut      = argv[4]; 
    TString outfile  = argv[5]; 
    
    bmass(datatype, label, cut, outfile); 

  } else if (func == "fl"){

    TString outfile  = argv[2];     
    fl(outfile);
 
  } else { 
    cerr << "No function available for: " << func.Data() << endl; 
  }
  gSystem->Exit(0);

  return 0 ;
}

