// -----------------------------------------------
//       Author: Xin Shi <Xin.Shi@cern.ch> 
//       Created:   [2013-08-15 Thu 14:54] 
// -----------------------------------------------
#include <TSystem.h>
#include <TH1.h>
#include <TFile.h>
#include <TPad.h> 
#include <TCanvas.h> 

#include <RooRealVar.h>
#include <RooGaussian.h>
#include <RooChebychev.h> 
#include <RooAddPdf.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooPlot.h>
// #include <RooArgList.h>

#include "tools.cc" 

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

  // Create two Gaussian PDFs g1(x,mean1,sigma) anf g2(x,mean2,sigma)
  // and their paramaters

  RooRealVar mean("mean","mean of gaussians", 5.27, 5.23, 5.32) ;
  RooRealVar sigma("sigma","width of gaussians", 0.0285, 0, 1) ;
  // RooRealVar sigma2("sigma2","width of gaussians", 0.0663, 0, 1) ;

  RooGaussian sig("sig","Signal component", x, mean, sigma) ;  
  // RooGaussian sig2("sig2","Signal component 2",x,mean,sigma2) ;  

  // // mean.setConstant(false); 
  // // sigma1.setConstant(false); 
  // // sigma2.setConstant(false); 
  
  // Build Chebychev polynomial p.d.f.  
  RooRealVar a0("a0", "constant", 0.5, -1, 1.) ;
  RooRealVar a1("a1", "linear", 0.6, -1, 1) ;
  RooRealVar a2("a2", "quadratic", 0.1, -1, 1) ;
  // // //   RooRealVar a3("a3", "oct", 0.1, -1, 1) ;
  
  RooChebychev bkg("bkg", "Background", x, RooArgSet(a0, a1, a2)) ;
  // // //  RooChebychev bkg("bkg", "Background", x, RooArgSet(a0, a1, a2, a3)) ;



  // // Build exponational background p.d.f.
  // stringstream myString;
  // myString << "exp(-" << x.getPlotLabel() << "/tau1)";

  // cout << ">>> myString: " << myString.str().c_str() << endl;
  
  // RooRealVar tau1("tau1","First background time constant",0.19, 0, 1);
  
  // RooGenericPdf BkgMassExp1("BkgMassExp1",myString.str().c_str(),
  // 			    RooArgSet(x,tau1));
  // RooAddPdf bkg("BkgMassComb","Background mass comb. bkg pdf",
  // 		RooArgList(BkgMassExp1), RooArgList());
  
  // Sum the signal components into a composite signal p.d.f.
  // RooRealVar sig1frac("sig1frac","fraction of component 1 in signal", 0.6, 0, 1) ;

  // RooAddPdf sig("sig", "Signal", RooArgList(sig1,sig2), sig1frac) ;

  // // Sum the composite signal and background 
  // // RooRealVar bkgfrac("bkgfrac","fraction of background",0.5,0.,1.) ;
  // // RooAddPdf  model("model","g1+g2+a",RooArgList(bkg,sig),bkgfrac) ;

  // Construct signal+background PDF
  RooRealVar nsig("nsig", "number of signal events", 4648, 0, 1000000); 
  RooRealVar nbkg("nbkg", "number of background events", 21472, 0, 1000000);
  RooAddPdf  model("model", "g+c", RooArgList(bkg, sig), RooArgList(nbkg, nsig)) ;
  
  // Print structure of composite p.d.f.
  model.Print("t") ;

  // // RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  // RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

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

  // // Overlay the background+sig2 components of model with a dotted line
  // // model.plotOn(xframe,Components("bkg,sig2"),LineStyle(kDotted)) ;

  // Draw the frame on the canvas
  // new TCanvas("BMass", "BMass", 400, 400) ;
  TCanvas* c = new TCanvas("c", "c", 400, 400); 
  set_root_style(); 
  c->UseCurrentStyle() ;

  gPad->SetLeftMargin(0.15) ;
  xframe->GetYaxis()->SetTitleOffset(1.7) ; 
  xframe->Draw();

  TString pdffile = outfile; 
  // gPad->Print(pdffile); 
  c->Print(pdffile); 
  
  // Persist fit result in root file 
  // -------------------------------------------------------------
  
  // Open new ROOT file save result 
  TString resfile = "test.root"; 
  TFile resf(resfile, "RECREATE") ;
  gPad->Write("plot"); 
  if (! test) fitres->Write("fitres") ;
  resf.Close() ;

  // In a clean ROOT session retrieve the persisted fit result as follows:
  // RooFitResult* r = gDirectory->Get("fitres") ;

}

int main(int argc, char** argv) {
  TString func     = argv[1]; 
  TString datatype = argv[2]; 
  TString label    = argv[3]; 
  TString cut      = argv[4]; 
  TString outfile  = argv[5]; 
  
  if (func == "bmass") 
    // bmass(datatype, label, cut, outfile); 
    bmass(datatype, label, cut, outfile); 
  else 
    cerr << "No function available for: " << func.Data() << endl; 

  gSystem->Exit(0);

  return 0 ;
}

