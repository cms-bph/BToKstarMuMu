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

void fl()
{
  // RooRealVar* FlS = new RooRealVar("FlS","F_{L}",0.0,0.0,1.0);
  // RooRealVar* AfbS = new RooRealVar("AfbS","A_{FB}",0.0,-1.0,1.0);
  // FlS->setConstant(false);
  // AfbS->setConstant(false);
  // RooArgSet* Vars; 
  // Vars->add(*FlS);
  // Vars->add(*AfbS);
  // stringstream myString; 
  // myString.str("");
  // S and P-wave decay rate
  // 
  // RooRealVar* x; 
  // myString << "(3/4 * FlS * (1 - " << x->getPlotLabel() << "*"
  // 	   << x->getPlotLabel() << ")";
  
  // from IterativeMassAngleFitq2Bins 
  // RooRealVar* fitVar;

  // Start with example from RooFit
  RooRealVar x("x", "x", -10, 10); 
  RooRealVar y("y", "y", -10, 10); 
  RooRealVar a("a", "a", 5); 
  RooRealVar b("b", "b", 5); 

  RooGenericPdf f("f", "a*x*x+b*y*y-0.3*y*y*y", RooArgSet(x,y,a,b));
  // Generate a 2-dimensional dataset data(x,y) 
  RooDataSet* data = f.generate(RooArgSet(x,y), 10000); 

  // fit the 2-dimensional model f(x,y) to data(x,y)
  f.fitTo(*data); 

  // plot the x distribution of data(x,y) and f(x,y)
  RooPlot* framex = x.frame(); 
  data->plotOn(framex); 
  f.plotOn(framex); 

  // plot the y distribution of data(x,y) and f(x,y)
  RooPlot* framey = y.frame(); 
  data->plotOn(framey); 
  f.plotOn(framey); 

  // Draw the frame on the canvas
  TCanvas* c = new TCanvas("c", "c", 400, 200); 
  
  c->Divide(2);
  c->cd(1); 
  framex->Draw(); 
  c->cd(2);
  framey->Draw();
  
  TString outfile = "test"; 
  
  TString pdffile = outfile + ".pdf"; 

 
  c->Print(pdffile); 
  
}

int main(int argc, char** argv) {
  TString func     = argv[1]; 
  TString datatype = argv[2]; 
  TString label    = argv[3]; 
  TString cut      = argv[4]; 
  TString outfile  = argv[5]; 
  
  if (func == "bmass") 
    bmass(datatype, label, cut, outfile); 

  else if (func == "fl")
    fl(); 

  else 
    cerr << "No function available for: " << func.Data() << endl; 

  gSystem->Exit(0);

  return 0 ;
}

