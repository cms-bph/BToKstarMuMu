//vim: sw=4 ts=4 fdm=marker et:

// -----------------------------------------------
//       Author: Xin Shi <Xin.Shi@cern.ch> 
//       Created:   [2013-08-15 Thu 14:54] 
// -----------------------------------------------
#include <sstream>
#include <math.h>

#include <TSystem.h>
#include <TF1.h>
#include <TF2.h>
#include <TH1.h>
#include <TH2.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMinuit.h>
#include <TFile.h>
#include <TPad.h> 
#include <TCanvas.h> 
#include <TChain.h> 
#include <TPaveText.h>
#include <TLatex.h>
#include <TString.h>
#include <TGraphAsymmErrors.h>
#include <TLorentzVector.h>

#include <RooRealVar.h>
#include <RooAbsPdf.h>
#include <RooAddPdf.h>
#include <RooGaussian.h>
#include <RooChebychev.h> 
#include <RooGenericPdf.h> 
#include <RooExponential.h>
#include <RooPolynomial.h>
#include <RooExtendPdf.h>
#include <RooProdPdf.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooAbsData.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooAddition.h>

#include "tools.h" 

using namespace std; 
using namespace RooFit;

TChain *ch=new TChain("tree");
//Constants, Fit results for efficiency, etc.. //{{{
char genQ2range[10][32] = {"genQ2 < 2.00 && genQ2 > 1.00",
                           "genQ2 < 4.30 && genQ2 > 2.00",
                           "genQ2 < 8.68 && genQ2 > 4.30",
                           "genQ2 <10.09 && genQ2 > 8.68",
                           "genQ2 <12.86 && genQ2 >10.09",
                           "genQ2 <14.18 && genQ2 >12.86",
                           "genQ2 <16.00 && genQ2 >14.18",
                           "genQ2 <19.00 && genQ2 >16.00",
                           "genQ2 <19.00 && genQ2 > 1.00",
                           "genQ2 < 6.00 && genQ2 > 1.00"};
char q2range[10][32] = {"Q2 < 2.00 && Q2 > 1.00",
                        "Q2 < 4.30 && Q2 > 2.00",
                        "Q2 < 8.68 && Q2 > 4.30",
                        "Q2 <10.09 && Q2 > 8.68",
                        "Q2 <12.86 && Q2 >10.09",
                        "Q2 <14.18 && Q2 >12.86",
                        "Q2 <16.00 && Q2 >14.18",
                        "Q2 <19.00 && Q2 >16.00",
                        "Q2 <19.00 && Q2 > 1.00",
                        "Q2 < 6.00 && Q2 > 1.00"};
double q2rangedn[10] = {1.00 , 2.00 , 4.30 , 8.68  , 10.09 , 12.86 , 14.18 , 16.00 ,  1.00 , 1.00};
double q2rangeup[10] = {2.00 , 4.30 , 8.68 , 10.09 , 12.86 , 14.18 , 16.00 , 19.00 , 19.00 , 6.00};
double arrAccPar2011   [8][20] = {{1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                              {1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                              {1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                              {1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                              {1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                              {1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                              {1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                              {1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}};
double arrAccParErr2011[8][20] = {{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                              {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                              {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                              {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                              {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                              {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                              {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                              {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}};
double arrRecPar2011   [8][20] = {
            {0.0017687,-8.73499e-06,-0.000935262,-0.000787674,
            -0.00529644,0.00155626,0.00356542,0.000171319,
            0,0,0,0,
            0.00551706,-0.00339678,-0.00463866,0.00240596,
            -0.00197731,0.00184597,0.00200529,-0.00178549},
            {0.00157318,7.37446e-05,-0.000749493,-0.000781159,
            -0.00288927,0.000261854,0.00187545,0.000883533,
            0,0,0,0,
            0.00132593,-0.000331488,-0.00112174,-0.000108788,
            0,0,0,0},
            {0.00118453,-1.23943e-05,-0.000454793,-0.000541742,
            -2.12267e-05,0.000714674,-9.69142e-05,-0.000768923,
            -7.59472e-05,0.000229203,5.36592e-05,-6.0488e-05,
            -0.00106615,-0.000896725,0.000710503,0.00156238,
            0,0,0,0},
            {0.00106852,-7.35469e-05,-0.000403409,-0.00035673,
            0.000637779,0.000294361,-0.000303145,-0.000885132,
            5.60979e-05,8.91354e-05,-0.000103992,-5.73688e-05,
            -0.00112888,9.3696e-05,0.000728873,0.000650714,
            0,0,0,0},
            {0.00102525,-0.000118438,-0.000346987,-0.000315243,
            0.000741172,0.000462559,-0.000313513,-0.00064603,
            -0.0007614,0.00153595,0.000840568,-0.0014875,
            -1.66019e-05,-0.00190496,-0.000546402,0.00220778,
            0,0,0,0},
            {0.00108772,0.000130305,-0.000457955,-0.000514284,
            0.000225576,-0.00108671,0.000283481,0.000592928,
            4.66766e-05,-0.00017112,4.27933e-05,0.000247925,
            7.59011e-05,0.000852078,-0.000514188,-8.52169e-05,
            0,0,0,0},
            {0.000964546,4.77019e-05,-0.000249802,-0.00035122,
            0.00120212,-0.0013481,-0.00081191,0.00109161,
            -0.000838523,0.00118858,0.000579591,-0.00101538,
            0,0,0,0,
            0,0,0,0},
            {0.00103644,7.27679e-05,-0.000366081,-0.000162294,
            0.000571208,-0.00089218,-1.99667e-05,0.00066918,
            -0.000423471,0.000360459,0.000228834,-0.000207738,
            0,0,0,0,
            0,0,0,0}
                            };
double arrRecParErr2011[8][20] = {
            {9.77352e-06,1.45679e-05,1.44267e-05,1.71613e-05,
            1.7851e-05,2.42707e-05,2.66439e-05,2.89112e-05,
            0,0,0,0,
            2.06968e-05,2.87918e-05,3.2287e-05,3.52259e-05,
            1.71744e-05,2.50134e-05,2.75925e-05,3.13877e-05},
            {1.93048e-05,4.54403e-05,3.19836e-05,6.34448e-05,
            6.97046e-05,0.000175846,0.000127026,0.000252818,
            0,0,0,0,
            5.60232e-05,0.000145372,0.000105324,0.000210834,
            0,0,0,0},
            {1.20764e-05,2.99275e-05,2.12433e-05,4.24699e-05,
            5.88911e-05,0.000155115,0.0001176,0.000229393,
            1.73594e-05,3.28355e-05,3.60904e-05,5.2698e-05,
            5.76856e-05,0.000147459,0.000126591,0.000229798,
            0,0,0,0},
            {1.07973e-05,2.65747e-05,1.89528e-05,3.78338e-05,
            8.26138e-05,0.000208233,0.000148573,0.000297846,
            1.46533e-05,4.07354e-05,3.25714e-05,6.24532e-05,
            0.000102431,0.000261492,0.00018865,0.000376784,
            0,0,0,0},
            {1.43812e-05,3.69822e-05,2.68571e-05,5.34828e-05,
            0.000109918,0.000289434,0.000212286,0.000422439,
            4.35946e-05,9.76997e-05,7.60719e-05,0.000139613,
            0.000142318,0.000374508,0.000278563,0.000548987,
            0,0,0,0},
            {4.57551e-05,0.000116528,7.84802e-05,0.000164459,
            0.000350888,0.000882504,0.000628837,0.00126522,
            6.95431e-05,0.000192566,0.000150432,0.000293712,
            0.000438904,0.00111066,0.000809841,0.00160928,
            0,0,0,0},
            {1.41088e-05,3.73896e-05,2.85699e-05,5.60619e-05,
            5.72502e-05,0.000142174,0.000109448,0.000209836,
            5.66788e-05,0.000139079,0.000107494,0.00020466,
            0,0,0,0,
            0,0,0,0},
            {1.32107e-05,3.6017e-05,2.71503e-05,5.4467e-05,
            3.99561e-05,0.000108688,8.82189e-05,0.000169816,
            3.68927e-05,0.000100146,8.35791e-05,0.000158447,
            0,0,0,0,
            0,0,0,0}
};

double arrAccPar2012   [8][20] = {{1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                              {1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                              {1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                              {1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                              {1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                              {1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                              {1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                              {1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}};
double arrAccParErr2012[8][20] = {{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                              {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                              {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                              {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                              {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                              {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                              {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                              {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}};
double arrRecPar2012   [8][20] = {{1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                              {1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                              {1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                              {1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                              {1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                              {1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                              {1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                              {1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}};
double arrRecParErr2012[8][20] = {{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                              {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                              {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                              {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                              {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                              {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                              {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                              {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}};
//}}}

//_________________________________________________________________________________

void bmass( const char outfile[] = "bmass")
{//{{{
  bool test = false; 

  RooRealVar x("Bmass", "B^{+/-} mass(GeV/c^{2})", 5.27953-0.28, 5.27953+0.28) ;
  RooDataSet data("data", "data", RooArgSet(x), Import(*ch)) ;
  data.Print();

  // Create model and dataset
  // -------------------------------------------------------------------------
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
  // ------------------------------------------------------------------------
  RooFitResult* fitres; 
  if (! test) {
    fitres = model.fitTo(data, Extended(kTRUE), Save(kTRUE)) ;
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
  paveText->SetBorderSize(0);
  paveText->SetFillColor(kWhite);
  paveText->AddText(Form("nsig = %.0f #pm %.0f ", nsig.getVal(), nsig.getError())); 
  paveText->AddText(Form("nbkg = %.0f #pm %.0f ", nbkg.getVal(), nbkg.getError())); 
  paveText->AddText(Form("mean = %.3f #pm %.3f ", mean.getVal(), mean.getError())); 
  paveText->AddText(Form("sigma = %.3f #pm %.3f ", sigma.getVal(), sigma.getError())); 
  paveText->Draw(); 

  c->Print(TString::Format("./plots/%s.pdf",outfile));
  
  // Persist fit result in root file 
  // -------------------------------------------------------------
  TFile resf(TString::Format("./plots/%s.root",outfile), "RECREATE") ;
  gPad->Write("plot"); 
  if (! test) fitres->Write("fitres") ;
  resf.Close() ;

  // In a clean ROOT session retrieve the persisted fit result as follows:
  // RooFitResult* r = gDirectory->Get("fitres") ;
  
  delete paveText; 
  delete c;

}//}}}

//_________________________________________________________________________________

std::vector<double> fl_bin(int iBin, const char outfile[] = "fl")
{//{{{
  // From fomula (2) in LHCb 2012 PRL108, 181806(2012)
  // integrated over theta_l and phi: 
  // 
  // 1/Gamma * d^2 Gamma/d cos(theta_K) dq^2 = 3/2 * F_L cos^2(theta_K)
  // + 3/4(1-F_L)(1-cos^2theta_K)
  // 

  bool test = false;

  RooRealVar genCosThetaK("genCosThetaK", "cos#theta_{K}", -1, 1);
  RooRealVar Q2("Q2","q^{2}",0.5,20.);
  RooRealVar fl("fl", "F_{L}", 0.5, -0.5, 1.5);

  RooGenericPdf f("f", "1.5*fl*genCosThetaK*genCosThetaK+0.75*(1-fl)*(1-genCosThetaK*genCosThetaK)", RooArgSet(genCosThetaK,fl));
  RooDataSet* data;
  
  if (test){
      fl.setVal(0.5);
      data = f.generate(RooArgSet(genCosThetaK,Q2), 10000);
  }else{
      data = new RooDataSet("data","data",ch,RooArgSet(genCosThetaK,Q2),q2range[iBin],0);
  }
  
  //f.fitTo(*data); 
  f.fitTo(*data,Extended(kTRUE)); 

  RooPlot* framecosk = genCosThetaK.frame(); 
  data->plotOn(framecosk); 
  f.plotOn(framecosk); 

  // Draw the frame on the canvas
  TCanvas *c = new TCanvas("c"); 
  framecosk->SetTitle("");
  framecosk->Draw();

  TLatex *t1 = new TLatex();
  t1->SetNDC();
  t1->DrawLatex(.40,.85,TString::Format("%s",q2range[iBin]));
  t1->DrawLatex(.40,.79,TString::Format("F_{L}=%5.3f#pm%5.3f",fl.getVal(),fl.getError()));

  c->Print(TString::Format("./plots/%s_bin%d.pdf",outfile,iBin));

  delete c;
  delete t1;
  delete data;

  std::vector<double> outvect;
  outvect.push_back(fl.getVal());
  outvect.push_back(fl.getError());
  return outvect;
}//}}}

void fl(const char outfile[] = "fl")
{//{{{

    TCanvas *c = new TCanvas();
    TH2F *frame = new TH2F("frame","",18,1,19,10,0,1.5);
    frame->SetStats(kFALSE);
    frame->SetXTitle("q^{2} [(GeV)^{2}]");
    frame->SetYTitle("F_{L}");
    frame->Draw();

    double x[8]={1.5,3.15,6.49,9.385,11.475,13.52,15.09,17.5};
    double xerr[8]={0.5,1.15,2.09,0.705,1.385,0.66,0.91,1.5};
    double yfl[8],yerrfl[8];

    std::vector<double> vbin;
    for(int ibin = 0; ibin < 8; ibin++){
        vbin = fl_bin(ibin);
        yfl[ibin]       =vbin.at(0);
        yerrfl[ibin]    =vbin.at(1);
    }
    
    // Check input data
    for(int ibin = 0; ibin < 8; ibin++){
        printf("yfl [%d]=%6.4f +- %6.4f\n",ibin,yfl[ibin],yerrfl[ibin]);
    }

    TGraphAsymmErrors *g_fl  = new TGraphAsymmErrors(8,x,yfl,xerr,xerr,yerrfl,yerrfl);
    g_fl->Draw("P*");
    c->Print(TString::Format("./plots/%s.pdf",outfile));

    delete g_fl;
    delete frame;
    delete c;
}//}}}

//_________________________________________________________________________________

std::vector<double> angular_gen_bin(int iBin, const char outfile[] = "angular_gen")
{//{{{

    RooRealVar genCosThetaK("genCosThetaK", "cos#theta_{K}", -1., 1.);
    RooRealVar genCosThetaL("genCosThetaL", "cos#theta_{L}", -1., 1.);
    RooRealVar genQ2("genQ2","q^{2}",0.5,20.);
    RooRealVar fl("fl", "F_{L}", 0.8, -0.2, 1.2);
    RooRealVar afb("afb", "A_{FB}", 0., -1., 1.);
    RooRealVar fs("fs","F_{S}",0.);//Derive from B0ToKstarJpsi
    RooRealVar as("as","A_{S}",0.);//Derive from B0ToKstarJpsi
    fs.setConstant(kTRUE);
    as.setConstant(kTRUE);

    RooRealVar nsig("nsig","nsig",1E4,1E2,1E8);
    RooRealVar nbkg("nbkg","nbkg",10,0.1,1E3);
    
    RooGenericPdf f_sig("f_sig", "9/16*((2/3*fs+4/3*as*genCosThetaK)*(1-genCosThetaL*genCosThetaL)+(1-fs)*(2*fl*genCosThetaK*genCosThetaK*(1-genCosThetaL*genCosThetaL)+1/2*(1-fl)*(1-genCosThetaK*genCosThetaK)*(1+genCosThetaL*genCosThetaL)+4/3*afb*(1-genCosThetaK*genCosThetaK)*genCosThetaL))", RooArgSet(genCosThetaK,genCosThetaL,fl,afb,fs,as));
    RooGenericPdf f_bkg("f_bkg", "1",RooArgSet());
    //nbkg.setConstant(kTRUE);
    RooAddPdf f("f","f",RooArgList(f_sig,f_bkg),RooArgList(nsig,nbkg));
    RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(genCosThetaK,genCosThetaL,genQ2),genQ2range[iBin],0);

    RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit2"));

    // Draw the frame on the canvas
    TCanvas* c = new TCanvas("c");
    RooPlot* framecosk = genCosThetaK.frame(); 
    data->plotOn(framecosk,Binning(100)); 
    f.plotOn(framecosk); 
    f.plotOn(framecosk,Components(f_bkg),LineStyle(2),LineColor(8),LineWidth(2));

    framecosk->SetTitle("");
    framecosk->SetMinimum(0);
    framecosk->Draw();

    TLatex *t1 = new TLatex();
    t1->SetNDC();
    double fixNDC = -0.5;
    if (iBin < 5) fixNDC = 0.;
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",genQ2range[iBin]));
    //t1->DrawLatex(.35,.80+fixNDC,TString::Format("F_{L}=%5.3f#pm%5.3f",fl.getVal(),fl.getError()));
    //t1->DrawLatex(.35,.74+fixNDC,TString::Format("A_{FB}=%5.3f#pm%5.3f",afb.getVal(),afb.getError()));
    if ( f_fitresult->status() == 0){
        t1->DrawLatex(.35,.68+fixNDC,TString::Format("Fit status: %s","GOOD"));
        //t1->DrawLatex(.50,.68,TString::Format("Fit status: %d",f_fitresult->covQual()));
    }else{
        t1->DrawLatex(.35,.68+fixNDC,TString::Format("Fit status: %s(%d)","BAD",f_fitresult->status()));
        //t1->DrawLatex(.50,.68,TString::Format("Fit status: %d",f_fitresult->covQual()));
    }
    //t1->DrawLatex(.35,.10,TString::Format("nbkg=%5.3f#pm%5.3f",nbkg.getVal(),nbkg.getError()));
    c->Print(TString::Format("./plots/%s_cosk_gen_bin%d.pdf",outfile,iBin));

    //
    RooPlot* framecosl = genCosThetaL.frame(); 
    data->plotOn(framecosl,Binning(100)); 
    f.plotOn(framecosl); 
    //f.plotOn(framecosk,Components(f_sig),LineColor(4),LineWidth(2));
    f.plotOn(framecosl,Components(f_bkg),LineStyle(2),LineColor(8),LineWidth(2));

    framecosl->SetTitle("");
    framecosl->SetMinimum(0);
    framecosl->Draw();

    if (iBin < 5) fixNDC = 0.;
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",genQ2range[iBin]));
    c->Update();
    c->Print(TString::Format("./plots/%s_cosl_gen_bin%d.pdf",outfile,iBin));

    // Make 2-D plot
    TH1 *h1 = data->createHistogram("genCosThetaL,genCosThetaK", 100, 100);
    h1->Draw("LEGO2");
    c->Update();
    c->Print(TString::Format("./plots/%s_2d_gen_bin%d.pdf",outfile,iBin));

    // clear
    delete t1;
    delete c;
    delete data;

    //write output
    std::vector<double> output;
    output.push_back(fl.getVal());
    output.push_back(fl.getError());
    output.push_back(afb.getVal());
    output.push_back(afb.getError());
    return output;

}//}}}

void angular_gen(const char outfile[] = "angular_gen")
{//{{{

    TCanvas *c = new TCanvas();
    TH2F *frame = new TH2F("frame","",18,1,19,10,-1,1);
    frame->SetStats(kFALSE);
    frame->SetXTitle("q^{2} [(GeV)^{2}]");
    frame->SetYTitle("F_{L}");
    frame->SetAxisRange(0,1,"Y");
    frame->Draw();

    double x[8]={1.5,3.15,6.49,9.385,11.475,13.52,15.09,17.5};
    double xerr[8]={0.5,1.15,2.09,0.705,1.385,0.66,0.91,1.5};
    double yafb[8],yerrafb[8],yfl[8],yerrfl[8];

    std::vector<double> vbin;
    for(int ibin = 0; ibin < 8; ibin++){
        vbin = angular_gen_bin(ibin);
        yfl[ibin]       =vbin.at(0);
        yerrfl[ibin]    =vbin.at(1);
        yafb[ibin]      =vbin.at(2);
        yerrafb[ibin]   =vbin.at(3);
    }
    
    // Check input data
    for(int ibin = 0; ibin < 8; ibin++){
        printf("yafb[%d]=%6.4f +- %6.4f\n",ibin,yafb[ibin],yerrafb[ibin]);
        printf("yfl [%d]=%6.4f +- %6.4f\n",ibin,yfl[ibin],yerrfl[ibin]);
    }

    TGraphAsymmErrors *g_fl  = new TGraphAsymmErrors(8,x,yfl,xerr,xerr,yerrfl,yerrfl);
    g_fl->SetFillColor(2);
    g_fl->SetFillStyle(3001);
    g_fl->Draw("P2");
    c->Print(TString::Format("./plots/%s_fl.pdf",outfile));
    c->Clear();

    frame->SetTitle("");
    frame->SetYTitle("A_{FB}");
    frame->SetXTitle("q^{2} [(GeV)^{2}]");
    frame->SetAxisRange(-1,1,"Y");
    frame->Draw();
    TGraphAsymmErrors *g_afb = new TGraphAsymmErrors(8,x,yafb,xerr,xerr,yerrafb,yerrafb);
    g_afb->SetFillColor(2);
    g_afb->SetFillStyle(3001);
    g_afb->Draw("P2");
    c->Print(TString::Format("./plots/%s_afb.pdf",outfile));
}//}}}

//_________________________________________________________________________________
//Fit parameters of acceptance and efficiency using TMinuit instead of RooFit.

TF2  *f2_fcn = NULL;
double model_2D(double *x, double *par)
{//{{{
    double xx = x[0];
    double yy = x[1];
    for (int i = 0; i < f2_fcn->GetNpar(); i++) f2_fcn->SetParameter(i,par[i]);
    return f2_fcn->Eval(xx,yy);
}//}}}

double model_2D_squared(double *x, double *par)
{//{{{
    double xx = x[0];
    double yy = x[1];
    for (int i = 0; i < f2_fcn->GetNpar(); i++) f2_fcn->SetParameter(i,par[i]);
    return pow(f2_fcn->Eval(xx,yy),2);
}//}}}

TH2F *h2_fcn = NULL;
// nParameters, ???, return fcn value, parameter array, strategy
void fcn_binnedChi2_2D(int &npar, double *gin, double &f, double *par, int iflag)
{//{{{
    f=0;
    for (int i = 1; i <= h2_fcn->GetNbinsX(); i++) {
        for (int j = 1; j <= h2_fcn->GetNbinsY(); j++) {
            int gBin = h2_fcn->GetBin(i,j);
            double x[2] = {h2_fcn->GetXaxis()->GetBinCenter(i),h2_fcn->GetYaxis()->GetBinCenter(j)};
            double measure  = h2_fcn->GetBinContent(gBin);
            double error    = h2_fcn->GetBinError(gBin);
            
            //// Naively test using center value
            //double func     = model_2D(x, par);//Take center value
            //double delta    = (measure-func)/error;
            //f+=delta*delta;
            
            //// Real run using integral
            TF2 *f2_square = new TF2("f2_square",model_2D_squared,h2_fcn->GetXaxis()->GetBinLowEdge(i),h2_fcn->GetXaxis()->GetBinUpEdge(i),h2_fcn->GetYaxis()->GetBinLowEdge(j),h2_fcn->GetYaxis()->GetBinUpEdge(j),f2_fcn->GetNpar());
            for (int i = 0; i < f2_fcn->GetNpar(); i++){//nPar MUST be the same value as f2_fcn
                f2_fcn->SetParameter(i,par[i]);
                f2_square->SetParameter(i,par[i]);
            }
            double xi = h2_fcn->GetXaxis()->GetBinLowEdge(i);
            double xf = h2_fcn->GetXaxis()->GetBinUpEdge(i);
            double yi = h2_fcn->GetYaxis()->GetBinLowEdge(j);
            double yf = h2_fcn->GetYaxis()->GetBinUpEdge(j);
            f += (f2_square->Integral(xi,xf,yi,yf)-2*f2_fcn->Integral(xi,xf,yi,yf)*measure+measure*measure*(xf-xi)*(yf-yi))/(error*error);

            f2_square = 0;
            delete f2_square;
        }
    }
    //printf("FCN in calls = %f\n",f);
    //printf("npar=%d ",npar);
}//}}}

std::vector<double> acceptance(int iBin) // acceptance
{//{{{

    double gQ2 = 0;
    double gCosThetaK = 0;
    double gCosThetaL = 0;
    double gmuppt = 0;
    double gmupeta= 0;
    double gmumpt = 0;
    double gmumeta= 0;

    ch->SetBranchStatus("*",0);
    ch->SetBranchStatus("genQ2"         , 1);
    ch->SetBranchStatus("genCosTheta*"  , 1);
    ch->SetBranchStatus("genMu*"        , 1);
    ch->SetBranchAddress("genQ2"        , &gQ2);
    ch->SetBranchAddress("genCosThetaK" , &gCosThetaK);
    ch->SetBranchAddress("genCosThetaL" , &gCosThetaL);
    ch->SetBranchAddress("genMupPt"     , &gmuppt);
    ch->SetBranchAddress("genMupEta"    , &gmupeta);
    ch->SetBranchAddress("genMumPt"     , &gmumpt);
    ch->SetBranchAddress("genMumEta"    , &gmumeta);

    // Fill histograms
    float thetaKBins[6]={-1,-0.7,0.,0.4,0.8,1};
    float thetaLBins[7]={-1,-0.7,-0.3,0.,0.3,0.7,1};
    TH2F h2_ngen("h2_ngen","h2_ngen",6,thetaLBins,5,thetaKBins);
    TH2F h2_nacc("h2_nacc" ,"h2_nacc" ,6,thetaLBins,5,thetaKBins); 
    for (int entry = 0; entry < ch->GetEntries(); entry++) {
        ch->GetEntry(entry);
        if (gQ2 > q2rangeup[iBin] || gQ2 < q2rangedn[iBin]) continue;
        h2_ngen.Fill(gCosThetaL,gCosThetaK);
        if ( fabs(gmumeta) < 2.3 && fabs(gmupeta) < 2.3 && gmumpt > 2.5 && gmuppt > 2.5 ) h2_nacc.Fill(gCosThetaL,gCosThetaK);
    }
    
    // Calculate acceptance
    TH2F h2_acc("h2_acc","h2_acc",6,thetaLBins,5,thetaKBins);
    h2_acc.SetAxisRange(0.,1.,"Z");
    for (int i = 1; i <= 6; i++) {
        for (int j = 1; j <= 5; j++) {
            if (h2_ngen.GetBinContent(i,j) == 0 || h2_nacc.GetBinContent(i,j) == 0) {
                printf("WARNING: Angular bin(%d,%d)=0, set error to be 1.\n",i,j);
                h2_acc.SetBinContent(i,j,0);
                h2_acc.SetBinError(i,j,1);
            }else{
                h2_acc.SetBinContent(i,j,h2_nacc.GetBinContent(i,j)/h2_ngen.GetBinContent(i,j));
                h2_acc.SetBinError(i,j,sqrt(h2_acc.GetBinContent(i,j)*(1-h2_acc.GetBinContent(i,j))));
            }
        }
    }
    
    // Using pure TMinuit
    int nPar = 20;
    TMinuit *gMinuit = new TMinuit(nPar);
    h2_fcn = &h2_acc;
    gMinuit->SetFCN(fcn_binnedChi2_2D);
    
    TF2 f2_model("f2_model","([0]+[1]*y+[2]*(3*y**2-1)/2+[3]*(5*y**3-3*y)/2)+([4]+[5]*y+[6]*(3*y**2-1)/2+[7]*(5*y**3-3*y)/2)*x**2+([8]+[9]*y+[10]*(3*y**2-1)/2+[11]*(5*y**3-3*y)/2)*x**3+([12]+[13]*y+[14]*(3*y**2-1)/2+[15]*(5*y**3-3*y)/2)*x**4+([16]+[17]*y+[18]*(3*y**2-1)/2+[19]*(5*y**3-3*y)/2)*x**6",-1.,1.,-1.,1.);
    f2_fcn = &f2_model;
    gMinuit->DefineParameter( 0, "k0l0",  .01,  1E-3,    -1E-2, 1E-1);
    gMinuit->DefineParameter( 1, "k1l0",   0.,  1E-3,    -1E-1, 1E-1);
    gMinuit->DefineParameter( 2, "k2l0",   0.,  1E-3,    -1E-1, 1E-1);
    gMinuit->DefineParameter( 3, "k3l0",   0.,  1E-3,    -1E-1, 1E-1);
    gMinuit->DefineParameter( 4, "k0l2", 1E-2,  1E-3,    -1E-1, 1E-1);
    gMinuit->DefineParameter( 5, "k1l2",   0.,  1E-3,    -1E-1, 1E-1);
    gMinuit->DefineParameter( 6, "k2l2",   0.,  1E-3,    -1E-1, 1E-1);
    gMinuit->DefineParameter( 7, "k3l2",   0.,  1E-3,    -1E-1, 1E-1);
    gMinuit->DefineParameter( 8, "k0l3",   0.,  1E-3,    -1E-1, 1E-1);
    gMinuit->DefineParameter( 9, "k1l3",   0.,  1E-3,    -1E-1, 1E-1);
    gMinuit->DefineParameter(10, "k2l3",   0.,  1E-3,    -1E-1, 1E-1);
    gMinuit->DefineParameter(11, "k3l3",   0.,  1E-3,    -1E-1, 1E-1);
    gMinuit->DefineParameter(12, "k0l4",-1E-2,  1E-3,    -1E-1, 1E-1);
    gMinuit->DefineParameter(13, "k1l4",   0.,  1E-3,    -1E-1, 1E-1);
    gMinuit->DefineParameter(14, "k2l4",   0.,  1E-3,    -1E-1, 1E-1);
    gMinuit->DefineParameter(15, "k3l4",   0.,  1E-3,    -1E-1, 1E-1);
    gMinuit->DefineParameter(16, "k0l6", 1E-2,  1E-3,    -1E+1, 1E+1);
    gMinuit->DefineParameter(17, "k1l6",   0.,  1E-3,    -1E+1, 1E+1);
    gMinuit->DefineParameter(18, "k2l6",   0.,  1E-3,    -1E+1, 1E+1);
    gMinuit->DefineParameter(19, "k3l6",   0.,  1E-3,    -1E+1, 1E+1);
    if (iBin == 0) {
        gMinuit->Command("SET PARM 9 0");
        gMinuit->Command("SET PARM 10 0");
        gMinuit->Command("SET PARM 11 0");
        gMinuit->Command("SET PARM 12 0");
        gMinuit->Command("FIX 9");
        gMinuit->Command("FIX 10");
        gMinuit->Command("FIX 11");
        gMinuit->Command("FIX 12");
    }else if (iBin == 1) {
        gMinuit->Command("SET PARM 9 0");
        gMinuit->Command("SET PARM 10 0");
        gMinuit->Command("SET PARM 11 0");
        gMinuit->Command("SET PARM 12 0");
        gMinuit->Command("SET PARM 17 0");
        gMinuit->Command("SET PARM 18 0");
        gMinuit->Command("SET PARM 19 0");
        gMinuit->Command("SET PARM 20 0");
        gMinuit->Command("FIX 9");
        gMinuit->Command("FIX 10");
        gMinuit->Command("FIX 11");
        gMinuit->Command("FIX 12");
        gMinuit->Command("FIX 17");
        gMinuit->Command("FIX 18");
        gMinuit->Command("FIX 19");
        gMinuit->Command("FIX 20");
    }else if (iBin > 1 && iBin < 6 ) {
        gMinuit->Command("SET PARM 17 0");
        gMinuit->Command("SET PARM 18 0");
        gMinuit->Command("SET PARM 19 0");
        gMinuit->Command("SET PARM 20 0");
        gMinuit->Command("FIX 17");
        gMinuit->Command("FIX 18");
        gMinuit->Command("FIX 19");
        gMinuit->Command("FIX 20");
    }else{
        gMinuit->Command("SET PARM 13 0");
        gMinuit->Command("SET PARM 14 0");
        gMinuit->Command("SET PARM 15 0");
        gMinuit->Command("SET PARM 16 0");
        gMinuit->Command("SET PARM 17 0");
        gMinuit->Command("SET PARM 18 0");
        gMinuit->Command("SET PARM 19 0");
        gMinuit->Command("SET PARM 20 0");
        gMinuit->Command("FIX 13");
        gMinuit->Command("FIX 14");
        gMinuit->Command("FIX 15");
        gMinuit->Command("FIX 16");
        gMinuit->Command("FIX 17");
        gMinuit->Command("FIX 18");
        gMinuit->Command("FIX 19");
        gMinuit->Command("FIX 20");
    }
    
    gMinuit->Command("MINI");
    gMinuit->Command("MINI");
    gMinuit->Command("IMPROVE");
    gMinuit->Command("MINOS");

    double arrPar[nPar];
    double arrParErr[nPar];
    for (int iPar = 0; iPar < nPar; iPar++) gMinuit->GetParameter(iPar,arrPar[iPar],arrParErr[iPar]);
    for (int iPar = 0; iPar < nPar; iPar++) f2_model.SetParameter(iPar,arrPar[iPar]);

    // Prepare draw
    TCanvas canvas("canvas");
    TLatex *latex = new TLatex();
    
    // Draw efficiency
    h2_acc.SetStats(0);
    h2_acc.SetMaximum(.02);
    h2_acc.Draw("LEGO2");
    latex->DrawLatexNDC(0.35,0.95,TString::Format("Acceptance_{RECO} in Bin%d",iBin));
    
    // Draw FitResult
    f2_model.SetTitle("");
    f2_model.SetMaximum(.02);
    f2_model.SetLineWidth(1);
    f2_model.Draw("SURF SAME ");
    canvas.Print(TString::Format("./plots/acceptance_2D_bin%d.pdf",iBin));
    

    //// Draw compare
    double chi2Val=0;
    fcn_binnedChi2_2D(nPar, 0, chi2Val, arrPar, 0);
    printf("Chi2=%f \n",chi2Val);
    TH2F h2_compFit("h2_compFit","",6,thetaLBins,5,thetaKBins);
    for (int i = 1; i <= 6; i++) {//thetaL
        for (int j = 1; j <= 5; j++) {//thetaK
            if (h2_acc.GetBinContent(i,j) != 0){
                h2_compFit.SetBinContent(i,j,f2_model.Eval(h2_acc.GetXaxis()->GetBinCenter(i),h2_acc.GetYaxis()->GetBinCenter(j))/h2_acc.GetBinContent(i,j));
            }else{
                h2_compFit.SetBinContent(i,j,0.);
            }
        }
    }
    h2_compFit.SetMinimum(0.);
    h2_compFit.SetStats(0);
    h2_compFit.Draw("LEGO2");
    latex->DrawLatexNDC(0.01,0.95,TString::Format("#chi^{2} = %f",chi2Val));
    latex->DrawLatexNDC(0.3,0.95,TString::Format("#varepsilon_{RECO,measured} / #varepsilon_{RECO,fit} in Bin%d",iBin));
    canvas.Update();
    canvas.Print(TString::Format("./plots/acceptance_compFit_2D_bin%d.pdf",iBin));

    // Draw projection to cosThetaK
    TH1D *h_cosk = new TH1D("h_cosk","",5,thetaKBins);
    for (int kBin = 1; kBin <= 5; kBin++) {
        float sumGen = 0.;
        float sumAcc = 0.;
        for ( int lBin = 1; lBin <= 6; lBin++) {
            sumGen+=h2_ngen.GetBinContent(lBin,kBin);
            sumAcc+=h2_nacc.GetBinContent(lBin,kBin);
        }
        h_cosk->SetBinContent(kBin,sumAcc/sumGen);
        h_cosk->SetBinError(kBin,sqrt(sumAcc*(sumGen-sumAcc)/pow(sumGen,3)));
    }
    h_cosk->SetStats(0);
    h_cosk->SetMinimum(0.);
    h_cosk->SetMaximum(0.02);
    h_cosk->Draw();
    latex->DrawLatexNDC(0.32,0.95,TString::Format("Projection to cos#theta_{k} in Bin%d",iBin));

    TF1  *f_cosk = new TF1("f_cosk"
                          ,"([0]+[1]*x+[2]*(3*x**2-1)/2+[3]*(5*x**3-3*x)/2)+([4]+[5]*x+[6]*(3*x**2-1)/2+[7]*(5*x**3-3*x)/2)/3+([8]+[9]*x+[10]*(3*x**2-1)/2+[11]*(5*x**3-3*x)/2)*0+([12]+[13]*x+[14]*(3*x**2-1)/2+[15]*(5*x**3-3*x)/2)/5+([16]+[17]*x+[18]*(3*x**2-1)/2+[19]*(5*x**3-3*x)/2)/7"
                          ,-1.,1.);
    for (int iPar = 0; iPar < nPar; iPar++) f_cosk->SetParameter(iPar,arrPar[iPar]);
    f_cosk->Draw("SAME");
    canvas.Update();
    canvas.Print(TString::Format("./plots/accptance_cosK_bin%d.pdf",iBin));
    
    // Draw projection to cosThetaL
    TH1D *h_cosl = new TH1D("h_cosl","",6,thetaLBins);
    for ( int lBin = 1; lBin <= 6; lBin++) {
        float sumGen = 0.;
        float sumAcc = 0.;
        for (int kBin = 1; kBin <= 5; kBin++) {
            sumAcc+=h2_nacc.GetBinContent(lBin,kBin);
            sumGen+=h2_ngen.GetBinContent(lBin,kBin);
        }
        h_cosl->SetBinContent(lBin,sumAcc/sumGen);
        h_cosl->SetBinError(lBin,sqrt(sumAcc*(sumGen-sumAcc)/pow(sumGen,3)));
    }
    h_cosl->SetStats(0);
    h_cosl->SetMinimum(0.);
    h_cosl->SetMaximum(0.02);
    h_cosl->Draw();
    latex->DrawLatexNDC(0.32,0.95,TString::Format("Projection to cos#theta_{l} in Bin%d",iBin));

    TF1  *f_cosl = new TF1("f_cosl"
                          ,"([0]+[1]*x*0+[2]*0+[3]*0)+([4]+[5]*0+[6]*0+[7]*0)*x**2+([8]+[9]*0+[10]*0+[11]*0)*x**3+([12]+[13]*0+[14]*0+[15]*0)*x**4+([16]+[17]*0+[18]*0+[19]*0)*x**6"
                          ,-1.,1.);
    for (int iPar = 0; iPar < nPar; iPar++) f_cosl->SetParameter(iPar,arrPar[iPar]);
    f_cosl->Draw("SAME");
    canvas.Update();
    canvas.Print(TString::Format("./plots/acceptance_cosL_bin%d.pdf",iBin));

    // Clear
    delete f_cosl;
    delete h_cosl;
    delete f_cosk;
    delete h_cosk;
    delete latex;
    delete gMinuit;

    //prepare output
    std::vector<double> output;
    for (int iPar = 0; iPar < nPar; iPar++){
        output.push_back(arrPar[iPar]);
        output.push_back(arrParErr[iPar]);
    }
    for (int i = 0; i < output.size(); i=i+2) {
        printf("%f +- %f\n",output[i],output[i+1]);
    }
    return output;
}//}}}

std::vector<double> efficiency(int iBin) // reconstruction efficiency
{//{{{
    printf("Evaluate reconstruction efficiency for bin#%d\n",iBin);
    double BMass = 0;
    double gQ2 = 0;
    double gCosThetaK = 0;
    double gCosThetaL = 0;
    double gmuppt = 0;
    double gmupeta= 0;
    double gmumpt = 0;
    double gmumeta= 0;

    ch->SetBranchStatus("*",0);
    ch->SetBranchStatus("Bmass"         , 1);
    ch->SetBranchStatus("genQ2"         , 1);
    ch->SetBranchStatus("genCosTheta*"  , 1);
    ch->SetBranchStatus("genMu*"        , 1);
    ch->SetBranchAddress("Bmass"        , &BMass);
    ch->SetBranchAddress("genQ2"        , &gQ2);
    ch->SetBranchAddress("genCosThetaK" , &gCosThetaK);
    ch->SetBranchAddress("genCosThetaL" , &gCosThetaL);
    ch->SetBranchAddress("genMupPt"     , &gmuppt);
    ch->SetBranchAddress("genMupEta"    , &gmupeta);
    ch->SetBranchAddress("genMumPt"     , &gmumpt);
    ch->SetBranchAddress("genMumEta"    , &gmumeta);
    RooRealVar genCosThetaL("genCosThetaL","genCosThetaL",-1,1);
    RooRealVar genCosThetaK("genCosThetaK","genCosThetaK",-1,1);

    // Fill histograms
    float thetaKBins[6]={-1,-0.7,0.,0.4,0.8,1};
    float thetaLBins[7]={-1,-0.7,-0.3,0.,0.3,0.7,1};
    TH2F h2_nacc("h2_nacc" ,"h2_nacc" ,6,thetaLBins,5,thetaKBins); 
    TH2F h2_nreco("h2_nreco","h2_nreco",6,thetaLBins,5,thetaKBins);
    TH2F h2_nrej("h2_nrej","h2_nrej",6,thetaLBins,5,thetaKBins);
    for (int entry = 0; entry < ch->GetEntries(); entry++) {
        ch->GetEntry(entry);
        if (gQ2 > q2rangeup[iBin] || gQ2 < q2rangedn[iBin]) continue;
        if ( fabs(gmumeta) < 2.3 && fabs(gmupeta) < 2.3 && gmumpt > 2.5 && gmuppt > 2.5 ) h2_nacc.Fill(gCosThetaL,gCosThetaK);
        if (BMass != 0){
            h2_nreco.Fill(gCosThetaL,gCosThetaK);
        }else{
            h2_nrej.Fill(gCosThetaL,gCosThetaK);
        }
    }
    
    // Calculate efficiency
    TH2F h2_rec("h2_rec","",6,thetaLBins,5,thetaKBins);
    h2_rec.SetMinimum(0.);
    h2_rec.SetXTitle("genCosThetaL");
    h2_rec.SetYTitle("genCosThetaK");
    for (int i = 1; i <= 6; i++) {
        for (int j = 1; j <= 5; j++) {
            // Build from MC samples
            if (h2_nacc.GetBinContent(i,j) == 0 || h2_nreco.GetBinContent(i,j) == 0) {
                printf("WARNING: Angular bin(%d,%d)=0, set error to be 1.\n",i,j);
                h2_rec.SetBinContent(i,j,0.);
                h2_rec.SetBinError(i,j,1.);
            }else{
                h2_rec.SetBinContent(i,j,h2_nreco.GetBinContent(i,j)/h2_nacc.GetBinContent(i,j));
                h2_rec.SetBinError(i,j,sqrt(h2_rec.GetBinContent(i,j)*(1-h2_rec.GetBinContent(i,j))/h2_nacc.GetBinContent(i,j)));
            }

            // Build a toy sample
            //h2_rec.SetBinContent(i,j,fabs(sin((thetaLBins[i]+thetaLBins[i-1])*3.1415926/4)/(20.-2*j))+0.1);
            //h2_rec.SetBinError(i,j,sqrt(h2_rec.GetBinContent(i,j)*(1-h2_rec.GetBinContent(i,j))/20));
        }
    }

    // Using pure TMinuit
    int nPar = 20;
    TMinuit *gMinuit = new TMinuit(nPar);
    h2_fcn = &h2_rec;
    gMinuit->SetFCN(fcn_binnedChi2_2D);
    
    // Use Legendre polynomial for better convergance
    // 1,x,(3x^2-1)/2,(5x^3-3x)/2
    //TF2 f2_model("f2_model","([0]+[1]*y+[2]*y**2+[3]*y**3)+([4]+[5]*y+[6]*y**2+[7]*y**3)*x**2+([8]+[9]*y+[10]*y**2+[11]*y**3)*x**3+([12]+[13]*y+[14]*y**2+[15]*y**3)*x**4+([16]+[17]*y+[18]*y**2+[19]*y**3)*x**6",-1.,1.,-1.,1.);
    TF2 f2_model("f2_model","([0]+[1]*y+[2]*(3*y**2-1)/2+[3]*(5*y**3-3*y)/2)+([4]+[5]*y+[6]*(3*y**2-1)/2+[7]*(5*y**3-3*y)/2)*x**2+([8]+[9]*y+[10]*(3*y**2-1)/2+[11]*(5*y**3-3*y)/2)*x**3+([12]+[13]*y+[14]*(3*y**2-1)/2+[15]*(5*y**3-3*y)/2)*x**4+([16]+[17]*y+[18]*(3*y**2-1)/2+[19]*(5*y**3-3*y)/2)*x**6",-1.,1.,-1.,1.);
    f2_fcn = &f2_model;
    gMinuit->DefineParameter( 0, "k0l0",  .01,  1E-3,    -1E-2, 1E-1);
    gMinuit->DefineParameter( 1, "k1l0",   0.,  1E-3,    -1E-1, 1E-1);
    gMinuit->DefineParameter( 2, "k2l0",   0.,  1E-3,    -1E-1, 1E-1);
    gMinuit->DefineParameter( 3, "k3l0",   0.,  1E-3,    -1E-1, 1E-1);
    gMinuit->DefineParameter( 4, "k0l2", 1E-2,  1E-3,    -1E-1, 1E-1);
    gMinuit->DefineParameter( 5, "k1l2",   0.,  1E-3,    -1E-1, 1E-1);
    gMinuit->DefineParameter( 6, "k2l2",   0.,  1E-3,    -1E-1, 1E-1);
    gMinuit->DefineParameter( 7, "k3l2",   0.,  1E-3,    -1E-1, 1E-1);
    gMinuit->DefineParameter( 8, "k0l3",   0.,  1E-3,    -1E-1, 1E-1);
    gMinuit->DefineParameter( 9, "k1l3",   0.,  1E-3,    -1E-1, 1E-1);
    gMinuit->DefineParameter(10, "k2l3",   0.,  1E-3,    -1E-1, 1E-1);
    gMinuit->DefineParameter(11, "k3l3",   0.,  1E-3,    -1E-1, 1E-1);
    gMinuit->DefineParameter(12, "k0l4",-1E-2,  1E-3,    -1E-1, 1E-1);
    gMinuit->DefineParameter(13, "k1l4",   0.,  1E-3,    -1E-1, 1E-1);
    gMinuit->DefineParameter(14, "k2l4",   0.,  1E-3,    -1E-1, 1E-1);
    gMinuit->DefineParameter(15, "k3l4",   0.,  1E-3,    -1E-1, 1E-1);
    gMinuit->DefineParameter(16, "k0l6", 1E-2,  1E-3,    -1E+1, 1E+1);
    gMinuit->DefineParameter(17, "k1l6",   0.,  1E-3,    -1E+1, 1E+1);
    gMinuit->DefineParameter(18, "k2l6",   0.,  1E-3,    -1E+1, 1E+1);
    gMinuit->DefineParameter(19, "k3l6",   0.,  1E-3,    -1E+1, 1E+1);
    if (iBin == 0) {
        gMinuit->Command("SET PARM 9 0");
        gMinuit->Command("SET PARM 10 0");
        gMinuit->Command("SET PARM 11 0");
        gMinuit->Command("SET PARM 12 0");
        gMinuit->Command("FIX 9");
        gMinuit->Command("FIX 10");
        gMinuit->Command("FIX 11");
        gMinuit->Command("FIX 12");
    }else if (iBin == 1) {
        gMinuit->Command("SET PARM 9 0");
        gMinuit->Command("SET PARM 10 0");
        gMinuit->Command("SET PARM 11 0");
        gMinuit->Command("SET PARM 12 0");
        gMinuit->Command("SET PARM 17 0");
        gMinuit->Command("SET PARM 18 0");
        gMinuit->Command("SET PARM 19 0");
        gMinuit->Command("SET PARM 20 0");
        gMinuit->Command("FIX 9");
        gMinuit->Command("FIX 10");
        gMinuit->Command("FIX 11");
        gMinuit->Command("FIX 12");
        gMinuit->Command("FIX 17");
        gMinuit->Command("FIX 18");
        gMinuit->Command("FIX 19");
        gMinuit->Command("FIX 20");
    }else if (iBin > 1 && iBin < 6 ) {
        gMinuit->Command("SET PARM 17 0");
        gMinuit->Command("SET PARM 18 0");
        gMinuit->Command("SET PARM 19 0");
        gMinuit->Command("SET PARM 20 0");
        gMinuit->Command("FIX 17");
        gMinuit->Command("FIX 18");
        gMinuit->Command("FIX 19");
        gMinuit->Command("FIX 20");
    }else{
        gMinuit->Command("SET PARM 13 0");
        gMinuit->Command("SET PARM 14 0");
        gMinuit->Command("SET PARM 15 0");
        gMinuit->Command("SET PARM 16 0");
        gMinuit->Command("SET PARM 17 0");
        gMinuit->Command("SET PARM 18 0");
        gMinuit->Command("SET PARM 19 0");
        gMinuit->Command("SET PARM 20 0");
        gMinuit->Command("FIX 13");
        gMinuit->Command("FIX 14");
        gMinuit->Command("FIX 15");
        gMinuit->Command("FIX 16");
        gMinuit->Command("FIX 17");
        gMinuit->Command("FIX 18");
        gMinuit->Command("FIX 19");
        gMinuit->Command("FIX 20");
    }
    
    gMinuit->Command("MINI");
    gMinuit->Command("MINI");
    gMinuit->Command("IMPROVE");
    gMinuit->Command("MINOS");

    double arrPar[nPar];
    double arrParErr[nPar];
    for (int iPar = 0; iPar < nPar; iPar++) gMinuit->GetParameter(iPar,arrPar[iPar],arrParErr[iPar]);
    for (int iPar = 0; iPar < nPar; iPar++) f2_model.SetParameter(iPar,arrPar[iPar]);
    
    // Prepare draw
    TCanvas canvas("canvas");
    TLatex *latex = new TLatex();
    
    // Draw efficiency
    h2_rec.SetStats(0);
    h2_rec.SetMaximum(.02);
    h2_rec.Draw("LEGO2");
    latex->DrawLatexNDC(0.35,0.95,TString::Format("#varepsilon_{RECO} in Bin%d",iBin));
    
    // Draw FitResult
    f2_model.SetTitle("");
    f2_model.SetMaximum(.02);
    f2_model.SetLineWidth(1);
    f2_model.Draw("SURF SAME ");
    canvas.Print(TString::Format("./plots/recoEff_2D_bin%d.pdf",iBin));

    //// Draw compare
    double chi2Val=0;
    fcn_binnedChi2_2D(nPar, 0, chi2Val, arrPar, 0);
    printf("Chi2(Bin center)=%f \n",chi2Val);
    
    TH2F h2_compFit("h2_compFit","",6,thetaLBins,5,thetaKBins);
    for (int i = 1; i <= 6; i++) {//thetaL
        for (int j = 1; j <= 5; j++) {//thetaK
            if (h2_rec.GetBinContent(i,j) != 0){
                h2_compFit.SetBinContent(i,j,f2_model.Eval(h2_rec.GetXaxis()->GetBinCenter(i),h2_rec.GetYaxis()->GetBinCenter(j))/h2_rec.GetBinContent(i,j));
            }else{
                h2_compFit.SetBinContent(i,j,0.);
            }
        }
    }
    h2_compFit.SetMinimum(0.);
    h2_compFit.SetStats(0);
    h2_compFit.Draw("LEGO2");
    latex->DrawLatexNDC(0.01,0.95,TString::Format("#chi^{2} = %f",chi2Val));
    latex->DrawLatexNDC(0.3,0.95,TString::Format("#varepsilon_{RECO,fit} / #varepsilon_{RECO,measured} in Bin%d",iBin));
    canvas.Update();
    canvas.Print(TString::Format("./plots/recoEff_compFit_2D_bin%d.pdf",iBin));

    // Draw projection to cosThetaK
    TH1D *h_cosk = new TH1D("h_cosk","",5,thetaKBins);
    for (int kBin = 1; kBin <= 5; kBin++) {
        float sumAcc = 0.;
        float sumRec = 0.;
        for ( int lBin = 1; lBin <= 6; lBin++) {
            sumAcc+=h2_nacc.GetBinContent(lBin,kBin);
            sumRec+=h2_nreco.GetBinContent(lBin,kBin);
        }
        h_cosk->SetBinContent(kBin,sumRec/sumAcc);
        h_cosk->SetBinError(kBin,sqrt(sumRec*(sumAcc-sumRec)/pow(sumAcc,3)));
    }
    h_cosk->SetStats(0);
    h_cosk->SetMinimum(0.);
    h_cosk->SetMaximum(0.02);
    h_cosk->Draw();
    latex->DrawLatexNDC(0.32,0.95,TString::Format("Projection to cos#theta_{k} in Bin%d",iBin));

    TF1  *f_cosk = new TF1("f_cosk"
                          ,"([0]+[1]*x+[2]*(3*x**2-1)/2+[3]*(5*x**3-3*x)/2)+([4]+[5]*x+[6]*(3*x**2-1)/2+[7]*(5*x**3-3*x)/2)/3+([8]+[9]*x+[10]*(3*x**2-1)/2+[11]*(5*x**3-3*x)/2)*0+([12]+[13]*x+[14]*(3*x**2-1)/2+[15]*(5*x**3-3*x)/2)/5+([16]+[17]*x+[18]*(3*x**2-1)/2+[19]*(5*x**3-3*x)/2)/7"
                          ,-1.,1.);
    for (int iPar = 0; iPar < nPar; iPar++) f_cosk->SetParameter(iPar,arrPar[iPar]);
    f_cosk->Draw("SAME");
    canvas.Update();
    canvas.Print(TString::Format("./plots/recoEff_cosK_bin%d.pdf",iBin));
    
    // Draw projection to cosThetaL
    TH1D *h_cosl = new TH1D("h_cosl","",6,thetaLBins);
    for ( int lBin = 1; lBin <= 6; lBin++) {
        float sumAcc = 0.;
        float sumRec = 0.;
        for (int kBin = 1; kBin <= 5; kBin++) {
            sumAcc+=h2_nacc.GetBinContent(lBin,kBin);
            sumRec+=h2_nreco.GetBinContent(lBin,kBin);
        }
        h_cosl->SetBinContent(lBin,sumRec/sumAcc);
        h_cosl->SetBinError(lBin,sqrt(sumRec*(sumAcc-sumRec)/pow(sumAcc,3)));
    }
    h_cosl->SetStats(0);
    h_cosl->SetMinimum(0.);
    h_cosl->SetMaximum(0.02);
    h_cosl->Draw();
    latex->DrawLatexNDC(0.32,0.95,TString::Format("Projection to cos#theta_{l} in Bin%d",iBin));

    TF1  *f_cosl = new TF1("f_cosl"
                          ,"([0]+[1]*x*0+[2]*0+[3]*0)+([4]+[5]*0+[6]*0+[7]*0)*x**2+([8]+[9]*0+[10]*0+[11]*0)*x**3+([12]+[13]*0+[14]*0+[15]*0)*x**4+([16]+[17]*0+[18]*0+[19]*0)*x**6"
                          ,-1.,1.);
    for (int iPar = 0; iPar < nPar; iPar++) f_cosl->SetParameter(iPar,arrPar[iPar]);
    f_cosl->Draw("SAME");
    canvas.Update();
    canvas.Print(TString::Format("./plots/recoEff_cosL_bin%d.pdf",iBin));

    // Clear
    delete f_cosl;
    delete h_cosl;
    delete f_cosk;
    delete h_cosk;
    delete latex;
    delete gMinuit;

    //prepare output
    std::vector<double> output;
    for (int iPar = 0; iPar < nPar; iPar++){
        output.push_back(arrPar[iPar]);
        output.push_back(arrParErr[iPar]);
    }
    for (int i = 0; i < output.size(); i=i+2) {
        printf("%f +- %f\n",output[i],output[i+1]);
    }
    return output;
}//}}}

//_________________________________________________________________________________
std::vector<double> angular2D_bin(int iBin, const char outfile[] = "angular2D")
{//{{{
    // Remark: You must use RooFit!! It's better in unbinned fit.
    //         Extended ML fitis adopted by Mauro, just follow!!
    
    // Tags
    bool is7TeVCheck = true;
    
    // Read data
    RooRealVar CosThetaK("CosThetaK", "cos#theta_{K}", -1., 1.);
    RooRealVar CosThetaL("CosThetaL", "cos#theta_{L}", -1., 1.);
    RooRealVar Q2("Q2","q^{2}",0.5,20.);
    RooRealVar fl("fl", "F_{L}", 0.8, -0.2, 1.2);
    RooRealVar afb("afb", "A_{FB}", 0., -1., 1.);
    RooRealVar fs("fs","F_{S}",0.0129254);//Derive from B0ToKstarJpsi
    fs.setAsymError(-0.00898344,0.0101371);
    //fs.setConstant(kTRUE);
    RooRealVar as("as","A_{S}",-0.0975919);//Derive from B0ToKstarJpsi
    as.setAsymError(-0.00490805,0.0049092);
    //as.setConstant(kTRUE);

    // Create physical distribution
    RooGenericPdf f_sig("f_sig", "9/16*((2/3*fs+4/3*as*CosThetaK)*(1-CosThetaL*CosThetaL)+(1-fs)*(2*fl*CosThetaK*CosThetaK*(1-CosThetaL*CosThetaL)+1/2*(1-fl)*(1-CosThetaK*CosThetaK)*(1+CosThetaL*CosThetaL)+4/3*afb*(1-CosThetaK*CosThetaK)*CosThetaL))", RooArgSet(CosThetaK,CosThetaL,fl,afb,fs,as));
    
    // Create acceptance*recoEfficiency map
    double *arrAccPar, *arrAccParErr, *arrRecPar, *arrRecParErr;
    if (is7TeVCheck){
        arrAccPar = arrAccPar2011[iBin];
        arrAccParErr = arrAccParErr2011[iBin];
        arrRecPar = arrRecPar2011[iBin];
        arrRecParErr = arrRecParErr2011[iBin];
    }else{
        arrAccPar = arrAccPar2012[iBin];
        arrAccParErr = arrAccParErr2012[iBin];
        arrRecPar = arrRecPar2012[iBin];
        arrRecParErr = arrRecParErr2012[iBin];
    }
    RooRealVar accK0L0("accK0L0","accK0L0",arrAccPar[ 0]);
    RooRealVar accK1L0("accK1L0","accK1L0",arrAccPar[ 1]);
    RooRealVar accK2L0("accK2L0","accK2L0",arrAccPar[ 2]);
    RooRealVar accK3L0("accK3L0","accK3L0",arrAccPar[ 3]);
    RooRealVar accK0L2("accK0L2","accK0L2",arrAccPar[ 4]);
    RooRealVar accK1L2("accK1L2","accK1L2",arrAccPar[ 5]);
    RooRealVar accK2L2("accK2L2","accK2L2",arrAccPar[ 6]);
    RooRealVar accK3L2("accK3L2","accK3L2",arrAccPar[ 7]);
    RooRealVar accK0L3("accK0L3","accK0L3",arrAccPar[ 8]);
    RooRealVar accK1L3("accK1L3","accK1L3",arrAccPar[ 9]);
    RooRealVar accK2L3("accK2L3","accK2L3",arrAccPar[10]);
    RooRealVar accK3L3("accK3L3","accK3L3",arrAccPar[11]);
    RooRealVar accK0L4("accK0L4","accK0L4",arrAccPar[12]);
    RooRealVar accK1L4("accK1L4","accK1L4",arrAccPar[13]);
    RooRealVar accK2L4("accK2L4","accK2L4",arrAccPar[14]);
    RooRealVar accK3L4("accK3L4","accK3L4",arrAccPar[15]);
    RooRealVar accK0L6("accK0L6","accK0L6",arrAccPar[16]);
    RooRealVar accK1L6("accK1L6","accK1L6",arrAccPar[17]);
    RooRealVar accK2L6("accK2L6","accK2L6",arrAccPar[18]);
    RooRealVar accK3L6("accK3L6","accK3L6",arrAccPar[19]);
    accK0L0.setError(arrAccParErr[ 0]);
    accK1L0.setError(arrAccParErr[ 1]);
    accK2L0.setError(arrAccParErr[ 2]);
    accK3L0.setError(arrAccParErr[ 3]);
    accK0L2.setError(arrAccParErr[ 4]);
    accK1L2.setError(arrAccParErr[ 5]);
    accK2L2.setError(arrAccParErr[ 6]);
    accK3L2.setError(arrAccParErr[ 7]);
    accK0L3.setError(arrAccParErr[ 8]);
    accK1L3.setError(arrAccParErr[ 9]);
    accK2L3.setError(arrAccParErr[10]);
    accK3L3.setError(arrAccParErr[11]);
    accK0L4.setError(arrAccParErr[12]);
    accK1L4.setError(arrAccParErr[13]);
    accK2L4.setError(arrAccParErr[14]);
    accK3L4.setError(arrAccParErr[15]);
    accK0L6.setError(arrAccParErr[16]);
    accK1L6.setError(arrAccParErr[17]);
    accK2L6.setError(arrAccParErr[18]);
    accK3L6.setError(arrAccParErr[19]);
    RooRealVar recK0L0("recK0L0","recK0L0",arrRecPar[ 0]);
    RooRealVar recK1L0("recK1L0","recK1L0",arrRecPar[ 1]);
    RooRealVar recK2L0("recK2L0","recK2L0",arrRecPar[ 2]);
    RooRealVar recK3L0("recK3L0","recK3L0",arrRecPar[ 3]);
    RooRealVar recK0L2("recK0L2","recK0L2",arrRecPar[ 4]);
    RooRealVar recK1L2("recK1L2","recK1L2",arrRecPar[ 5]);
    RooRealVar recK2L2("recK2L2","recK2L2",arrRecPar[ 6]);
    RooRealVar recK3L2("recK3L2","recK3L2",arrRecPar[ 7]);
    RooRealVar recK0L3("recK0L3","recK0L3",arrRecPar[ 8]);
    RooRealVar recK1L3("recK1L3","recK1L3",arrRecPar[ 9]);
    RooRealVar recK2L3("recK2L3","recK2L3",arrRecPar[10]);
    RooRealVar recK3L3("recK3L3","recK3L3",arrRecPar[11]);
    RooRealVar recK0L4("recK0L4","recK0L4",arrRecPar[12]);
    RooRealVar recK1L4("recK1L4","recK1L4",arrRecPar[13]);
    RooRealVar recK2L4("recK2L4","recK2L4",arrRecPar[14]);
    RooRealVar recK3L4("recK3L4","recK3L4",arrRecPar[15]);
    RooRealVar recK0L6("recK0L6","recK0L6",arrRecPar[16]);
    RooRealVar recK1L6("recK1L6","recK1L6",arrRecPar[17]);
    RooRealVar recK2L6("recK2L6","recK2L6",arrRecPar[18]);
    RooRealVar recK3L6("recK3L6","recK3L6",arrRecPar[19]);
    recK0L0.setError(arrRecParErr[ 0]);
    recK1L0.setError(arrRecParErr[ 1]);
    recK2L0.setError(arrRecParErr[ 2]);
    recK3L0.setError(arrRecParErr[ 3]);
    recK0L2.setError(arrRecParErr[ 4]);
    recK1L2.setError(arrRecParErr[ 5]);
    recK2L2.setError(arrRecParErr[ 6]);
    recK3L2.setError(arrRecParErr[ 7]);
    recK0L3.setError(arrRecParErr[ 8]);
    recK1L3.setError(arrRecParErr[ 9]);
    recK2L3.setError(arrRecParErr[10]);
    recK3L3.setError(arrRecParErr[11]);
    recK0L4.setError(arrRecParErr[12]);
    recK1L4.setError(arrRecParErr[13]);
    recK2L4.setError(arrRecParErr[14]);
    recK3L4.setError(arrRecParErr[15]);
    recK0L6.setError(arrRecParErr[16]);
    recK1L6.setError(arrRecParErr[17]);
    recK2L6.setError(arrRecParErr[18]);
    recK3L6.setError(arrRecParErr[19]);
    RooArgSet f_acc_argset(CosThetaL,CosThetaK);
    f_acc_argset.add(RooArgSet(accK0L0,accK1L0,accK2L0,accK3L0));
    f_acc_argset.add(RooArgSet(accK0L2,accK1L2,accK2L2,accK3L2));
    f_acc_argset.add(RooArgSet(accK0L3,accK1L3,accK2L3,accK3L3));
    f_acc_argset.add(RooArgSet(accK0L4,accK1L4,accK2L4,accK3L4));
    f_acc_argset.add(RooArgSet(accK0L6,accK1L6,accK2L6,accK3L6));
    RooArgSet f_rec_argset(CosThetaL,CosThetaK);
    f_rec_argset.add(RooArgSet(recK0L0,recK1L0,recK2L0,recK3L0));
    f_rec_argset.add(RooArgSet(recK0L2,recK1L2,recK2L2,recK3L2));
    f_rec_argset.add(RooArgSet(recK0L3,recK1L3,recK2L3,recK3L3));
    f_rec_argset.add(RooArgSet(recK0L4,recK1L4,recK2L4,recK3L4));
    f_rec_argset.add(RooArgSet(recK0L6,recK1L6,recK2L6,recK3L6));
    TString f_acc_format, f_acc_L0, f_acc_L2, f_acc_L3, f_acc_L4, f_acc_L6;
    TString f_rec_format, f_rec_L0, f_rec_L2, f_rec_L3, f_rec_L4, f_rec_L6;
    if (is7TeVCheck){
        f_acc_L0 = "(accK0L0+accK1L0*CosThetaK+accK2L0*CosThetaK**2+accK3L0*CosThetaK**3)";
        f_acc_L2 = "(accK0L2+accK1L2*CosThetaK+accK2L2*CosThetaK**2+accK3L2*CosThetaK**3)*CosThetaL**2";
        f_acc_L3 = "(accK0L3+accK1L3*CosThetaK+accK2L3*CosThetaK**2+accK3L3*CosThetaK**3)*CosThetaL**3";
        f_acc_L4 = "(accK0L4+accK1L4*CosThetaK+accK2L4*CosThetaK**2+accK3L4*CosThetaK**3)*CosThetaL**4";
        f_acc_L6 = "(accK0L6+accK1L6*CosThetaK+accK2L6*CosThetaK**2+accK3L6*CosThetaK**3)*CosThetaL**6";
        f_rec_L0 = "(recK0L0+recK1L0*CosThetaK+recK2L0*CosThetaK**2+recK3L0*CosThetaK**3)";
        f_rec_L2 = "(recK0L2+recK1L2*CosThetaK+recK2L2*CosThetaK**2+recK3L2*CosThetaK**3)*CosThetaL**2";
        f_rec_L3 = "(recK0L3+recK1L3*CosThetaK+recK2L3*CosThetaK**2+recK3L3*CosThetaK**3)*CosThetaL**3";
        f_rec_L4 = "(recK0L4+recK1L4*CosThetaK+recK2L4*CosThetaK**2+recK3L4*CosThetaK**3)*CosThetaL**4";
        f_rec_L6 = "(recK0L6+recK1L6*CosThetaK+recK2L6*CosThetaK**2+recK3L6*CosThetaK**3)*CosThetaL**6";
    }else{
        f_acc_L0 = "(accK0L0+accK1L0*CosThetaK+accK2L0*(3*CosThetaK**2-1)/2+accK3L0*(5*CosThetaK**3-3*CosThetaK)/2)";
        f_acc_L2 = "(accK0L2+accK1L2*CosThetaK+accK2L2*(3*CosThetaK**2-1)/2+accK3L3*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**2";
        f_acc_L3 = "(accK0L3+accK1L3*CosThetaK+accK2L3*(3*CosThetaK**2-1)/2+accK3L3*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**3";
        f_acc_L4 = "(accK0L4+accK1L4*CosThetaK+accK2L4*(3*CosThetaK**2-1)/2+accK3L3*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**4";
        f_acc_L6 = "(accK0L6+accK1L6*CosThetaK+accK2L6*(3*CosThetaK**2-1)/2+accK3L3*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**6";
        f_rec_L0 = "(recK0L0+recK1L0*CosThetaK+recK2L0*(3*CosThetaK**2-1)/2+recK3L0*(5*CosThetaK**3-3*CosThetaK)/2)";
        f_rec_L2 = "(recK0L2+recK1L2*CosThetaK+recK2L2*(3*CosThetaK**2-1)/2+recK3L3*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**2";
        f_rec_L3 = "(recK0L3+recK1L3*CosThetaK+recK2L3*(3*CosThetaK**2-1)/2+recK3L3*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**3";
        f_rec_L4 = "(recK0L4+recK1L4*CosThetaK+recK2L4*(3*CosThetaK**2-1)/2+recK3L3*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**4";
        f_rec_L6 = "(recK0L6+recK1L6*CosThetaK+recK2L6*(3*CosThetaK**2-1)/2+recK3L3*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**6";
    }

    if (iBin == 0) {
        f_acc_format = TString::Format("%s+%s+%s+%s",f_acc_L0.Data(),f_acc_L2.Data(),f_acc_L4.Data(),f_acc_L6.Data());
        f_rec_format = TString::Format("%s+%s+%s+%s",f_rec_L0.Data(),f_rec_L2.Data(),f_rec_L4.Data(),f_rec_L6.Data());
    }else if (iBin == 1) {
        f_acc_format = TString::Format("%s+%s+%s",f_acc_L0.Data(),f_acc_L2.Data(),f_acc_L4.Data());
        f_rec_format = TString::Format("%s+%s+%s",f_rec_L0.Data(),f_rec_L2.Data(),f_rec_L4.Data());
    }else if (iBin > 1 && iBin < 6) {
        f_acc_format = TString::Format("%s+%s+%s+%s",f_acc_L0.Data(),f_acc_L2.Data(),f_acc_L3.Data(),f_acc_L4.Data());
        f_rec_format = TString::Format("%s+%s+%s+%s",f_rec_L0.Data(),f_rec_L2.Data(),f_rec_L3.Data(),f_rec_L4.Data());
    }else{
        f_acc_format = TString::Format("%s+%s+%s",f_acc_L0.Data(),f_acc_L2.Data(),f_acc_L3.Data());
        f_rec_format = TString::Format("%s+%s+%s",f_rec_L0.Data(),f_rec_L2.Data(),f_rec_L3.Data());
    }
    RooGenericPdf f_acc("f_acc", f_acc_format,f_acc_argset);
    RooGenericPdf f_rec("f_rec", f_rec_format,f_rec_argset);
    RooProdPdf f_eff("f_eff","f_eff",f_acc,f_rec);

    // Observed spectrum = model*fullEfficiency
    RooRealVar nsig("nsig","nsig",1E4,1E2,1E8);
    RooProdPdf f("f","f",f_sig,f_eff);
    
    // Get data and apply unbinned fit
    RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(CosThetaK,CosThetaL,Q2),q2range[iBin],0);
    //RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"));
    RooFitResult *f_fitresult = f.fitTo(*data,Save(kTRUE),Minimizer("Minuit"));

    // Draw the frame on the canvas
    TCanvas* c = new TCanvas("c");
    RooPlot* framecosk = CosThetaK.frame(); 
    data->plotOn(framecosk,Binning(20)); 
    f.plotOn(framecosk); 
    //f.plotOn(framecosk,Components(f_sig),LineColor(4),LineWidth(2));
    //f.plotOn(framecosk,Components(f_bkg),LineStyle(2),LineColor(8),LineWidth(2));

    framecosk->SetTitle("");
    framecosk->SetMinimum(0);
    framecosk->Draw();

    TLatex *t1 = new TLatex();
    t1->SetNDC();
    double fixNDC = 0;
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    //t1->DrawLatex(.35,.80+fixNDC,TString::Format("F_{L}=%5.3f#pm%5.3f",fl.getVal(),fl.getError()));
    //t1->DrawLatex(.35,.74+fixNDC,TString::Format("A_{FB}=%5.3f#pm%5.3f",afb.getVal(),afb.getError()));
    if ( f_fitresult->status() == 0){
        //t1->DrawLatex(.35,.68+fixNDC,TString::Format("Fit status: %s","GOOD"));
        //t1->DrawLatex(.50,.68,TString::Format("Fit status: %d",f_fitresult->covQual()));
    }else{
        //t1->DrawLatex(.35,.68+fixNDC,TString::Format("Fit status: %s(%d)","BAD",f_fitresult->status()));
        //t1->DrawLatex(.50,.68,TString::Format("Fit status: %d",f_fitresult->covQual()));
    }
    //t1->DrawLatex(.35,.10,TString::Format("nbkg=%5.3f#pm%5.3f",nbkg.getVal(),nbkg.getError()));
    c->Print(TString::Format("./plots/%s_cosk_bin%d.pdf",outfile,iBin));

    // Draw projection to CosThetaL
    RooPlot* framecosl = CosThetaL.frame(); 
    data->plotOn(framecosl,Binning(20)); 
    f.plotOn(framecosl); 

    framecosl->SetTitle("");
    framecosl->SetMinimum(0);
    framecosl->Draw();

    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    c->Update();
    c->Print(TString::Format("./plots/%s_cosl_bin%d.pdf",outfile,iBin));

    // Make 2-D plot
    TH1 *h1 = data->createHistogram("CosThetaL,CosThetaK", 6, 5);
    h1->Draw("LEGO2");
    c->Update();
    c->Print(TString::Format("./plots/%s_bin%d.pdf",outfile,iBin));

    // clear
    delete t1;
    delete c;
    delete data;

    //write output
    std::vector<double> output;
    output.push_back(fl.getVal());
    output.push_back(fl.getError());
    output.push_back(afb.getVal());
    output.push_back(afb.getError());
    return output;
}//}}}

double readParam(int iBin, const char parName[], int iColumn)
{//{{{
    std::vector<double> output;
    char lineBuff[100];
    memset(lineBuff,' ',100*sizeof(char));
    FILE *fp = fopen(TString::Format("fitParameters%d.txt",iBin),"r");
    while(fgets(lineBuff,100,fp) != NULL ){
        if ( strstr(lineBuff,parName) != NULL ) break;
    }
    char *valBuff;
    valBuff = strtok(lineBuff," ");
    valBuff = strtok(NULL," ");
    while(valBuff != NULL){
        output.push_back(stof(valBuff));
        valBuff = strtok(NULL," ");
    }
    fclose(fp);
    return output.at(iColumn);
}//}}}

void angular3D_1a_Sm(int iBin, const char outfile[] = "angular3D_1a_Sm", bool keepParam = false)
{//{{{
    // Fit to signal simulation by YsSm+YcCm to determine Sm
    RooRealVar Bmass("Bmass","M_{K^{*}#Mu#Mu}",5.,5.56);
    RooRealVar Q2("Q2","q^{2}",0.5,20.);
    
    // Create parameters and PDFs
        // Signal double gaussian
    RooRealVar sigGauss_mean("sigGauss_mean","M_{K*#Mu#Mu}",5.28,5.25,5.30);
    RooRealVar sigGauss1_sigma("sigGauss1_sigma","#sigma_{1}",.03,.01,.08);
    RooRealVar sigGauss2_sigma("sigGauss2_sigma","#sigma_{2}",.15,.08,.35);
    RooRealVar sigM_frac("sigM_frac","sigM_frac",.5,0.,1.);
    
    // Create signal distribution
        // mass distro of signal
    RooGaussian f_sigMGauss1("f_sigMGauss1","f_sigMGauss1", Bmass, sigGauss_mean, sigGauss1_sigma);//double gaussian with shared mean
    RooGaussian f_sigMGauss2("f_sigMGauss2","f_sigMGauss2", Bmass, sigGauss_mean, sigGauss2_sigma);//double gaussian with shared mean
    RooAddPdf f_sigM("f_sigM","f_sigM", f_sigMGauss1, f_sigMGauss2, sigM_frac);
    
    // Create combinatorial background distribution
    RooRealVar bkgCombM_c("bkgCombM_c","c1",0,-30,50);
    RooRealVar offset("offset","offset",-5.);
    RooAddition Bmass_offset("Bmass_offset","Bmass_offset",RooArgList(Bmass,offset));
    RooExponential f_bkgCombM("f_bkgCombM","f_bkgCombM",Bmass_offset,bkgCombM_c);// exponential decay
    
    RooRealVar nsig("nsig","nsig",0,1E5);
    RooRealVar nbkg("nbkg","nbkg",0,1E5);
    RooAddPdf f("f", "f",RooArgList(f_sigM,f_bkgCombM),RooArgList(nsig,nbkg));

    // Get data and apply unbinned fit
    RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass),q2range[iBin],0);
    RooFitResult *f_fitresult = f.fitTo(*data,Save(kTRUE),Minimizer("Minuit"));

    // Draw the frame on the canvas
    TCanvas* c = new TCanvas("c");
    RooPlot* frame = Bmass.frame(); 
    data->plotOn(frame,Binning(20)); 
    f.plotOn(frame); 
    f.plotOn(frame,Components(f_sigM),LineColor(2),LineWidth(2));
    f.plotOn(frame,Components(f_bkgCombM),LineColor(3),LineStyle(2),LineWidth(2));

    frame->SetTitle("");
    frame->SetMinimum(0);
    frame->Draw();

    TLatex *t1 = new TLatex();
    t1->SetNDC();
    double fixNDC = 0;
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    //t1->DrawLatex(.35,.10,TString::Format("nbkg=%5.3f#pm%5.3f",nbkg.getVal(),nbkg.getError()));
    c->Print(TString::Format("./plots/%s_bin%d.pdf",outfile,iBin));

    // clear
    delete t1;
    delete c;
    delete data;

    // Prepare datacard
    if (keepParam){
        FILE *fp = fopen( TString::Format("fitParameters%d.txt",iBin) ,"a");
        fprintf(fp,"iBin %d\n",iBin);
        fprintf(fp,"mode 2011\n");
        fprintf(fp,"sigGauss1_sigma %f %f\n",sigGauss1_sigma.getVal(),sigGauss1_sigma.getError());
        fprintf(fp,"sigGauss2_sigma %f %f\n",sigGauss2_sigma.getVal(),sigGauss2_sigma.getError());
        fprintf(fp,"sigM_frac %f %f\n",sigM_frac.getVal(),sigM_frac.getError());
        fclose(fp);
    }
}//}}}
void angular3D_1b_YpPm(int iBin, const char outfile[] = "angular3D_1b_YpPm", bool keepParam = false)
{//{{{
    if (iBin ==0 || iBin%2 == 1){
        if (keepParam){
            FILE *fp = fopen( TString::Format("fitParameters%d.txt",iBin) ,"a");
            fprintf(fp,"bkgGauss1_mean 1.000000 0.000000\n");
            fprintf(fp,"bkgGauss1_sigma1 1.000000 0.000000\n");
            fprintf(fp,"bkgGauss1_sigma2 1.000000 0.000000\n");
            fprintf(fp,"bkgM_frac1 1.000000 0.000000\n");
            fprintf(fp,"bkgGauss2_mean 1.000000 0.000000\n");
            fprintf(fp,"bkgGauss2_sigma1 1.000000 0.000000\n");
            fprintf(fp,"bkgGauss2_sigma2 1.000000 0.000000\n");
            fprintf(fp,"bkgM_frac2 1.000000 0.000000\n");
            fprintf(fp,"nbkgPeak 0.000000 0.000000\n");
        }
        return;
    }

    // Fit to control channel simulations by YpPm to determine Yp,Pm.
    RooRealVar Bmass("Bmass","M_{K^{*}#Mu#Mu}",5.,5.56);
    RooRealVar Q2("Q2","q^{2}",0.5,20.);
    
    // Create peak background distribution
    RooRealVar bkgGauss1_mean("bkgGauss1_mean","M_{K*#Mu#Mu}",5.05,5.,5.12);
    RooRealVar bkgGauss2_mean("bkgGauss2_mean","M_{K*#Mu#Mu}",5.40,5.35,5.45);
    RooRealVar bkgGauss1_sigma1("bkgGauss1_sigma1","#sigma_{11}",.03,.01,.08);
    RooRealVar bkgGauss1_sigma2("bkgGauss1_sigma2","#sigma_{12}",.12,.08,.30);
    RooRealVar bkgGauss2_sigma1("bkgGauss2_sigma1","#sigma_{21}",.03,.01,.08);
    RooRealVar bkgGauss2_sigma2("bkgGauss2_sigma2","#sigma_{22}",.12,.08,.30);
    RooRealVar bkgM_frac1("bkgM_frac1","bkgM_frac1",1.,0.,1.);
    RooRealVar bkgM_frac2("bkgM_frac2","bkgM_frac2",1.,0.,1.);
    RooRealVar bkgM_frac12("bkgM_frac12","bkgM_frac12",0.,0.,1.);
    RooGaussian f_bkgPeakMGauss11("f_bkgPeakMGauss11","f_bkgPeakMGauss11", Bmass, bkgGauss1_mean, bkgGauss1_sigma1);
    RooGaussian f_bkgPeakMGauss12("f_bkgPeakMGauss12","f_bkgPeakMGauss12", Bmass, bkgGauss1_mean, bkgGauss1_sigma2);
    RooGaussian f_bkgPeakMGauss21("f_bkgPeakMGauss21","f_bkgPeakMGauss21", Bmass, bkgGauss2_mean, bkgGauss2_sigma1);
    RooGaussian f_bkgPeakMGauss22("f_bkgPeakMGauss22","f_bkgPeakMGauss22", Bmass, bkgGauss2_mean, bkgGauss2_sigma2);
    RooAddPdf f_bkgPeakM1("f_bkgPeakM1","f_bkgPeakM1", RooArgList(f_bkgPeakMGauss11, f_bkgPeakMGauss12), bkgM_frac1);
    RooAddPdf f_bkgPeakM2("f_bkgPeakM2","f_bkgPeakM2", RooArgList(f_bkgPeakMGauss21, f_bkgPeakMGauss22), bkgM_frac2);
    RooAddPdf f_bkgPeakM12("f_bkgPeakM12","f_bkgPeakM12", RooArgList(f_bkgPeakM1,f_bkgPeakM2), bkgM_frac12);
    
    RooRealVar nbkgPeak("nbkgPeak","nbkgPeak",1E2,1E1,1E7);
    RooExtendPdf *f = 0;
    switch (iBin) {
        case 2:
            //1 double guassian ,4+4 deg. ploy
            f = new RooExtendPdf("f","f",f_bkgPeakM1,nbkgPeak);
            break;
        case 4:
            //2 double guassian ,4+4 deg. ploy
            f = new RooExtendPdf("f","f",f_bkgPeakM12,nbkgPeak);
            break;
        case 6:
            //1 guassian ,2+2 deg. ploy
            f = new RooExtendPdf("f","f",f_bkgPeakMGauss21,nbkgPeak);
            break;
        default:
            break;
    }

    // Get data and apply unbinned fit
    RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass),q2range[iBin],0);
    RooFitResult *f_fitresult = f->fitTo(*data,Save(kTRUE),Minimizer("Minuit"),Extended());

    // Draw the frame on the canvas
    TCanvas* c = new TCanvas("c");
    RooPlot* frame = Bmass.frame(); 
    data->plotOn(frame,Binning(20)); 
    f->plotOn(frame); 

    frame->SetTitle("");
    frame->SetMinimum(0);
    frame->Draw();

    TLatex *t1 = new TLatex();
    t1->SetNDC();
    double fixNDC = 0;
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    c->Print(TString::Format("./plots/%s_bin%d.pdf",outfile,iBin));

    // clear
    delete t1;
    delete c;
    delete data;

    if (keepParam){
        FILE *fp = fopen( TString::Format("fitParameters%d.txt",iBin) ,"a");
        fprintf(fp,"bkgGauss1_mean %f %f\n",bkgGauss1_mean.getVal(),bkgGauss1_mean.getError());
        fprintf(fp,"bkgGauss1_sigma1 %f %f\n",bkgGauss1_sigma1.getVal(),bkgGauss1_sigma1.getError());
        fprintf(fp,"bkgGauss1_sigma2 %f %f\n",bkgGauss1_sigma2.getVal(),bkgGauss1_sigma2.getError());
        fprintf(fp,"bkgM_frac1 %f %f\n",bkgM_frac1.getVal(),bkgM_frac1.getError());
        fprintf(fp,"bkgGauss2_mean %f %f\n",bkgGauss2_mean.getVal(),bkgGauss2_mean.getError());
        fprintf(fp,"bkgGauss2_sigma1 %f %f\n",bkgGauss2_sigma1.getVal(),bkgGauss2_sigma1.getError());
        fprintf(fp,"bkgGauss2_sigma2 %f %f\n",bkgGauss2_sigma2.getVal(),bkgGauss2_sigma2.getError());
        fprintf(fp,"bkgM_frac2 %f %f\n",bkgM_frac2.getVal(),bkgM_frac2.getError());
        fprintf(fp,"nbkgPeak %f %f\n",nbkgPeak.getVal(),nbkgPeak.getError());
        fclose(fp);
    }
}//}}}
void angular3D_2a_PkPl(int iBin, const char outfile[] = "angular3D_2a_PkPl", bool keepParam = false)
{//{{{
    if (iBin ==0 || iBin%2 == 1){
        // Pm is flat(and the yield is 0) for bins other than 2,4,6
        if (keepParam){
            FILE *fp = fopen( TString::Format("fitParameters%d.txt",iBin) ,"a");
            fprintf(fp,"bkgPeakL_c1 0.000000 0.000000\n");
            fprintf(fp,"bkgPeakL_c2 0.000000 0.000000\n");
            fprintf(fp,"bkgPeakL_c3 0.000000 0.000000\n");
            fprintf(fp,"bkgPeakL_c4 0.000000 0.000000\n");
            fprintf(fp,"bkgPeakK_c1 0.000000 0.000000\n");
            fprintf(fp,"bkgPeakK_c2 0.000000 0.000000\n");
            fprintf(fp,"bkgPeakK_c3 0.000000 0.000000\n");
            fprintf(fp,"bkgPeakK_c4 0.000000 0.000000\n");
            fclose(fp);
        }
        return;
    }

    RooRealVar Q2("Q2","q^{2}",0.5,20.);
    RooRealVar CosThetaK("CosThetaK", "cos#theta_{K}", -1., 1.);
    RooRealVar CosThetaL("CosThetaL", "cos#theta_{L}", -1., 1.);
    
    RooArgSet f_bkgPeakL_argset;
    RooArgSet f_bkgPeakK_argset;
    RooRealVar bkgPeakL_c1("bkgPeakL_c1","c1",0,-5,5);
    RooRealVar bkgPeakL_c2("bkgPeakL_c2","c2",0,-5,5);
    RooRealVar bkgPeakL_c3("bkgPeakL_c3","c3",0,-5,5);
    RooRealVar bkgPeakL_c4("bkgPeakL_c4","c4",0,-5,5);
    RooRealVar bkgPeakK_c1("bkgPeakK_c1","c1",0,-5,5);
    RooRealVar bkgPeakK_c2("bkgPeakK_c2","c2",0,-5,5);
    RooRealVar bkgPeakK_c3("bkgPeakK_c3","c3",0,-5,5);
    RooRealVar bkgPeakK_c4("bkgPeakK_c4","c4",0,-5,5);
    switch (iBin) {
        case 2:
            //1 double guassian ,4+4 deg. ploy
            f_bkgPeakL_argset.add(RooArgSet(bkgPeakL_c1,bkgPeakL_c2,bkgPeakL_c3,bkgPeakL_c4));
            f_bkgPeakK_argset.add(RooArgSet(bkgPeakK_c1,bkgPeakK_c2,bkgPeakK_c3,bkgPeakK_c4));
            
            break;
        case 4:
            //2 double guassian ,4+4 deg. ploy
            f_bkgPeakL_argset.add(RooArgSet(bkgPeakL_c1,bkgPeakL_c2,bkgPeakL_c3,bkgPeakL_c4));
            f_bkgPeakK_argset.add(RooArgSet(bkgPeakK_c1,bkgPeakK_c2,bkgPeakK_c3,bkgPeakK_c4));
            break;
        case 6:
            //1 guassian ,2+2 deg. ploy
            f_bkgPeakL_argset.add(RooArgSet(bkgPeakL_c1,bkgPeakL_c2));
            f_bkgPeakK_argset.add(RooArgSet(bkgPeakK_c1,bkgPeakK_c2));
            break;
        default:
            break;
    }
    RooPolynomial f_bkgPeakL("f_bkgPeakL","f_bkgPeakL",CosThetaL,f_bkgPeakL_argset);
    RooPolynomial f_bkgPeakK("f_bkgPeakK","f_bkgPeakK",CosThetaK,f_bkgPeakK_argset);
    RooProdPdf f_bkgPeakA("f_bkgPeakA", "f_bckPeakA",f_bkgPeakK,f_bkgPeakL);
    
    // Get data
    RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(CosThetaK,CosThetaL,Q2),q2range[iBin],0);
    RooFitResult *f_fitresult = f_bkgPeakA.fitTo(*data,Save(kTRUE),Minimizer("Minuit"));

    // Draw CosThetaK
    TCanvas* c = new TCanvas("c");
    RooPlot* framecosk = CosThetaK.frame(); 
    data->plotOn(framecosk,Binning(20)); 
    f_bkgPeakK.plotOn(framecosk); 

    framecosk->SetTitle("");
    framecosk->SetMinimum(0);
    framecosk->Draw();

    TLatex *t1 = new TLatex();
    t1->SetNDC();
    double fixNDC = 0;
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    c->Print(TString::Format("./plots/%s_cosk_bin%d.pdf",outfile,iBin));
    
    // Draw CosThetaL
    RooPlot* framecosl = CosThetaL.frame(); 
    data->plotOn(framecosl,Binning(20)); 
    f_bkgPeakL.plotOn(framecosl); 

    framecosl->SetTitle("");
    framecosl->SetMinimum(0);
    framecosl->Draw();

    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    c->Update();
    c->Print(TString::Format("./plots/%s_cosl_bin%d.pdf",outfile,iBin));

    // Make 2-D plot
    TH1 *h1 = data->createHistogram("CosThetaL,CosThetaK", 6, 5);
    h1->Draw("LEGO2");
    c->Update();
    c->Print(TString::Format("./plots/%s_bin%d.pdf",outfile,iBin));

    // clear
    delete t1;
    delete c;
    delete data;

    if (keepParam){
        FILE *fp = fopen( TString::Format("fitParameters%d.txt",iBin) ,"a");
        fprintf(fp,"bkgPeakL_c1 %f %f\n",bkgPeakL_c1.getVal(),bkgPeakL_c1.getError());
        fprintf(fp,"bkgPeakL_c2 %f %f\n",bkgPeakL_c2.getVal(),bkgPeakL_c2.getError());
        fprintf(fp,"bkgPeakL_c3 %f %f\n",bkgPeakL_c3.getVal(),bkgPeakL_c3.getError());
        fprintf(fp,"bkgPeakL_c4 %f %f\n",bkgPeakL_c4.getVal(),bkgPeakL_c4.getError());
        fprintf(fp,"bkgPeakK_c1 %f %f\n",bkgPeakK_c1.getVal(),bkgPeakK_c1.getError());
        fprintf(fp,"bkgPeakK_c2 %f %f\n",bkgPeakK_c2.getVal(),bkgPeakK_c2.getError());
        fprintf(fp,"bkgPeakK_c3 %f %f\n",bkgPeakK_c3.getVal(),bkgPeakK_c3.getError());
        fprintf(fp,"bkgPeakK_c4 %f %f\n",bkgPeakK_c4.getVal(),bkgPeakK_c4.getError());
        fclose(fp);
    }
}//}}}
void angular3D_prior(int iBin, const char outfile[] = "angular3D_prior", bool keepParam = false)
{//{{{
    // Fit to signal simulation by YsSm+YcCm to determine Sm
    RooRealVar Bmass("Bmass","M_{K^{*}#Mu#Mu}",5.,5.56);
    RooRealVar Q2("Q2","q^{2}",0.5,20.);
    RooRealVar CosThetaK("CosThetaK", "cos#theta_{K}", -1., 1.);
    RooRealVar CosThetaL("CosThetaL", "cos#theta_{L}", -1., 1.);
    
    // Create combinatorial background distribution
    RooRealVar bkgCombL_c1("bkgCombL_c1","c1",0.,-1.5,1.5);
    RooRealVar bkgCombL_c2("bkgCombL_c2","c2",0.,-1.5,1.5);
    RooRealVar bkgCombL_c3("bkgCombL_c3","c3",0.,-1.5,1.5);
    RooRealVar bkgCombL_c4("bkgCombL_c4","c4",0.,-1.5,1.5);
    RooArgSet f_bkgCombL_argset;
    switch (iBin) {
        case 7:
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1));
            bkgCombL_c2.setConstant(kTRUE);
            bkgCombL_c3.setConstant(kTRUE);
            bkgCombL_c4.setConstant(kTRUE);
            break;
        case 0:
        case 1:
        case 4:
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1,bkgCombL_c2));
            bkgCombL_c3.setConstant(kTRUE);
            bkgCombL_c4.setConstant(kTRUE);
            break;
        case 2:
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1,bkgCombL_c2,bkgCombL_c3));
            bkgCombL_c4.setConstant(kTRUE);
            break;
        case 3:
        case 5:
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1,bkgCombL_c2,bkgCombL_c3,bkgCombL_c4));
            break;
        default:
            bkgCombL_c1.setConstant(kTRUE);
            bkgCombL_c2.setConstant(kTRUE);
            bkgCombL_c3.setConstant(kTRUE);
            bkgCombL_c4.setConstant(kTRUE);
            break;
    }
    RooPolynomial f_bkgCombL("f_bkgCombL","f_bkgCombL",CosThetaL,f_bkgCombL_argset);
    RooRealVar bkgCombK_c1("bkgCombK_c1","c1",0.,-1.5,1.5);
    RooRealVar bkgCombK_c2("bkgCombK_c2","c2",0.,-1.5,1.5);
    RooRealVar bkgCombK_c3("bkgCombK_c3","c3",0.,-1.5,1.5);
    RooRealVar bkgCombK_c4("bkgCombK_c4","c4",0.,-1.5,1.5);
    RooArgSet f_bkgCombK_argset;
    switch (iBin) {
        case 2:
            f_bkgCombK_argset.add(RooArgSet(bkgCombK_c1));
            bkgCombK_c2.setConstant(kTRUE);
            bkgCombK_c3.setConstant(kTRUE);
            bkgCombK_c4.setConstant(kTRUE);
            break;
        case 0:
        case 1:
            f_bkgCombK_argset.add(RooArgSet(bkgCombK_c1,bkgCombK_c2));
            bkgCombK_c3.setConstant(kTRUE);
            bkgCombK_c4.setConstant(kTRUE);
            break;
        case 3:
        case 4:
        case 5:
            f_bkgCombK_argset.add(RooArgSet(bkgCombK_c1,bkgCombK_c2,bkgCombK_c3,bkgCombK_c4));
            break;
        default:
            bkgCombK_c1.setConstant(kTRUE);
            bkgCombK_c2.setConstant(kTRUE);
            bkgCombK_c3.setConstant(kTRUE);
            bkgCombK_c4.setConstant(kTRUE);
            break;
    }
    RooPolynomial f_bkgCombK("f_bkgCombK","f_bkgCombK",CosThetaK,f_bkgCombK_argset);
    RooProdPdf f_bkgCombA("f_bkgCombA", "f_bckCombA",f_bkgCombK,f_bkgCombL);
    
    // Get data and apply unbinned fit
    //RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass),q2range[iBin],0);
    RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass,CosThetaK,CosThetaL),TString::Format("%s && (Bmass > 5.38 || Bmass < 5.18)",q2range[iBin]),0);
    RooFitResult *f_fitresult = f_bkgCombA.fitTo(*data,Save(kTRUE),Minimizer("Minuit"));

    // Draw the frame on the canvas
    TCanvas* c = new TCanvas("c");
    RooPlot* framecosk = CosThetaK.frame(); 
    data->plotOn(framecosk,Binning(20)); 
    f_bkgCombA.plotOn(framecosk); 

    framecosk->SetTitle("");
    framecosk->SetMinimum(0);
    framecosk->Draw();

    TLatex *t1 = new TLatex();
    t1->SetNDC();
    double fixNDC = 0;
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    c->Print(TString::Format("./plots/%s_cosk_bin%d.pdf",outfile,iBin));
    
    // 
    RooPlot* framecosl = CosThetaL.frame(); 
    data->plotOn(framecosl,Binning(20)); 
    f_bkgCombA.plotOn(framecosl); 

    framecosl->SetTitle("");
    framecosl->SetMinimum(0);
    framecosl->Draw();

    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    c->Update();
    c->Print(TString::Format("./plots/%s_cosl_bin%d.pdf",outfile,iBin));

    // clear
    delete t1;
    delete c;
    delete data;

    // Prepare datacard
    if (keepParam){
        FILE *fp = fopen( TString::Format("fitParameters%d.txt",iBin) ,"a");
        fprintf(fp,"bkgCombL_c1 %f %f\n",bkgCombL_c1.getVal(),bkgCombL_c1.getError());
        fprintf(fp,"bkgCombL_c2 %f %f\n",bkgCombL_c2.getVal(),bkgCombL_c2.getError());
        fprintf(fp,"bkgCombL_c3 %f %f\n",bkgCombL_c3.getVal(),bkgCombL_c3.getError());
        fprintf(fp,"bkgCombL_c4 %f %f\n",bkgCombL_c4.getVal(),bkgCombL_c4.getError());
        fprintf(fp,"bkgCombK_c1 %f %f\n",bkgCombK_c1.getVal(),bkgCombK_c1.getError());
        fprintf(fp,"bkgCombK_c2 %f %f\n",bkgCombK_c2.getVal(),bkgCombK_c2.getError());
        fprintf(fp,"bkgCombK_c3 %f %f\n",bkgCombK_c3.getVal(),bkgCombK_c3.getError());
        fprintf(fp,"bkgCombK_c4 %f %f\n",bkgCombK_c4.getVal(),bkgCombK_c4.getError());
        fclose(fp);
    }
}//}}}

std::vector<double> angular3D_bin(int iBin, const char outfile[] = "angular3D")
{//{{{
    // Remark: You must use RooFit!! It's better in unbinned fit.
    //         Extended ML fit is adopted by Mauro, just follow!!
    
    // Tags and Buff
    bool is7TeVCheck = true;
        // Efficiency map for 2011 result...
    // Read data
    RooRealVar CosThetaK("CosThetaK", "cos#theta_{K}", -1., 1.);
    RooRealVar CosThetaL("CosThetaL", "cos#theta_{L}", -1., 1.);
    RooRealVar Bmass("Bmass","M_{K^{*}#Mu#Mu}",5.,5.56);
    RooRealVar Q2("Q2","q^{2}",0.5,20.);

    // Create parameters and PDFs
        // Signal double gaussian
    RooRealVar sigGauss_mean("sigGauss_mean","M_{K*#Mu#Mu}",5.26,5.20,5.30);
    RooRealVar sigGauss1_sigma("sigGauss1_sigma","#sigma_{1}",readParam(iBin,"sigGauss1_sigma ",0));
    sigGauss_mean.setError(readParam(iBin,"sigGauss1_sigma",1));
    RooRealVar sigGauss2_sigma("sigGauss2_sigma","#sigma_{2}",readParam(iBin,"sigGauss2_sigma ",0));
    sigGauss_mean.setError(readParam(iBin,"sigGauss2_sigma ",1));
    RooRealVar sigM_frac("sigM_frac","sigM_frac",readParam(iBin,"sigM_frac ",0));
    sigGauss_mean.setError(readParam(iBin,"sigM_frac ",1));
        // Angular parameters
    RooRealVar afb("afb", "A_{FB}", 0., -1., 1.);
    RooRealVar fl("fl", "F_{L}", 0.8, -0.2, 1.2);
    RooRealVar fs("fs","F_{S}",0.01,0.005,0.2);//Derive from B0ToKstarJpsi, Bin3
    RooRealVar as("as","A_{S}",-0.1,-0.2,0.2);//Derive from B0ToKstarJpsi, Bin3
    if (iBin != 3 && iBin != 5){
        fs.setVal(0.0129254);
        fs.setAsymError(-0.00898344,0.0101371);
        as.setVal(-0.0975919);
        as.setAsymError(-0.00490805,0.0049092);
        //fs.setVal(readParam(3,"fs ",0));
        //fs.setAsymError(readParam(3,"fs ",1),readParam(3,"fs ",2));
        //as.setVal(readParam(3,"as ",0));
        //fs.setAsymError(readParam(3,"as ",1),readParam(3,"as ",2));
        fs.setConstant(kTRUE);
        as.setConstant(kTRUE);
    }
        // Efficiency and acceptance
    double *arrAccPar, *arrAccParErr, *arrRecPar, *arrRecParErr;
    if (is7TeVCheck){
        arrAccPar = arrAccPar2011[iBin];
        arrAccParErr = arrAccParErr2011[iBin];
        arrRecPar = arrRecPar2011[iBin];
        arrRecParErr = arrRecParErr2011[iBin];
    }else{
        arrAccPar = arrAccPar2012[iBin];
        arrAccParErr = arrAccParErr2012[iBin];
        arrRecPar = arrRecPar2012[iBin];
        arrRecParErr = arrRecParErr2012[iBin];
    }
    RooRealVar accK0L0("accK0L0","accK0L0",arrAccPar[ 0]);
    RooRealVar accK1L0("accK1L0","accK1L0",arrAccPar[ 1]);
    RooRealVar accK2L0("accK2L0","accK2L0",arrAccPar[ 2]);
    RooRealVar accK3L0("accK3L0","accK3L0",arrAccPar[ 3]);
    RooRealVar accK0L2("accK0L2","accK0L2",arrAccPar[ 4]);
    RooRealVar accK1L2("accK1L2","accK1L2",arrAccPar[ 5]);
    RooRealVar accK2L2("accK2L2","accK2L2",arrAccPar[ 6]);
    RooRealVar accK3L2("accK3L2","accK3L2",arrAccPar[ 7]);
    RooRealVar accK0L3("accK0L3","accK0L3",arrAccPar[ 8]);
    RooRealVar accK1L3("accK1L3","accK1L3",arrAccPar[ 9]);
    RooRealVar accK2L3("accK2L3","accK2L3",arrAccPar[10]);
    RooRealVar accK3L3("accK3L3","accK3L3",arrAccPar[11]);
    RooRealVar accK0L4("accK0L4","accK0L4",arrAccPar[12]);
    RooRealVar accK1L4("accK1L4","accK1L4",arrAccPar[13]);
    RooRealVar accK2L4("accK2L4","accK2L4",arrAccPar[14]);
    RooRealVar accK3L4("accK3L4","accK3L4",arrAccPar[15]);
    RooRealVar accK0L6("accK0L6","accK0L6",arrAccPar[16]);
    RooRealVar accK1L6("accK1L6","accK1L6",arrAccPar[17]);
    RooRealVar accK2L6("accK2L6","accK2L6",arrAccPar[18]);
    RooRealVar accK3L6("accK3L6","accK3L6",arrAccPar[19]);
    accK0L0.setError(arrAccParErr[ 0]);
    accK1L0.setError(arrAccParErr[ 1]);
    accK2L0.setError(arrAccParErr[ 2]);
    accK3L0.setError(arrAccParErr[ 3]);
    accK0L2.setError(arrAccParErr[ 4]);
    accK1L2.setError(arrAccParErr[ 5]);
    accK2L2.setError(arrAccParErr[ 6]);
    accK3L2.setError(arrAccParErr[ 7]);
    accK0L3.setError(arrAccParErr[ 8]);
    accK1L3.setError(arrAccParErr[ 9]);
    accK2L3.setError(arrAccParErr[10]);
    accK3L3.setError(arrAccParErr[11]);
    accK0L4.setError(arrAccParErr[12]);
    accK1L4.setError(arrAccParErr[13]);
    accK2L4.setError(arrAccParErr[14]);
    accK3L4.setError(arrAccParErr[15]);
    accK0L6.setError(arrAccParErr[16]);
    accK1L6.setError(arrAccParErr[17]);
    accK2L6.setError(arrAccParErr[18]);
    accK3L6.setError(arrAccParErr[19]);
    RooRealVar recK0L0("recK0L0","recK0L0",arrRecPar[ 0]);
    RooRealVar recK1L0("recK1L0","recK1L0",arrRecPar[ 1]);
    RooRealVar recK2L0("recK2L0","recK2L0",arrRecPar[ 2]);
    RooRealVar recK3L0("recK3L0","recK3L0",arrRecPar[ 3]);
    RooRealVar recK0L2("recK0L2","recK0L2",arrRecPar[ 4]);
    RooRealVar recK1L2("recK1L2","recK1L2",arrRecPar[ 5]);
    RooRealVar recK2L2("recK2L2","recK2L2",arrRecPar[ 6]);
    RooRealVar recK3L2("recK3L2","recK3L2",arrRecPar[ 7]);
    RooRealVar recK0L3("recK0L3","recK0L3",arrRecPar[ 8]);
    RooRealVar recK1L3("recK1L3","recK1L3",arrRecPar[ 9]);
    RooRealVar recK2L3("recK2L3","recK2L3",arrRecPar[10]);
    RooRealVar recK3L3("recK3L3","recK3L3",arrRecPar[11]);
    RooRealVar recK0L4("recK0L4","recK0L4",arrRecPar[12]);
    RooRealVar recK1L4("recK1L4","recK1L4",arrRecPar[13]);
    RooRealVar recK2L4("recK2L4","recK2L4",arrRecPar[14]);
    RooRealVar recK3L4("recK3L4","recK3L4",arrRecPar[15]);
    RooRealVar recK0L6("recK0L6","recK0L6",arrRecPar[16]);
    RooRealVar recK1L6("recK1L6","recK1L6",arrRecPar[17]);
    RooRealVar recK2L6("recK2L6","recK2L6",arrRecPar[18]);
    RooRealVar recK3L6("recK3L6","recK3L6",arrRecPar[19]);
    recK0L0.setError(arrRecParErr[ 0]);
    recK1L0.setError(arrRecParErr[ 1]);
    recK2L0.setError(arrRecParErr[ 2]);
    recK3L0.setError(arrRecParErr[ 3]);
    recK0L2.setError(arrRecParErr[ 4]);
    recK1L2.setError(arrRecParErr[ 5]);
    recK2L2.setError(arrRecParErr[ 6]);
    recK3L2.setError(arrRecParErr[ 7]);
    recK0L3.setError(arrRecParErr[ 8]);
    recK1L3.setError(arrRecParErr[ 9]);
    recK2L3.setError(arrRecParErr[10]);
    recK3L3.setError(arrRecParErr[11]);
    recK0L4.setError(arrRecParErr[12]);
    recK1L4.setError(arrRecParErr[13]);
    recK2L4.setError(arrRecParErr[14]);
    recK3L4.setError(arrRecParErr[15]);
    recK0L6.setError(arrRecParErr[16]);
    recK1L6.setError(arrRecParErr[17]);
    recK2L6.setError(arrRecParErr[18]);
    recK3L6.setError(arrRecParErr[19]);
    RooArgSet f_acc_argset(CosThetaL,CosThetaK);
    f_acc_argset.add(RooArgSet(accK0L0,accK1L0,accK2L0,accK3L0));
    f_acc_argset.add(RooArgSet(accK0L2,accK1L2,accK2L2,accK3L2));
    f_acc_argset.add(RooArgSet(accK0L3,accK1L3,accK2L3,accK3L3));
    f_acc_argset.add(RooArgSet(accK0L4,accK1L4,accK2L4,accK3L4));
    f_acc_argset.add(RooArgSet(accK0L6,accK1L6,accK2L6,accK3L6));
    RooArgSet f_rec_argset(CosThetaL,CosThetaK);
    f_rec_argset.add(RooArgSet(recK0L0,recK1L0,recK2L0,recK3L0));
    f_rec_argset.add(RooArgSet(recK0L2,recK1L2,recK2L2,recK3L2));
    f_rec_argset.add(RooArgSet(recK0L3,recK1L3,recK2L3,recK3L3));
    f_rec_argset.add(RooArgSet(recK0L4,recK1L4,recK2L4,recK3L4));
    f_rec_argset.add(RooArgSet(recK0L6,recK1L6,recK2L6,recK3L6));
    TString f_acc_format, f_acc_L0, f_acc_L2, f_acc_L3, f_acc_L4, f_acc_L6;
    TString f_rec_format, f_rec_L0, f_rec_L2, f_rec_L3, f_rec_L4, f_rec_L6;
    if (is7TeVCheck){
        f_acc_L0 = "(accK0L0+accK1L0*CosThetaK+accK2L0*CosThetaK**2+accK3L0*CosThetaK**3)";
        f_acc_L2 = "(accK0L2+accK1L2*CosThetaK+accK2L2*CosThetaK**2+accK3L2*CosThetaK**3)*CosThetaL**2";
        f_acc_L3 = "(accK0L3+accK1L3*CosThetaK+accK2L3*CosThetaK**2+accK3L3*CosThetaK**3)*CosThetaL**3";
        f_acc_L4 = "(accK0L4+accK1L4*CosThetaK+accK2L4*CosThetaK**2+accK3L4*CosThetaK**3)*CosThetaL**4";
        f_acc_L6 = "(accK0L6+accK1L6*CosThetaK+accK2L6*CosThetaK**2+accK3L6*CosThetaK**3)*CosThetaL**6";
        f_rec_L0 = "(recK0L0+recK1L0*CosThetaK+recK2L0*CosThetaK**2+recK3L0*CosThetaK**3)";
        f_rec_L2 = "(recK0L2+recK1L2*CosThetaK+recK2L2*CosThetaK**2+recK3L2*CosThetaK**3)*CosThetaL**2";
        f_rec_L3 = "(recK0L3+recK1L3*CosThetaK+recK2L3*CosThetaK**2+recK3L3*CosThetaK**3)*CosThetaL**3";
        f_rec_L4 = "(recK0L4+recK1L4*CosThetaK+recK2L4*CosThetaK**2+recK3L4*CosThetaK**3)*CosThetaL**4";
        f_rec_L6 = "(recK0L6+recK1L6*CosThetaK+recK2L6*CosThetaK**2+recK3L6*CosThetaK**3)*CosThetaL**6";
    }else{
        f_acc_L0 = "(accK0L0+accK1L0*CosThetaK+accK2L0*(3*CosThetaK**2-1)/2+accK3L0*(5*CosThetaK**3-3*CosThetaK)/2)";
        f_acc_L2 = "(accK0L2+accK1L2*CosThetaK+accK2L2*(3*CosThetaK**2-1)/2+accK3L3*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**2";
        f_acc_L3 = "(accK0L3+accK1L3*CosThetaK+accK2L3*(3*CosThetaK**2-1)/2+accK3L3*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**3";
        f_acc_L4 = "(accK0L4+accK1L4*CosThetaK+accK2L4*(3*CosThetaK**2-1)/2+accK3L3*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**4";
        f_acc_L6 = "(accK0L6+accK1L6*CosThetaK+accK2L6*(3*CosThetaK**2-1)/2+accK3L3*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**6";
        f_rec_L0 = "(recK0L0+recK1L0*CosThetaK+recK2L0*(3*CosThetaK**2-1)/2+recK3L0*(5*CosThetaK**3-3*CosThetaK)/2)";
        f_rec_L2 = "(recK0L2+recK1L2*CosThetaK+recK2L2*(3*CosThetaK**2-1)/2+recK3L3*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**2";
        f_rec_L3 = "(recK0L3+recK1L3*CosThetaK+recK2L3*(3*CosThetaK**2-1)/2+recK3L3*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**3";
        f_rec_L4 = "(recK0L4+recK1L4*CosThetaK+recK2L4*(3*CosThetaK**2-1)/2+recK3L3*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**4";
        f_rec_L6 = "(recK0L6+recK1L6*CosThetaK+recK2L6*(3*CosThetaK**2-1)/2+recK3L3*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**6";
    }

    if (iBin == 0) {
        f_acc_format = TString::Format("%s+%s+%s+%s",f_acc_L0.Data(),f_acc_L2.Data(),f_acc_L4.Data(),f_acc_L6.Data());
        f_rec_format = TString::Format("%s+%s+%s+%s",f_rec_L0.Data(),f_rec_L2.Data(),f_rec_L4.Data(),f_rec_L6.Data());
    }else if (iBin == 1) {
        f_acc_format = TString::Format("%s+%s+%s",f_acc_L0.Data(),f_acc_L2.Data(),f_acc_L4.Data());
        f_rec_format = TString::Format("%s+%s+%s",f_rec_L0.Data(),f_rec_L2.Data(),f_rec_L4.Data());
    }else if (iBin > 1 && iBin < 6) {
        f_acc_format = TString::Format("%s+%s+%s+%s",f_acc_L0.Data(),f_acc_L2.Data(),f_acc_L3.Data(),f_acc_L4.Data());
        f_rec_format = TString::Format("%s+%s+%s+%s",f_rec_L0.Data(),f_rec_L2.Data(),f_rec_L3.Data(),f_rec_L4.Data());
    }else{
        f_acc_format = TString::Format("%s+%s+%s",f_acc_L0.Data(),f_acc_L2.Data(),f_acc_L3.Data());
        f_rec_format = TString::Format("%s+%s+%s",f_rec_L0.Data(),f_rec_L2.Data(),f_rec_L3.Data());
    }

    // Create signal distribution
        // mass distro of signal
    RooGaussian f_sigMGauss1("f_sigMGauss1","f_sigMGauss1", Bmass, sigGauss_mean, sigGauss1_sigma);//double gaussian with shared mean
    RooGaussian f_sigMGauss2("f_sigMGauss2","f_sigMGauss2", Bmass, sigGauss_mean, sigGauss2_sigma);//double gaussian with shared mean
    RooAddPdf f_sigM("f_sigM","f_sigM", RooArgList(f_sigMGauss1, f_sigMGauss2), sigM_frac);
        // angular distro of signal
    RooGenericPdf f_sigA("f_sigA", "9/16*((2/3*fs+4/3*as*CosThetaK)*(1-CosThetaL*CosThetaL)+(1-fs)*(2*fl*CosThetaK*CosThetaK*(1-CosThetaL*CosThetaL)+1/2*(1-fl)*(1-CosThetaK*CosThetaK)*(1+CosThetaL*CosThetaL)+4/3*afb*(1-CosThetaK*CosThetaK)*CosThetaL))", RooArgSet(CosThetaK,CosThetaL,fl,afb,fs,as));
        // efficiency map of signal
    RooGenericPdf f_acc("f_acc", f_acc_format,f_acc_argset);
    RooGenericPdf f_rec("f_rec", f_rec_format,f_rec_argset);
    RooProdPdf f_eff("f_eff","f_eff",f_acc,f_rec);
        // Mass*Angular*efficiency
    RooProdPdf f_sigMA("f_sigMA","f_sigMA",f_sigM,f_sigA);
    RooProdPdf f_sig("f_sig","f_sig",f_sigMA,f_eff);
    printf("INFO: f_sig prepared.\n");

    // Create combinatorial background distribution
    RooRealVar bkgCombM_c("bkgCombM_c","c1",0,-10,10);
    RooRealVar offset("offset","offset",-5.);
    RooAddition Bmass_offset("Bmass_offset","Bmass_offset",RooArgList(Bmass,offset));
    RooExponential f_bkgCombM("f_bkgCombM","f_bkgCombM",Bmass_offset,bkgCombM_c);// exponential decay
    RooRealVar bkgCombL_c1("bkgCombL_c1","c1",readParam(iBin,"bkgCombL_c1",0),-1.5,1.5);
    RooRealVar bkgCombL_c2("bkgCombL_c2","c2",readParam(iBin,"bkgCombL_c2",0),-1.5,1.5);
    RooRealVar bkgCombL_c3("bkgCombL_c3","c3",readParam(iBin,"bkgCombL_c3",0),-1.5,1.5);
    RooRealVar bkgCombL_c4("bkgCombL_c4","c4",readParam(iBin,"bkgCombL_c4",0),-1.5,1.5);
    RooArgSet f_bkgCombL_argset;
    switch (iBin) {
        case 7:
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1));
            bkgCombL_c2.setConstant(kTRUE);
            bkgCombL_c3.setConstant(kTRUE);
            bkgCombL_c4.setConstant(kTRUE);
            break;
        case 0:
        case 1:
        case 4:
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1,bkgCombL_c2));
            bkgCombL_c3.setConstant(kTRUE);
            bkgCombL_c4.setConstant(kTRUE);
            break;
        case 2:
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1,bkgCombL_c2,bkgCombL_c3));
            bkgCombL_c4.setConstant(kTRUE);
            break;
        case 3:
        case 5:
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1,bkgCombL_c2,bkgCombL_c3,bkgCombL_c4));
            break;
        default:
            bkgCombL_c1.setConstant(kTRUE);
            bkgCombL_c2.setConstant(kTRUE);
            bkgCombL_c3.setConstant(kTRUE);
            bkgCombL_c4.setConstant(kTRUE);
            break;
    }
    RooPolynomial f_bkgCombL("f_bkgCombL","f_bkgCombL",CosThetaL,f_bkgCombL_argset);
    RooRealVar bkgCombK_c1("bkgCombK_c1","c1",readParam(iBin,"bkgCombK_c1",0),-1.5,1.5);
    RooRealVar bkgCombK_c2("bkgCombK_c2","c2",readParam(iBin,"bkgCombK_c2",0),-1.5,1.5);
    RooRealVar bkgCombK_c3("bkgCombK_c3","c3",readParam(iBin,"bkgCombK_c3",0),-1.5,1.5);
    RooRealVar bkgCombK_c4("bkgCombK_c4","c4",readParam(iBin,"bkgCombK_c4",0),-1.5,1.5);
    RooArgSet f_bkgCombK_argset;
    switch (iBin) {
        case 2:
            f_bkgCombK_argset.add(RooArgSet(bkgCombK_c1));
            bkgCombK_c2.setConstant(kTRUE);
            bkgCombK_c3.setConstant(kTRUE);
            bkgCombK_c4.setConstant(kTRUE);
            break;
        case 0:
        case 1:
            f_bkgCombK_argset.add(RooArgSet(bkgCombK_c1,bkgCombK_c2));
            bkgCombK_c3.setConstant(kTRUE);
            bkgCombK_c4.setConstant(kTRUE);
            break;
        case 3:
        case 4:
        case 5:
            f_bkgCombK_argset.add(RooArgSet(bkgCombK_c1,bkgCombK_c2,bkgCombK_c3,bkgCombK_c4));
            break;
        default:
            bkgCombK_c1.setConstant(kTRUE);
            bkgCombK_c2.setConstant(kTRUE);
            bkgCombK_c3.setConstant(kTRUE);
            bkgCombK_c4.setConstant(kTRUE);
            break;
    }
    RooPolynomial f_bkgCombK("f_bkgCombK","f_bkgCombK",CosThetaK,f_bkgCombK_argset);
    RooProdPdf f_bkgCombA("f_bkgCombA", "f_bckCombA",f_bkgCombK,f_bkgCombL);
    RooProdPdf f_bkgComb("f_bkgComb", "f_bckComb",f_bkgCombA,f_bkgCombM);
    printf("INFO: f_bkgComb prepared.\n");
    
    // Create peak background distribution
    RooRealVar bkgGauss1_mean("bkgGauss1_mean","M_{K*#Mu#Mu}",readParam(iBin,"bkgGauss1_mean",0));
    bkgGauss1_mean.setError(readParam(iBin,"bkgGauss1_mean",1));
    RooRealVar bkgGauss1_sigma1("bkgGauss1_sigma1","#sigma_{11}",readParam(iBin,"bkgGauss1_sigma1",0));
    bkgGauss1_sigma1.setError(readParam(iBin,"bkgGauss1_sigma1",1));
    RooRealVar bkgGauss1_sigma2("bkgGauss1_sigma2","#sigma_{12}",readParam(iBin,"bkgGauss1_sigma2",0));
    bkgGauss1_sigma2.setError(readParam(iBin,"bkgGauss1_sigma2",1));
    RooRealVar bkgM_frac1("bkgM_frac1","bkgM_frac1",readParam(iBin,"bkgM_frac1",0));
    bkgM_frac1.setError(readParam(iBin,"bkgM_frac1",1));
    RooRealVar bkgGauss2_mean("bkgGauss2_mean","M_{K*#Mu#Mu}",readParam(iBin,"bkgGauss2_mean",0));
    bkgGauss2_mean.setError(readParam(iBin,"bkgGauss1_mean",1));
    RooRealVar bkgGauss2_sigma1("bkgGauss2_sigma1","#sigma_{21}",readParam(iBin,"bkgGauss2_sigma1",0));
    bkgGauss2_sigma1.setError(readParam(iBin,"bkgGauss1_sigma1",1));
    RooRealVar bkgGauss2_sigma2("bkgGauss2_sigma2","#sigma_{22}",readParam(iBin,"bkgGauss2_sigma2",0));
    bkgGauss2_sigma2.setError(readParam(iBin,"bkgGauss2_sigma2",1));
    RooRealVar bkgM_frac2("bkgM_frac2","bkgM_frac2",readParam(iBin,"bkgM_frac2",0));
    bkgM_frac2.setError(readParam(iBin,"bkgM_frac2",1));
    RooRealVar bkgM_frac12("bkgM_frac12","bkgM_frac12",readParam(iBin,"bkgM_frac12",0));
    bkgM_frac12.setError(readParam(iBin,"bkgM_frac12",1));
    RooGaussian f_bkgPeakMGauss11("f_bkgPeakMGauss11","f_bkgPeakMGauss11", Bmass, bkgGauss1_mean, bkgGauss1_sigma1);
    RooGaussian f_bkgPeakMGauss12("f_bkgPeakMGauss12","f_bkgPeakMGauss12", Bmass, bkgGauss1_mean, bkgGauss1_sigma2);
    RooGaussian f_bkgPeakMGauss21("f_bkgPeakMGauss21","f_bkgPeakMGauss21", Bmass, bkgGauss2_mean, bkgGauss2_sigma1);
    RooGaussian f_bkgPeakMGauss22("f_bkgPeakMGauss22","f_bkgPeakMGauss22", Bmass, bkgGauss2_mean, bkgGauss2_sigma2);
    RooAddPdf f_bkgPeakM1("f_bkgPeakM1","f_bkgPeakM1", RooArgList(f_bkgPeakMGauss11, f_bkgPeakMGauss12), bkgM_frac1);
    RooAddPdf f_bkgPeakM2("f_bkgPeakM2","f_bkgPeakM2", RooArgList(f_bkgPeakMGauss21, f_bkgPeakMGauss22), bkgM_frac2);
    RooAddPdf f_bkgPeakM12("f_bkgPeakM12","f_bkgPeakM12", RooArgList(f_bkgPeakM1, f_bkgPeakM2), bkgM_frac12);
    RooRealVar bkgPeakL_c1("bkgPeakL_c1","c1",readParam(iBin,"bkgPeakL_c1",0));
    bkgPeakL_c1.setError(readParam(iBin,"bkgPeakL_c1",1));
    RooRealVar bkgPeakL_c2("bkgPeakL_c2","c2",readParam(iBin,"bkgPeakL_c2",0));
    bkgPeakL_c2.setError(readParam(iBin,"bkgPeakL_c2",1));
    RooRealVar bkgPeakL_c3("bkgPeakL_c3","c3",readParam(iBin,"bkgPeakL_c3",0));
    bkgPeakL_c3.setError(readParam(iBin,"bkgPeakL_c3",1));
    RooRealVar bkgPeakL_c4("bkgPeakL_c4","c4",readParam(iBin,"bkgPeakL_c4",0));
    bkgPeakL_c4.setError(readParam(iBin,"bkgPeakL_c4",1));
    RooRealVar bkgPeakK_c1("bkgPeakK_c1","c1",readParam(iBin,"bkgPeakK_c1",0));
    bkgPeakL_c1.setError(readParam(iBin,"bkgPeakK_c1",1));
    RooRealVar bkgPeakK_c2("bkgPeakK_c2","c2",readParam(iBin,"bkgPeakK_c2",0));
    bkgPeakL_c2.setError(readParam(iBin,"bkgPeakK_c2",1));
    RooRealVar bkgPeakK_c3("bkgPeakK_c3","c3",readParam(iBin,"bkgPeakK_c3",0));
    bkgPeakL_c3.setError(readParam(iBin,"bkgPeakK_c3",1));
    RooRealVar bkgPeakK_c4("bkgPeakK_c4","c4",readParam(iBin,"bkgPeakK_c4",0));
    bkgPeakL_c4.setError(readParam(iBin,"bkgPeakK_c4",1));
    RooArgSet f_bkgPeakL_argset(bkgPeakL_c1,bkgPeakL_c2,bkgPeakL_c3,bkgPeakL_c4);
    RooArgSet f_bkgPeakK_argset(bkgPeakK_c1,bkgPeakK_c2,bkgPeakK_c3,bkgPeakK_c4);
    switch (iBin) {// Should be fixed constants already.
        case 2:
            //1 double guassian ,4+4 deg. ploy
            bkgM_frac12.setVal(1.);
            bkgM_frac12.setConstant(kTRUE);
            break;
        case 4:
            //2 double guassian ,4+4 deg. ploy
            bkgM_frac12.setVal(0.);
            bkgM_frac1.setVal(1.);
            bkgM_frac12.setConstant(kTRUE);
            bkgM_frac2.setConstant(kTRUE);
            break;
        case 6:
            //1 guassian ,2+2 deg. ploy
            bkgPeakK_c3.setConstant(kTRUE);
            bkgPeakL_c3.setConstant(kTRUE);
            bkgPeakK_c4.setConstant(kTRUE);
            bkgPeakL_c4.setConstant(kTRUE);
            bkgM_frac12.setVal(0.);
            bkgM_frac1.setVal(1.);
            bkgM_frac12.setConstant(kTRUE);
            bkgM_frac2.setConstant(kTRUE);
            break;
        default:
            bkgPeakK_c1.setConstant(kTRUE);
            bkgPeakL_c1.setConstant(kTRUE);
            bkgPeakK_c2.setConstant(kTRUE);
            bkgPeakL_c2.setConstant(kTRUE);
            bkgPeakK_c3.setConstant(kTRUE);
            bkgPeakL_c3.setConstant(kTRUE);
            bkgPeakK_c4.setConstant(kTRUE);
            bkgPeakL_c4.setConstant(kTRUE);
            break;
    }
    RooPolynomial f_bkgPeakL("f_bkgPeakL","f_bkgPeakL",Bmass,f_bkgPeakL_argset);
    RooPolynomial f_bkgPeakK("f_bkgPeakK","f_bkgPeakK",Bmass,f_bkgPeakK_argset);
    RooProdPdf f_bkgPeakA("f_bkgPeakA", "f_bckPeakA",f_bkgPeakK,f_bkgPeakL);
    RooProdPdf f_bkgPeak("f_bkgPeak", "f_bkgPeak",f_bkgPeakA,f_bkgPeakM12);
    printf("INFO: f_bkgPeak prepared.\n");

    // Observed spectrum = model*fullEfficiency
    RooRealVar nsig("nsig","nsig",1E4,1E2,1E8);
    RooRealVar nbkgComb("nbkgComb","nbkgComb",10,0.1,1E3);
    RooRealVar nbkgPeak("nbkgPeak","nbkgPeak",readParam(iBin,"nbkgPeak",0));
    nbkgPeak.setError(readParam(iBin,"nbkgPeak",1));
    RooAddPdf f_model("f_model","f_model",RooArgList(f_sig,f_bkgComb,f_bkgPeak),RooArgList(nsig,nbkgComb,nbkgPeak));

    // Extra penalty term to confine As, Fs, Fl, Afb.
    RooRealVar t_penalty("t_penalty","t",0.01);
    RooGenericPdf f_penalty("f_penalty","(1-TMath::Erf((afb-0.75*(1-fl))/(1.5*t_penalty*(1-fl))))*(1-TMath::Erf((-afb-0.75*(1-fl))/(1.5*t_penalty*(1-fl))))",RooArgSet(afb,fl,t_penalty));
    RooProdPdf f("f","f",f_model,f_penalty);
    printf("INFO: f_penalty prepared.\n");
    
    // Get data and apply unbinned fit
    RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Bmass,CosThetaK,CosThetaL,Q2),q2range[iBin],0);
    RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"));

    // Draw the frame on the canvas
    TCanvas* c = new TCanvas("c");
    RooPlot* framemass = Bmass.frame();
    data->plotOn(framemass,Binning(20));
    f.plotOn(framemass);
    f.plotOn(framemass,Components(f_sig),LineColor(4),LineWidth(2));
    f.plotOn(framemass,Components(f_bkgPeak),LineColor(3),LineWidth(2),LineStyle(2));
    f.plotOn(framemass,Components(f_bkgComb),LineColor(8),LineWidth(2),LineStyle(2));

    framemass->SetTitle("");
    framemass->SetMinimum(0);
    framemass->Draw();
    
    TLatex *t1 = new TLatex();
    t1->SetNDC();
    double fixNDC = 0;
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    //t1->DrawLatex(.35,.80+fixNDC,TString::Format("F_{L}=%5.3f#pm%5.3f",fl.getVal(),fl.getError()));
    //t1->DrawLatex(.35,.74+fixNDC,TString::Format("A_{FB}=%5.3f#pm%5.3f",afb.getVal(),afb.getError()));
    if ( f_fitresult->status() == 0){
        //t1->DrawLatex(.35,.68+fixNDC,TString::Format("Fit status: %s","GOOD"));
        //t1->DrawLatex(.50,.68,TString::Format("Fit status: %d",f_fitresult->covQual()));
    }else{
        //t1->DrawLatex(.35,.68+fixNDC,TString::Format("Fit status: %s(%d)","BAD",f_fitresult->status()));
        //t1->DrawLatex(.50,.68,TString::Format("Fit status: %d",f_fitresult->covQual()));
    }
    //t1->DrawLatex(.35,.10,TString::Format("nbkg=%5.3f#pm%5.3f",nbkg.getVal(),nbkg.getError()));
    c->Print(TString::Format("./plots/%s_bin%d.pdf",outfile,iBin));

    // Draw projection to CosThetaK
    RooPlot* framecosk = CosThetaK.frame(); 
    data->plotOn(framecosk,Binning(20)); 
    f.plotOn(framecosk); 

    framecosk->SetTitle("");
    framecosk->SetMinimum(0);
    framecosk->Draw();

    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    c->Update();
    c->Print(TString::Format("./plots/%s_cosk_bin%d.pdf",outfile,iBin));

    // Draw projection to CosThetaL
    RooPlot* framecosl = CosThetaL.frame(); 
    data->plotOn(framecosl,Binning(20)); 
    f.plotOn(framecosl); 

    framecosl->SetTitle("");
    framecosl->SetMinimum(0);
    framecosl->Draw();

    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    c->Update();
    c->Print(TString::Format("./plots/%s_cosl_bin%d.pdf",outfile,iBin));

    // clear
    delete t1;
    delete c;
    delete data;

    // write output
    //FILE *fp = fopen( TString::Format("fitParameters%d.txt",iBin) ,"a");
    //fprintf(fp,"fl %f %f\n" ,fl.getVal(),fl.getError());
    //fprintf(fp,"afb %f %f\n",afb.getVal(),afb.getError());
    //fprintf(fp,"fs %f %f %f\n" ,fs.getVal(),fs.getErrorHi(),fs.getErrorLo());
    //fprintf(fp,"as %f %f %f\n" ,as.getVal(),as.getErrorHi(),fs.getErrorLo());
    //fclose(fp);

    std::vector<double> output;
    output.push_back(fl.getVal());
    output.push_back(fl.getError());
    output.push_back(afb.getVal());
    output.push_back(afb.getError());
    //output.push_back(fs.getVal());
    //output.push_back(fs.getError());
    //output.push_back(as.getVal());
    //output.push_back(as.getError());
    return output;
}//}}}

void angular(const char outfile[] = "angular")
{//{{{

    TCanvas *c = new TCanvas();
    TH2F *frame = new TH2F("frame","",18,1,19,10,-1,1);
    frame->SetStats(kFALSE);
    frame->SetXTitle("q^{2} [(GeV)^{2}]");
    frame->SetYTitle("F_{L}");
    frame->SetAxisRange(0,1,"Y");
    frame->Draw();

    double x[8]={1.5,3.15,6.49,9.385,11.475,13.52,15.09,17.5};
    double xerr[8]={0.5,1.15,2.09,0.705,1.385,0.66,0.91,1.5};
    double yafb[8],yerrafb[8],yfl[8],yerrfl[8];

    std::vector<double> vbin;
    angular3D_bin(3);
    angular3D_bin(5);
    for(int ibin = 0; ibin < 8; ibin++){
        vbin = angular3D_bin(ibin);
        yfl[ibin]       =vbin.at(0);
        yerrfl[ibin]    =vbin.at(1);
        yafb[ibin]      =vbin.at(2);
        yerrafb[ibin]   =vbin.at(3);
    }
    
    // Check input data
    for(int ibin = 0; ibin < 8; ibin++){
        printf("yafb[%d]=%6.4f +- %6.4f\n",ibin,yafb[ibin],yerrafb[ibin]);
        printf("yfl [%d]=%6.4f +- %6.4f\n",ibin,yfl[ibin],yerrfl[ibin]);
    }

    TGraphAsymmErrors *g_fl  = new TGraphAsymmErrors(8,x,yfl,xerr,xerr,yerrfl,yerrfl);
    g_fl->SetFillColor(2);
    g_fl->SetFillStyle(3001);
    g_fl->Draw("P2");
    c->Print(TString::Format("./plots/%s_fl.pdf",outfile));
    c->Clear();

    frame->SetTitle("");
    frame->SetYTitle("A_{FB}");
    frame->SetXTitle("q^{2} [(GeV)^{2}]");
    frame->SetAxisRange(-1,1,"Y");
    frame->Draw();
    TGraphAsymmErrors *g_afb = new TGraphAsymmErrors(8,x,yafb,xerr,xerr,yerrafb,yerrafb);
    g_afb->SetFillColor(2);
    g_afb->SetFillStyle(3001);
    g_afb->Draw("P2");
    c->Print(TString::Format("./plots/%s_afb.pdf",outfile));
}//}}}

//_________________________________________________________________________________

int main(int argc, char** argv) {
    if (argc <= 2) {
        printf("Need at least 2 arguments.\n");
        printf("./fit [bmass, fl, angular_gen, acceptance, efficiency, angular3D_1a_Sm, angular3D_1b_YpPm, angular3D_2a_PkPl, angular] infile\n");
        printf("Outputs will be stored in ./plots.\n");
        return 0;
    }

    TString func    = argv[1];
    TString infile  = argv[2];
  
    if (func == "bmass") {
        if (argc < 6){
            printf("./fit bmass infile datatype lable cut\n");
        }
        TString datatype = argv[3]; 
        TString label    = argv[4]; 
        TString cut      = argv[5]; 

        TChain* ch = add_chain(datatype, label, cut); 
        if (ch == NULL) gSystem->Exit(0);
        
        const char outfile[]="bmass";
        bmass(outfile); 
    }else if (func == "fl"){
        ch->Add(infile.Data());
        if (ch == NULL) gSystem->Exit(0);
        fl("fl");
    }else if (func == "angular_gen"){
        ch->Add(infile.Data());
        if (ch == NULL) gSystem->Exit(0);
        const char outfile[]="angular_gen";
        angular_gen(outfile);
    }else if (func == "acceptance") {
        ch->Add(infile.Data());
        if (ch == NULL) gSystem->Exit(0);
        for (int iBin = 0; iBin < 8; iBin++) {
            acceptance(iBin);
        }
    }else if (func == "efficiency") {
        ch->Add(infile.Data());
        if (ch == NULL) gSystem->Exit(0);
        for (int iBin = 0; iBin < 8; iBin++) {
            efficiency(iBin);
        }
    }else if (func == "angular"){
        ch->Add(infile.Data());
        if (ch == NULL) gSystem->Exit(0);
        const char outfile[]="angular3D";
        angular(outfile);
        //angular3D_bin(3);
    }else if (func == "angular3D_1a_Sm" || func == "angular3D_1b_YpPm" || func == "angular3D_2a_PkPl" || func == "angular3D_prior"){
        void (*fx)(int, const char*, bool);
        if ( func == "angular3D_1a_Sm" ){
            fx = angular3D_1a_Sm;
        }else if (func == "angular3D_1b_YpPm"){
            fx = angular3D_1b_YpPm;
        }else if (func == "angular3D_2a_PkPl"){
            fx = angular3D_2a_PkPl;
        }else{
            fx = angular3D_prior;
        }
        ch->Add(infile.Data());
        if (ch == NULL) gSystem->Exit(0);
        for (int iBin = 0; iBin < 8; iBin++) {
            fx(iBin,func,false);
        }
    }else{ 
        cerr << "No function available for: " << func.Data() << endl; 
    }
    printf("%lld entries processed.\n",ch->GetEntries());
    gSystem->Exit(0);

    return 0 ;
}
