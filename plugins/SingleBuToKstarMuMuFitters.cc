//vim: sw=4 ts=4 fdm=marker et:

// -----------------------------------------------
//       Author: Xin Shi <Xin.Shi@cern.ch> 
//       Created:   [2013-08-15 Thu 14:54] 
// -----------------------------------------------
#include <stdio.h>
#include <sstream>
#include <sys/stat.h>
#include <math.h>

#include <TSystem.h>
#include <TStyle.h>
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

#include <RooConstVar.h>
#include <RooRealVar.h>
#include <RooAbsPdf.h>
#include <RooAddPdf.h>
#include <RooGaussian.h>
#include <RooBifurGauss.h>
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

// Tags configration
bool is7TeVCheck = false; // Efficiency map for 2011 result...
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
double arrAccParErr2011[8][20]= {{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
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
            0,0,0,0}};

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
            0,0,0,0}};

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
double arrRecPar2012   [8][20] = {{0.302556,-0.063084,-0.032737,-0.022350,0.405606,0.109493,-0.009563,0.099872,0,0,0,0,-0.880160,-0.193247,0.092921,-0.248739,0.440497,0.215391,0.077059,0.059313},
                                  {0.290427,-0.049869,-0.032928,-0.010764,0.075667,-0.015028,0.070311,-0.050658,0.000000,0.000000,0.000000,0.000000,-0.009506,0.088915,-0.054895,0.012487,0.000000,0.000000,0.000000,0.000000},
                                  {0.263419,-0.045855,-0.013515,-0.010111,0.079976,0.030516,-0.021958,-0.040610,0.023867,0.013435,-0.003391,-0.005476,-0.038232,-0.007666,0.030161,0.048279,0.000000,0.000000,0.000000,0.000000},
                                  {0.253280,-0.030755,0.004367,-0.008005,0.086521,0.003089,-0.075045,-0.032558,0.022282,-0.004811,-0.003856,-0.006578,-0.055950,-0.009123,0.080280,0.066337,0.000000,0.000000,0.000000,0.000000},
                                  {0.241065,-0.026626,0.009456,-0.007636,0.094950,-0.003958,-0.050344,-0.017850,0.014244,0.011032,0.002682,-0.006849,-0.047625,0.002464,0.041776,0.034477,0.000000,0.000000,0.000000,0.000000},
                                  {0.242378,-0.016753,0.007127,-0.002982,0.081482,-0.011376,-0.010285,-0.049778,0.023257,-0.004237,0.000603,-0.013435,-0.035768,-0.006341,-0.020483,0.062262,0.000000,0.000000,0.000000,0.000000},
                                  {0.245440,-0.020216,0.003496,-0.005518,0.042054,-0.015924,-0.006110,0.016815,0.011219,0.018159,0.003690,0.007730,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                                  {0.253529,-0.014504,0.003312,-0.005950,0.030220,-0.026910,-0.021029,0.015787,-0.002315,0.002565,0.003259,-0.013746,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000}};
double arrRecParErr2012[8][20] = {{0.018651,0.029866,0.041196,0.048162,0.346595,0.570466,0.779432,0.925547,0,0,0,0,1.235893,2.047440,2.765063,3.285229,1.060858,1.786039,2.403616,2.880058},
                                  {0.010365,0.016132,0.022594,0.026082,0.087765,0.141006,0.194628,0.225204,0.000000,0.000000,0.000000,0.000000,0.133357,0.220169,0.300464,0.347696,0.000000,0.000000,0.000000,0.000000},
                                  {0.006283,0.010361,0.014140,0.016368,0.044351,0.076192,0.102044,0.118048,0.016454,0.028588,0.037906,0.043885,0.056317,0.099729,0.131673,0.152091,0.000000,0.000000,0.000000,0.000000},
                                  {0.009747,0.017043,0.022631,0.026201,0.063636,0.116011,0.151871,0.174953,0.022712,0.040156,0.054132,0.061852,0.076303,0.143072,0.185733,0.213276,0.000000,0.000000,0.000000,0.000000},
                                  {0.006370,0.011494,0.015087,0.017411,0.040269,0.075979,0.097937,0.112485,0.014202,0.025177,0.033383,0.038450,0.047104,0.091007,0.116142,0.133182,0.000000,0.000000,0.000000,0.000000},
                                  {0.008782,0.016211,0.021113,0.024290,0.054285,0.104656,0.133530,0.153182,0.018777,0.033356,0.043661,0.050680,0.062484,0.122775,0.155011,0.177873,0.000000,0.000000,0.000000,0.000000},
                                  {0.005799,0.010933,0.014111,0.016200,0.016482,0.030753,0.039878,0.045879,0.014331,0.026705,0.034612,0.039876,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                                  {0.005179,0.009880,0.012654,0.014518,0.013748,0.026752,0.033939,0.039006,0.011658,0.022922,0.028923,0.033282,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000}};
//}}}

//_________________________________________________________________________________

void bmass(int iBin, const char outfile[] = "bmass")
{//{{{
  bool test = false; 

  RooRealVar Bmass("Bmass", "B^{+/-} mass(GeV/c^{2})", 5.27953-0.28, 5.27953+0.28) ;
  RooRealVar Q2("Q2","q^{2}",0.5,20.);
  RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass),q2range[iBin],0);

  data->Print();

  // Create model and dataset
  // -------------------------------------------------------------------------
  // Gaussian signal 
  RooRealVar mean("mean","mean of gaussians", 5.27, 5.23, 5.32) ;
  RooRealVar sigma1("sigma1","width of Gaussian1", 0.0285, 0.01, 0.05) ;
  RooRealVar sigma2("sigma2","width of Gaussian2", 0.08, 0.05, 0.35) ;
  RooRealVar sigM_frac("sigM_frac","fraction of Gaussians",0,0.,1.);
  RooGaussian sigGauss1("siggauss1","Signal component", Bmass, mean, sigma1) ;  
  RooGaussian sigGauss2("siggauss2","Signal component", Bmass, mean, sigma2) ;  
  RooAddPdf sig("sig","sig",RooArgList(sigGauss1,sigGauss2),RooArgList(sigM_frac));

  // Build Chebychev polynomial p.d.f.  
  RooRealVar a0("a0", "constant", 0.5, -1, 1.) ;
  RooRealVar a1("a1", "linear", 0.6, -1, 1) ;
  RooRealVar a2("a2", "quadratic", 0.1, -1, 1) ;
  RooChebychev bkg("bkg", "Background", Bmass, RooArgSet(a0, a1, a2)) ;

  // Construct signal+background PDF
  RooRealVar nsig("nsig", "number of signal events", 4648, 0, 1E8); 
  RooRealVar nbkg("nbkg", "number of background events", 21472, 0, 1E8);
  RooAddPdf  model("model", "g+c", RooArgList(bkg, sig), RooArgList(nbkg, nsig)) ;
  
  // Print structure of composite p.d.f.
  model.Print("t") ;

  // Fit model to data, save fitresult 
  // ------------------------------------------------------------------------
  RooFitResult* fitres; 
  if (! test) {
    fitres = model.fitTo(*data, Extended(kTRUE), Save(kTRUE)) ;
    fitres->Print("v"); 
  }
  
  // Plot model 
  // ---------------------------------------------------------
  TString title = "B^{+/-} mass";
  int nbins = 50; 
  RooPlot* frame = Bmass.frame(Title(title), Bins(nbins));
  data->plotOn(frame) ;
  model.plotOn(frame) ;

  // Overlay the background component of model with a dashed line
  model.plotOn(frame,Components("bkg"), LineStyle(kDashed)) ;

  // Draw the frame on the canvas
  TCanvas* c = new TCanvas("c", "c", 400, 400); 
  set_root_style(); 
  c->UseCurrentStyle() ;

  gPad->SetLeftMargin(0.15) ;
  frame->GetYaxis()->SetTitleOffset(1.7) ; 
  frame->Draw();

  TPaveText* paveText = new TPaveText(0.17, 0.70, 0.41, 0.88, "NDC"); 
  paveText->SetBorderSize(0);
  paveText->SetFillColor(kWhite);
  paveText->AddText(Form("nsig = %.0f #pm %.0f " , nsig.getVal()  , nsig.getError())); 
  paveText->AddText(Form("nbkg = %.0f #pm %.0f " , nbkg.getVal()  , nbkg.getError())); 
  paveText->AddText(Form("mean = %.3f #pm %.3f " , mean.getVal()  , mean.getError())); 
  paveText->AddText(Form("sigma1 = %.3f #pm %.3f ", sigma1.getVal(), sigma1.getError())); 
  paveText->AddText(Form("sigma2 = %.3f #pm %.3f ", sigma2.getVal(), sigma2.getError())); 
  paveText->AddText(Form("frac = %.3f #pm %.3f ", sigM_frac.getVal(), sigM_frac.getError())); 
  paveText->Draw(); 

  c->Print(TString::Format("./plots/%s_bin%d.pdf",outfile,iBin));
  
  // Persist fit result in root file 
  // -------------------------------------------------------------
  //TFile resf(TString::Format("./plots/%s.root",outfile), "RECREATE") ;
  //gPad->Write("plot"); 
  //if (! test) fitres->Write("fitres") ;
  //resf.Close() ;

  // In a clean ROOT session retrieve the persisted fit result as follows:
  // RooFitResult* r = gDirectory->Get("fitres") ;
  
  delete paveText; 
  delete c;

}//}}}

//_________________________________________________________________________________

std::vector<double> fl_gen_bin(int iBin, const char outfile[] = "fl_gen")
{//{{{
  // From fomula (2) in LHCb 2012 PRL108, 181806(2012)
  // integrated over theta_l and phi: 
  // 
  // 1/Gamma * d^2 Gamma/d cos(theta_K) dq^2 = 3/2 * F_L cos^2(theta_K)
  // + 3/4(1-F_L)(1-cos^2theta_K)
  // 

  bool test = false;

  RooRealVar genCosThetaK("genCosThetaK", "cos#theta_{K}", -1, 1);
  RooRealVar genQ2("genQ2","q^{2}",0.5,20.);
  RooRealVar Q2("Q2","q^{2}",0.5,20.);
  RooRealVar fl("fl", "F_{L}", 0.5, -0.5, 1.5);

  RooGenericPdf f("f", "1.5*fl*genCosThetaK*genCosThetaK+0.75*(1-fl)*(1-genCosThetaK*genCosThetaK)", RooArgSet(genCosThetaK,fl));
  RooDataSet* data;
  
  if (test){
      fl.setVal(0.5);
      data = f.generate(RooArgSet(genCosThetaK,Q2), 10000);
  }else{
      data = new RooDataSet("data","data",ch,RooArgSet(genCosThetaK,genQ2),genQ2range[iBin],0);
  }
  
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
  t1->DrawLatex(.40,.85,TString::Format("%s",genQ2range[iBin]));
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

void fl_gen(const char outfile[] = "fl_gen")
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
        vbin = fl_gen_bin(ibin);
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
    RooRealVar fl("fl", "F_{L}", 0.8, 0., 1.);
    RooRealVar afb("afb", "A_{FB}", 0., -1., 1.);
    RooRealVar fs("fs","F_{S}",0.,0.,1.);//Derive from B0ToKstarJpsi
    RooRealVar as("as","A_{S}",0.,-1,1.);//Derive from B0ToKstarJpsi
    fs.setConstant(kTRUE);
    as.setConstant(kTRUE);

    RooRealVar nsig("nsig","nsig",1E4,1E2,1E8);
    RooRealVar nbkg("nbkg","nbkg",10,0.1,1E3);
    
    RooGenericPdf f_sig("f_sig", "9/16*((2/3*fs+4/3*as*genCosThetaK)*(1-genCosThetaL*genCosThetaL)+(1-fs)*(2*fl*genCosThetaK*genCosThetaK*(1-genCosThetaL*genCosThetaL)+1/2*(1-fl)*(1-genCosThetaK*genCosThetaK)*(1+genCosThetaL*genCosThetaL)+4/3*afb*(1-genCosThetaK*genCosThetaK)*genCosThetaL))", RooArgSet(genCosThetaK,genCosThetaL,fl,afb,fs,as));
    RooGenericPdf f_bkg("f_bkg", "1",RooArgSet());
    //nbkg.setConstant(kTRUE);
    RooAddPdf f("f","f",RooArgList(f_sig,f_bkg),RooArgList(nsig,nbkg));
    RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(genCosThetaK,genCosThetaL,genQ2),genQ2range[iBin],0);

    RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"));

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
    //t1->DrawLatex(.35,.10,TString::Format("nbkg=%5.3f#pm%5.3f",nbkg.getVal(),nbkg.getError()));
    c->Print(TString::Format("./plots/%s_cosk_bin%d.pdf",outfile,iBin));

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
    c->Print(TString::Format("./plots/%s_cosl_bin%d.pdf",outfile,iBin));

    // Make 2-D plot
    TH1 *h1 = data->createHistogram("genCosThetaL,genCosThetaK", 100, 100);
    h1->Draw("LEGO2");
    c->Update();
    c->Print(TString::Format("./plots/%s_2D_bin%d.pdf",outfile,iBin));

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
            //if (measure != 0) 
            //    f+=delta*delta;
            
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
            f += ( f2_square->Integral(xi,xf,yi,yf) - 2*f2_fcn->Integral(xi,xf,yi,yf)*measure + measure*measure*(xf-xi)*(yf-yi) )/(error*error);

            f2_square = 0;
            delete f2_square;
        }
    }
    //printf("FCN in calls = %f\n",f);
    //printf("npar=%d ",npar);
}//}}}

std::vector<double> acceptance(int iBin) // acceptance
{//{{{

    double accUpperBound = 0.09;
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
    int nbinsL = 20;
    int nbinsK = 20;
    float thetaKBins[6]={-1,-0.7,0.,0.4,0.8,1};
    float thetaLBins[7]={-1,-0.7,-0.3,0.,0.3,0.7,1};
    //TH2F h2_ngen("h2_ngen","h2_ngen",nbinsL,thetaLBins,nbinsK,thetaKBins);
    //TH2F h2_nacc("h2_nacc" ,"h2_nacc" ,nbinsL,thetaLBins,nbinsK,thetaKBins); 
    TH2F h2_ngen("h2_ngen","h2_ngen",nbinsL,-1.,1,nbinsK,-1,1);
    TH2F h2_nacc("h2_nacc" ,"h2_nacc" ,nbinsL,-1,1,nbinsK,-1,1); 
    for (int entry = 0; entry < ch->GetEntries(); entry++) {
        ch->GetEntry(entry);
        if (gQ2 > q2rangeup[iBin] || gQ2 < q2rangedn[iBin]) continue;
        h2_ngen.Fill(gCosThetaL,gCosThetaK);
        if ( fabs(gmumeta) < 2.3 && fabs(gmupeta) < 2.3 && gmumpt > 2.5 && gmuppt > 2.5 ) h2_nacc.Fill(gCosThetaL,gCosThetaK);
    }
    
    // Calculate acceptance
    //TH2F h2_acc("h2_acc","",nbinsL,thetaLBins,nbinsK,thetaKBins);
    TH2F h2_acc("h2_acc","",nbinsL,-1,1,nbinsK,-1,1);
    h2_acc.SetAxisRange(0.,1.,"Z");
    for (int i = 1; i <= nbinsL; i++) {
        for (int j = 1; j <= nbinsK; j++) {
            if (h2_ngen.GetBinContent(i,j) == 0 || h2_nacc.GetBinContent(i,j) == 0) {
                printf("WARNING: Acceptance(%d,%d)=%f/%f\n",i,j,h2_nacc.GetBinContent(i,j),h2_ngen.GetBinContent(i,j));
                h2_acc.SetBinContent(i,j,0.);
                h2_acc.SetBinError(i,j,1.);//Set typical value 5%
            }else{
                h2_acc.SetBinContent(i,j,h2_nacc.GetBinContent(i,j)/h2_ngen.GetBinContent(i,j));
                h2_acc.SetBinError(i,j,sqrt(h2_acc.GetBinContent(i,j)*(1.-h2_acc.GetBinContent(i,j))/h2_ngen.GetBinContent(i,j)));
                //printf("INFO: Angular bin(%d,%d)= %f +- %f ( %f / %f).\n",i,j,h2_acc.GetBinContent(i,j),h2_acc.GetBinError(i,j),h2_nacc.GetBinContent(i,j),h2_ngen.GetBinContent(i,j));
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
    gMinuit->DefineParameter( 0, "k0l0",  .01,  1E-3,    -1E+0, 1E+2);
    gMinuit->DefineParameter( 1, "k1l0",   0.,  1E-3,    -1E+1, 1E+2);
    gMinuit->DefineParameter( 2, "k2l0",   0.,  1E-3,    -1E+1, 1E+2);
    gMinuit->DefineParameter( 3, "k3l0",   0.,  1E-3,    -1E+1, 1E+2);
    gMinuit->DefineParameter( 4, "k0l2", 1E-2,  1E-3,    -1E+1, 1E+2);
    gMinuit->DefineParameter( 5, "k1l2",   0.,  1E-3,    -1E+1, 1E+2);
    gMinuit->DefineParameter( 6, "k2l2",   0.,  1E-3,    -1E+1, 1E+2);
    gMinuit->DefineParameter( 7, "k3l2",   0.,  1E-3,    -1E+1, 1E+2);
    gMinuit->DefineParameter( 8, "k0l3",   0.,  1E-3,    -1E+1, 1E+2);
    gMinuit->DefineParameter( 9, "k1l3",   0.,  1E-3,    -1E+1, 1E+2);
    gMinuit->DefineParameter(10, "k2l3",   0.,  1E-3,    -1E+1, 1E+2);
    gMinuit->DefineParameter(11, "k3l3",   0.,  1E-3,    -1E+1, 1E+2);
    gMinuit->DefineParameter(12, "k0l4",-1E-2,  1E-3,    -1E+1, 1E+2);
    gMinuit->DefineParameter(13, "k1l4",   0.,  1E-3,    -1E+1, 1E+2);
    gMinuit->DefineParameter(14, "k2l4",   0.,  1E-3,    -1E+1, 1E+2);
    gMinuit->DefineParameter(15, "k3l4",   0.,  1E-3,    -1E+1, 1E+2);
    gMinuit->DefineParameter(16, "k0l6", 1E-2,  1E-3,    -1E+2, 1E+4);
    gMinuit->DefineParameter(17, "k1l6",   0.,  1E-3,    -1E+2, 1E+4);
    gMinuit->DefineParameter(18, "k2l6",   0.,  1E-3,    -1E+2, 1E+4);
    gMinuit->DefineParameter(19, "k3l6",   0.,  1E-3,    -1E+2, 1E+4);
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
    h2_acc.SetMaximum(accUpperBound);
    h2_acc.Draw("LEGO2");
    latex->DrawLatexNDC(0.35,0.95,TString::Format("Acceptance in Bin%d",iBin));
    
    // Draw FitResult
    f2_model.SetTitle("");
    f2_model.SetMaximum(accUpperBound);
    f2_model.SetLineWidth(1);
    f2_model.Draw("SURF SAME ");
    canvas.Print(TString::Format("./plots/acceptance_2D_bin%d.pdf",iBin));

    //// Draw compare
    double chi2Val=0;
    fcn_binnedChi2_2D(nPar, 0, chi2Val, arrPar, 0);
    printf("Chi2=%f \n",chi2Val);
    //TH2F h2_compFit("h2_compFit","",nbinsL,thetaLBins,nbinsK,thetaKBins);
    TH2F h2_compFit("h2_compFit","",nbinsL,-1,1,nbinsK,-1,1);
    for (int i = 1; i <= nbinsL; i++) {//thetaL
        for (int j = 1; j <= nbinsK; j++) {//thetaK
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
    latex->DrawLatexNDC(0.3,0.95,TString::Format("acceptance_{measured} / acceptance_{fit} in Bin%d",iBin));
    canvas.Update();
    canvas.Print(TString::Format("./plots/acceptance_compFit_2D_bin%d.pdf",iBin));
    
    // Draw pull
    TH1F h_pull("h_pull","",15,-3.,3.);
    h_pull.SetXTitle("Pull");
    for (int i = 1; i <= nbinsL; i++) {//thetaL
        for (int j = 1; j <= nbinsK; j++) {//thetaK
            double _xlo = h2_acc.GetXaxis()->GetBinLowEdge(i);
            double _xhi = h2_acc.GetXaxis()->GetBinUpEdge(i);
            double _ylo = h2_acc.GetYaxis()->GetBinLowEdge(j);
            double _yhi = h2_acc.GetYaxis()->GetBinUpEdge(j);
            if (h2_nacc.GetBinContent(i,j) != 0){
                h_pull.Fill((f2_model.Integral(_xlo,_xhi,_ylo,_yhi)/(_xhi-_xlo)/(_yhi-_ylo)-h2_acc.GetBinContent(i,j))/h2_acc.GetBinError(i,j));
            }
        }
    }
    h_pull.Draw();
    canvas.Update();
    canvas.Print(TString::Format("./plots/acceptance_pull_bin%d.pdf",iBin));

    // Draw projection to cosThetaK
    //TH1D *h_cosk = new TH1D("h_cosk","",nbinsK,thetaKBins);
    TH1D *h_cosk = new TH1D("h_cosk","",nbinsK,-1,1);
    for (int kBin = 1; kBin <= nbinsK; kBin++) {
        float sumGen = 0.;
        float sumAcc = 0.;
        for ( int lBin = 1; lBin <= nbinsL; lBin++) {
            sumGen+=h2_ngen.GetBinContent(lBin,kBin);
            sumAcc+=h2_nacc.GetBinContent(lBin,kBin);
        }
        h_cosk->SetBinContent(kBin,sumAcc/sumGen);
        h_cosk->SetBinError(kBin,sqrt(sumAcc*(sumGen-sumAcc)/pow(sumGen,3)));
    }
    h_cosk->SetStats(0);
    h_cosk->SetMinimum(0.);
    h_cosk->SetMaximum(accUpperBound);
    h_cosk->Draw();
    latex->DrawLatexNDC(0.32,0.95,TString::Format("Projection to cos#theta_{k} in Bin%d",iBin));

    TF1  *f_cosk = new TF1("f_cosk"
                          ,"([0]+[1]*x+[2]*(3*x**2-1)/2+[3]*(5*x**3-3*x)/2)+([4]+[5]*x+[6]*(3*x**2-1)/2+[7]*(5*x**3-3*x)/2)/3+([8]+[9]*x+[10]*(3*x**2-1)/2+[11]*(5*x**3-3*x)/2)*0+([12]+[13]*x+[14]*(3*x**2-1)/2+[15]*(5*x**3-3*x)/2)/5+([16]+[17]*x+[18]*(3*x**2-1)/2+[19]*(5*x**3-3*x)/2)/7"
                          ,-1.,1.);
    for (int iPar = 0; iPar < nPar; iPar++) f_cosk->SetParameter(iPar,arrPar[iPar]);
    f_cosk->Draw("SAME");
    canvas.Update();
    canvas.Print(TString::Format("./plots/acceptance_cosK_bin%d.pdf",iBin));
    
    // Draw projection to cosThetaL
    //TH1D *h_cosl = new TH1D("h_cosl","",nbinsL,thetaLBins);
    TH1D *h_cosl = new TH1D("h_cosl","",nbinsL,-1,1);
    for ( int lBin = 1; lBin <= nbinsL; lBin++) {
        float sumGen = 0.;
        float sumAcc = 0.;
        for (int kBin = 1; kBin <= nbinsK; kBin++) {
            sumAcc+=h2_nacc.GetBinContent(lBin,kBin);
            sumGen+=h2_ngen.GetBinContent(lBin,kBin);
        }
        h_cosl->SetBinContent(lBin,sumAcc/sumGen);
        h_cosl->SetBinError(lBin,sqrt(sumAcc*(sumGen-sumAcc)/pow(sumGen,3)));
    }
    h_cosl->SetStats(0);
    h_cosl->SetMinimum(0.);
    h_cosl->SetMaximum(accUpperBound);
    h_cosl->Draw();
    latex->DrawLatexNDC(0.32,0.95,TString::Format("Projection to cos#theta_{l} in Bin%d",iBin));

    TF1  *f_cosl = new TF1("f_cosl"
                          ,"([0]+[1]*0+[2]*0+[3]*0)+([4]+[5]*0+[6]*0+[7]*0)*x**2+([8]+[9]*0+[10]*0+[11]*0)*x**3+([12]+[13]*0+[14]*0+[15]*0)*x**4+([16]+[17]*0+[18]*0+[19]*0)*x**6"
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
        
        printf("%f,",arrPar[iPar]);
        if (iPar+1 >= nPar) printf("\n");
    }
    for (int i = 0; i < output.size(); i=i+2) {
        printf("%f,",output[i+1]);
        if (i+2 >= output.size()) printf("\n");
    }
    return output;
}//}}}

std::vector<double> efficiency(int iBin) // reconstruction efficiency
{//{{{
    printf("Evaluate reconstruction efficiency for bin#%d\n",iBin);
    double effUpperBound = 0.5;
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
    h2_rec.SetTitleOffset(2,"XY");
    h2_rec.SetXTitle("genCosThetaL");
    h2_rec.SetYTitle("genCosThetaK");
    for (int i = 1; i <= 6; i++) {
        for (int j = 1; j <= 5; j++) {
            // Build from MC samples
            if (h2_nacc.GetBinContent(i,j) == 0 || h2_nreco.GetBinContent(i,j) == 0) {
                printf("WARNING: Efficiency(%d,%d)=0, set error to be 1.\n",i,j);
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
    if (is7TeVCheck){
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
    }else{
        gMinuit->DefineParameter( 0, "k0l0",  .01,  1E-3,    -1E+1, 1E+1);
        gMinuit->DefineParameter( 1, "k1l0",   0.,  1E-3,    -1E+1, 1E+1);
        gMinuit->DefineParameter( 2, "k2l0",   0.,  1E-3,    -1E+1, 1E+1);
        gMinuit->DefineParameter( 3, "k3l0",   0.,  1E-3,    -1E+1, 1E+1);
        gMinuit->DefineParameter( 4, "k0l2", 1E-2,  1E-3,    -1E+1, 1E+1);
        gMinuit->DefineParameter( 5, "k1l2",   0.,  1E-3,    -1E+1, 1E+1);
        gMinuit->DefineParameter( 6, "k2l2",   0.,  1E-3,    -1E+1, 1E+1);
        gMinuit->DefineParameter( 7, "k3l2",   0.,  1E-3,    -1E+1, 1E+1);
        gMinuit->DefineParameter( 8, "k0l3",   0.,  1E-3,    -1E+1, 1E+1);
        gMinuit->DefineParameter( 9, "k1l3",   0.,  1E-3,    -1E+1, 1E+1);
        gMinuit->DefineParameter(10, "k2l3",   0.,  1E-3,    -1E+1, 1E+1);
        gMinuit->DefineParameter(11, "k3l3",   0.,  1E-3,    -1E+1, 1E+1);
        gMinuit->DefineParameter(12, "k0l4",-1E-2,  1E-3,    -1E+1, 1E+1);
        gMinuit->DefineParameter(13, "k1l4",   0.,  1E-3,    -1E+1, 1E+1);
        gMinuit->DefineParameter(14, "k2l4",   0.,  1E-3,    -1E+1, 1E+1);
        gMinuit->DefineParameter(15, "k3l4",   0.,  1E-3,    -1E+1, 1E+1);
        gMinuit->DefineParameter(16, "k0l6", 1E-2,  1E-3,    -1E+2, 1E+3);
        gMinuit->DefineParameter(17, "k1l6",   0.,  1E-3,    -1E+2, 1E+3);
        gMinuit->DefineParameter(18, "k2l6",   0.,  1E-3,    -1E+2, 1E+3);
        gMinuit->DefineParameter(19, "k3l6",   0.,  1E-3,    -1E+2, 1E+3);
    }

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
    h2_rec.SetMaximum(effUpperBound);
    h2_rec.Draw("LEGO2");
    latex->DrawLatexNDC(0.35,0.95,TString::Format("#varepsilon_{RECO} in Bin%d",iBin));
    
    // Draw FitResult
    f2_model.SetTitle("");
    f2_model.SetMaximum(effUpperBound);
    f2_model.SetLineWidth(1);
    f2_model.Draw("SURF SAME ");
    canvas.Print(TString::Format("./plots/recoEff_2D_bin%d.pdf",iBin));

    //// Draw compare
    double chi2Val=0;
    fcn_binnedChi2_2D(nPar, 0, chi2Val, arrPar, 0);
    printf("Chi2(Bin center)=%f \n",chi2Val);
    
    TH2F h2_compFit("h2_compFit","",6,thetaLBins,5,thetaKBins);
    h2_compFit.SetTitleOffset(2,"XY");
    h2_compFit.SetXTitle("genCosThetaL");
    h2_compFit.SetYTitle("genCosThetaK");
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

    // Draw pull
    TH1F h_pull("h_pull","",15,-3.,3.);
    h_pull.SetXTitle("Pull");
    for (int i = 1; i <= 6; i++) {//thetaL
        for (int j = 1; j <= 5; j++) {//thetaK
            double _xlo = h2_rec.GetXaxis()->GetBinLowEdge(i);
            double _xhi = h2_rec.GetXaxis()->GetBinUpEdge(i);
            double _ylo = h2_rec.GetYaxis()->GetBinLowEdge(j);
            double _yhi = h2_rec.GetYaxis()->GetBinUpEdge(j);
            if (h2_rec.GetBinContent(i,j) != 0){
                h_pull.Fill((f2_model.Integral(_xlo,_xhi,_ylo,_yhi)/(_xhi-_xlo)/(_yhi-_ylo)-h2_rec.GetBinContent(i,j))/h2_rec.GetBinError(i,j));
            }
        }
    }
    h_pull.Draw();
    canvas.Update();
    canvas.Print(TString::Format("./plots/recoEff_pull_bin%d.pdf",iBin));
    

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
    h_cosk->SetMaximum(effUpperBound);
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
    //h_cosl->SetMaximum(0.02);
    h_cosl->Draw();
    latex->DrawLatexNDC(0.32,0.95,TString::Format("Projection to cos#theta_{l} in Bin%d",iBin));

    TF1  *f_cosl = new TF1("f_cosl"
                          ,"([0]+[1]*0+[2]*0+[3]*0)+([4]+[5]*0+[6]*0+[7]*0)*x**2+([8]+[9]*0+[10]*0+[11]*0)*x**3+([12]+[13]*0+[14]*0+[15]*0)*x**4+([16]+[17]*0+[18]*0+[19]*0)*x**6"
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
        
        printf("%f,",arrPar[iPar]);
        if (iPar+1 >= nPar) printf("\n");
    }
    for (int i = 0; i < output.size(); i=i+2) {
        printf("%f,",output[i+1]);
        if (i+2 >= output.size()) printf("\n");
    }
    return output;
}//}}}

void drawEfficiency(int iBin)
{//{{{
    
    char *fcn = 0;
    char fcn2011[] = "([0]+[1]*y+[2]*y**2+[3]*y**3)+([4]+[5]*y+[6]*y**2+[7]*y**3)*x**2+([8]+[9]*y+[10]*y**2+[11]*y**3)*x**3+([12]+[13]*y+[14]*y**2+[15]*y**3)*x**4+([1     6]+[17]*y+[18]*y**2+[19]*y**3)*x**6";
    char fcn2012[] = "([0]+[1]*y+[2]*(3*y**2-1)/2+[3]*(5*y**3-3*y)/2)+([4]+[5]*y+[6]*(3*y**2-1)/2+[7]*(5*y**3-3*y)/2)*x**2+([8]+[9]*y+[10]*(3*y**2-1)/2+[11]*(5*y**3-3*y     )/2)*x**3+([12]+[13]*y+[14]*(3*y**2-1)/2+[15]*(5*y**3-3*y)/2)*x**4+([16]+[17]*y+[18]*(3*y**2-1)/2+[19]*(5*y**3-3*y)/2)*x**6";
    if (is7TeVCheck) fcn = fcn2011;
    
    double xx,yy;
    TF2 f2_eff("f2_eff",fcn,-1.,1.,-1.,1.);
    
    // Prepare canvas
    TCanvas *canvas = new TCanvas("canvas");
    f2_eff.SetTitle(";cos#theta_{L};cos#theta_{L}");
    gStyle->SetTitleOffset(2.,"XYZ");
    f2_eff.SetLineWidth(1);
    TLatex *latex = new TLatex();
    
    f2_eff.SetParameters(arrAccPar2011[iBin]);
    f2_eff.Draw("SURF2 Z");
    f2_eff.GetMinimumXY(xx,yy);
    latex->DrawLatexNDC(0.35,0.95,TString::Format("Bin\#%d, #varepsilon_{acceptance} #geq %.2e",iBin, f2_eff.Eval(xx,yy)));
    canvas->Print(TString::Format("./plots/drawAcceptence2011_bin%d.pdf",iBin));

    f2_eff.SetParameters(arrRecPar2011[iBin]);
    f2_eff.SetMaximum(.004);
    f2_eff.SetMinimum(0.);
    f2_eff.Draw("SURF2 Z");
    f2_eff.GetMinimumXY(xx,yy);
    latex->DrawLatexNDC(0.35,0.95,TString::Format("Bin\#%d, #varepsilon_{reconstruction} #geq %.2e",iBin, f2_eff.Eval(xx,yy)));
    canvas->Update();
    canvas->Print(TString::Format("./plots/drawEfficiency2011_bin%d.pdf",iBin));

}//}}}

//_________________________________________________________________________________
std::vector<double> angular2D_bin(int iBin, const char outfile[] = "angular2D")
{//{{{
    // Remark: You must use RooFit!! It's better in unbinned fit.
    //         Extended ML fitis adopted by Mauro, just follow!!
    
    is7TeVCheck = true;// 8TeV efficiency NOT OK.

    // Read data
    RooRealVar CosThetaK("CosThetaK", "cos#theta_{K}", -1., 1.);
    RooRealVar CosThetaL("CosThetaL", "cos#theta_{L}", -1., 1.);
    RooRealVar Q2("Q2","q^{2}",0.5,20.);
    RooRealVar fl("fl", "F_{L}", 0.5, 0., 1.0);
    RooRealVar afb("afb", "A_{FB}", 0., -1., 1.);
    RooRealVar fs("fs","F_{S}",0.0129254,0.,1.);
    RooRealVar as("as","A_{S}",-0.0975919,-1.,1.);

    if (iBin != 3 && iBin != 5){
        // 2011 cross check
        fs.setVal(0.0129254);
        fs.setAsymError(-0.00898344,0.0101371);
        as.setVal(-0.0975919);
        as.setAsymError(-0.00490805,0.0049092);

        // read parameter from datacard
        //fs.setVal(readParam(3,"fs ",0));
        //fs.setAsymError(readParam(3,"fs ",1),readParam(3,"fs ",2));
        //as.setVal(readParam(3,"as ",0));
        //fs.setAsymError(readParam(3,"as ",1),readParam(3,"as ",2));
    }
        // Efficiency and acceptance
    double *arrAccPar, *arrAccParErr, *arrRecPar, *arrRecParErr;
    if (is7TeVCheck){
        arrAccPar       = arrAccPar2011[iBin];
        arrAccParErr    = arrAccParErr2011[iBin];
        arrRecPar       = arrRecPar2011[iBin];
        arrRecParErr    = arrRecParErr2011[iBin];
    }else{
        arrAccPar       = arrAccPar2012[iBin];
        arrAccParErr    = arrAccParErr2012[iBin];
        arrRecPar       = arrRecPar2012[iBin];
        arrRecParErr    = arrRecParErr2012[iBin];
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
    RooArgSet f_sigA_argset(CosThetaL,CosThetaK);
    f_sigA_argset.add(RooArgSet(fl,afb,fs,as));
    f_sigA_argset.add(RooArgSet(accK0L0,accK1L0,accK2L0,accK3L0));
    f_sigA_argset.add(RooArgSet(accK0L2,accK1L2,accK2L2,accK3L2));
    f_sigA_argset.add(RooArgSet(recK0L0,recK1L0,recK2L0,recK3L0));
    f_sigA_argset.add(RooArgSet(recK0L2,recK1L2,recK2L2,recK3L2));
    TString f_sigA_format;
    TString f_ang_format = "9/16*((2/3*fs+4/3*as*CosThetaK)*(1-CosThetaL**2)+(1-fs)*(2*fl*CosThetaK**2*(1-CosThetaL**2)+1/2*(1-fl)*(1-CosThetaK**2)*(1+CosThetaL**2)+4/3*afb*(1-CosThetaK**2)*CosThetaL))";
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
        f_acc_L2 = "(accK0L2+accK1L2*CosThetaK+accK2L2*(3*CosThetaK**2-1)/2+accK3L2*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**2";
        f_acc_L3 = "(accK0L3+accK1L3*CosThetaK+accK2L3*(3*CosThetaK**2-1)/2+accK3L3*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**3";
        f_acc_L4 = "(accK0L4+accK1L4*CosThetaK+accK2L4*(3*CosThetaK**2-1)/2+accK3L4*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**4";
        f_acc_L6 = "(accK0L6+accK1L6*CosThetaK+accK2L6*(3*CosThetaK**2-1)/2+accK3L6*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**6";
        f_rec_L0 = "(recK0L0+recK1L0*CosThetaK+recK2L0*(3*CosThetaK**2-1)/2+recK3L0*(5*CosThetaK**3-3*CosThetaK)/2)";
        f_rec_L2 = "(recK0L2+recK1L2*CosThetaK+recK2L2*(3*CosThetaK**2-1)/2+recK3L2*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**2";
        f_rec_L3 = "(recK0L3+recK1L3*CosThetaK+recK2L3*(3*CosThetaK**2-1)/2+recK3L3*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**3";
        f_rec_L4 = "(recK0L4+recK1L4*CosThetaK+recK2L4*(3*CosThetaK**2-1)/2+recK3L4*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**4";
        f_rec_L6 = "(recK0L6+recK1L6*CosThetaK+recK2L6*(3*CosThetaK**2-1)/2+recK3L6*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**6";
    }

    if (iBin == 0) {
        f_sigA_argset.add(RooArgSet(accK0L4,accK1L4,accK2L4,accK3L4));
        f_sigA_argset.add(RooArgSet(accK0L6,accK1L6,accK2L6,accK3L6));
        f_sigA_argset.add(RooArgSet(recK0L4,recK1L4,recK2L4,recK3L4));
        f_sigA_argset.add(RooArgSet(recK0L6,recK1L6,recK2L6,recK3L6));
        f_sigA_format = TString::Format("(%s+%s+%s+%s)*(%s+%s+%s+%s)*%s",f_acc_L0.Data(),f_acc_L2.Data(),f_acc_L4.Data(),f_acc_L6.Data(),f_rec_L0.Data(),f_rec_L2.Data(),f_rec_L4.Data(),f_rec_L6.Data(),f_ang_format.Data());
    }else if (iBin == 1) {
        f_sigA_argset.add(RooArgSet(accK0L4,accK1L4,accK2L4,accK3L4));
        f_sigA_argset.add(RooArgSet(recK0L4,recK1L4,recK2L4,recK3L4));
        f_sigA_format = TString::Format("(%s+%s+%s)*(%s+%s+%s)*%s",f_acc_L0.Data(),f_acc_L2.Data(),f_acc_L4.Data(),f_rec_L0.Data(),f_rec_L2.Data(),f_rec_L4.Data(),f_ang_format.Data());
    }else if (iBin > 1 && iBin < 6) {
        f_sigA_argset.add(RooArgSet(accK0L3,accK1L3,accK2L3,accK3L3));
        f_sigA_argset.add(RooArgSet(accK0L4,accK1L4,accK2L4,accK3L4));
        f_sigA_argset.add(RooArgSet(recK0L3,recK1L3,recK2L3,recK3L3));
        f_sigA_argset.add(RooArgSet(recK0L4,recK1L4,recK2L4,recK3L4));
        f_sigA_format = TString::Format("(%s+%s+%s+%s)*(%s+%s+%s+%s)*%s",f_acc_L0.Data(),f_acc_L2.Data(),f_acc_L3.Data(),f_acc_L4.Data(),f_rec_L0.Data(),f_rec_L2.Data(),f_rec_L3.Data(),f_rec_L4.Data(),f_ang_format.Data());
    }else{
        f_sigA_argset.add(RooArgSet(accK0L3,accK1L3,accK2L3,accK3L3));
        f_sigA_argset.add(RooArgSet(recK0L3,recK1L3,recK2L3,recK3L3));
        f_sigA_format = TString::Format("(%s+%s+%s)*(%s+%s+%s)*%s",f_acc_L0.Data(),f_acc_L2.Data(),f_acc_L3.Data(),f_rec_L0.Data(),f_rec_L2.Data(),f_rec_L3.Data(),f_ang_format.Data());
    }
        // angular map of signal
    RooGenericPdf f_sigA("f_sigA", f_sigA_format,f_sigA_argset);
    RooRealVar nsig("nsig","nsig",1E4,1E2,1E8);
    RooExtendPdf f_ext("f_ext","f_ext",f_sigA,nsig);
    
    // Get data and apply unbinned fit
    RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(CosThetaK,CosThetaL,Q2),q2range[iBin],0);
    RooFitResult *f_fitresult = f_ext.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"));

    // Draw the frame on the canvas
    TCanvas* c = new TCanvas("c");
    RooPlot* framecosk = CosThetaK.frame(); 
    data->plotOn(framecosk,Binning(20)); 
    f_ext.plotOn(framecosk); 

    framecosk->SetTitle("");
    framecosk->SetMinimum(0);
    framecosk->Draw();

    TLatex *t1 = new TLatex();
    t1->SetNDC();
    double fixNDC = 0;
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    c->Print(TString::Format("./plots/%s_cosk_bin%d.pdf",outfile,iBin));

    // Draw projection to CosThetaL
    RooPlot* framecosl = CosThetaL.frame(); 
    data->plotOn(framecosl,Binning(20)); 
    f_ext.plotOn(framecosl); 

    framecosl->SetTitle("");
    framecosl->SetMinimum(0);
    framecosl->Draw();

    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    c->Update();
    c->Print(TString::Format("./plots/%s_cosl_bin%d.pdf",outfile,iBin));

    // Make 2-D plot
    TH1 *h1 = data->createHistogram("CosThetaL,CosThetaK", 6, 5);
    h1->SetXTitle("CosThetaL");
    h1->SetYTitle("CosThetaK");
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
void writeParam(int iBin, const char parName[], double *val, int nVal=2, bool overwrite=true)
{//{{{
    struct stat fiBuff;
    FILE *fi = 0;
    if (stat(TString::Format("fitParameters%d.txt",iBin),&fiBuff) == 0){
        rename(TString::Format("fitParameters%d.txt",iBin),TString::Format("fitParameters%d.txt.temp",iBin));
        fi = fopen(TString::Format("fitParameters%d.txt.temp",iBin),"r");
    }else{
        fi = fopen(TString::Format("fitParameters%d.txt.temp",iBin),"w");
    }
    
    bool parExist = false;
    char lineBuff[100];
    memset(lineBuff,' ',100*sizeof(char));
    FILE *fp = fopen(TString::Format("fitParameters%d.txt",iBin),"w");
    while(fgets(lineBuff,100,fi) != NULL ){
        if ( strstr(lineBuff,parName) != NULL ){
            fprintf(fp,"%s",parName);
            int iVal = 0;
            while(iVal < nVal){
                fprintf(fp," %f",val[iVal]);
                iVal++;
            }
            fprintf(fp,"\n");
            parExist = true;
        }else{
            fprintf(fp,"%s",lineBuff);
        }
    }
    if (parExist == false){
        fprintf(fp,"%s",parName);
        int iVal = 0;
        while(iVal < nVal){
            fprintf(fp," %f",val[iVal]);
            iVal++;
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
    fclose(fi);
    remove(TString::Format("fitParameters%d.txt.temp",iBin));
}//}}}

void angular3D_1a_Sm(int iBin, const char outfile[] = "angular3D_1a_Sm", bool keepParam = false)
{//{{{
    // Fit to signal simulation by YsSm+YcCm to determine Sm
    RooRealVar Bmass("Bmass","M_{K^{*}#Mu#Mu}",5.,5.56);
    RooRealVar Q2("Q2","q^{2}",0.5,20.);
    
    // Create parameters and PDFs
        // Signal double gaussian
    RooRealVar sigGauss_mean("sigGauss_mean","M_{K*#Mu#Mu}",5.28,5.25,5.30);
    RooRealVar sigGauss1_sigma("sigGauss1_sigma","#sigma_{1}",.03,.01,.05);
    RooRealVar sigGauss2_sigma("sigGauss2_sigma","#sigma_{2}",.08,.05,.35);
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
    
    RooRealVar nsig("nsig","nsig",0,1E8);
    RooRealVar nbkg("nbkg","nbkg",0,1E8);
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
        double val[3]={0,0,0};
        writeParam(iBin, "iBin", new double((double)iBin), 1);
        if (is7TeVCheck){
            writeParam(iBin, "mode", new double(2011), 1);
        }else{
            writeParam(iBin, "mode", new double(2012), 1);
        }
        val[0]=sigGauss1_sigma.getVal();val[1]=sigGauss1_sigma.getError();
        writeParam(iBin, "sigGauss1_sigma", val);
        val[0]=sigGauss2_sigma.getVal();val[1]=sigGauss2_sigma.getError();
        writeParam(iBin, "sigGauss2_sigma", val);
        val[0]=sigM_frac.getVal();val[1]=sigM_frac.getError();
        writeParam(iBin, "sigM_frac", val);
    }
}//}}}
void angular3D_1b_YpPm(int iBin, const char outfile[] = "angular3D_1b_YpPm", bool keepParam = false)
{//{{{
    if (iBin ==0 || iBin%2 == 1){
        if (keepParam){
            double val[3]={1,0,0};
            writeParam(iBin, "bkgGauss1_mean1", val);
            writeParam(iBin, "bkgGauss1_mean2", val);
            writeParam(iBin, "bkgGauss1_sigma1", val);
            writeParam(iBin, "bkgGauss1_sigma2", val);
            writeParam(iBin, "bkgM_frac1", val);
            writeParam(iBin, "bkgGauss2_mean1", val);
            writeParam(iBin, "bkgGauss2_mean2", val);
            writeParam(iBin, "bkgGauss2_sigma1", val);
            writeParam(iBin, "bkgGauss2_sigma2", val);
            writeParam(iBin, "bkgM_frac2", val);
            writeParam(iBin, "bkgM_frac12", val);
            val[0]=0;
            writeParam(iBin, "nbkgPeak", val);
        }
        return;
    }

    // Fit to control channel simulations by YpPm to determine Yp,Pm.
    RooRealVar Bmass("Bmass","M_{K^{*}#Mu#Mu}",5.,5.56);
    RooRealVar Q2("Q2","q^{2}",0.5,20.);
    
    // Create peak background distribution
    RooRealVar bkgGauss1_mean1("bkgGauss1_mean1","M_{K*#Mu#Mu}",5.05,5.,5.12);
    RooRealVar bkgGauss1_mean2("bkgGauss1_mean2","M_{K*#Mu#Mu}",5.0,4.8,5.20);
    RooRealVar bkgGauss2_mean1("bkgGauss2_mean1","M_{K*#Mu#Mu}",5.40,5.35,5.45);
    RooRealVar bkgGauss2_mean2("bkgGauss2_mean2","M_{K*#Mu#Mu}",5.40,5.35,5.45);
    RooRealVar bkgGauss1_sigma1("bkgGauss1_sigma1","#sigma_{11}",.03,.01,.08);
    RooRealVar bkgGauss1_sigma2("bkgGauss1_sigma2","#sigma_{12}",.12,.08,.50);
    RooRealVar bkgGauss2_sigma1("bkgGauss2_sigma1","#sigma_{21}",.03,.01,.05);
    RooRealVar bkgGauss2_sigma2("bkgGauss2_sigma2","#sigma_{22}",.12,.05,.50);
    RooRealVar bkgM_frac1("bkgM_frac1","bkgM_frac1",1.,0.,1.);
    RooRealVar bkgM_frac2("bkgM_frac2","bkgM_frac2",1.,0.,1.);
    RooRealVar bkgM_frac12("bkgM_frac12","bkgM_frac12",0.,0.,1.);
    RooGaussian f_bkgPeakMGauss11("f_bkgPeakMGauss11","f_bkgPeakMGauss11", Bmass, bkgGauss1_mean1, bkgGauss1_sigma1);
    RooGaussian f_bkgPeakMGauss12("f_bkgPeakMGauss12","f_bkgPeakMGauss12", Bmass, bkgGauss1_mean2, bkgGauss1_sigma2);
    RooGaussian f_bkgPeakMGauss21("f_bkgPeakMGauss21","f_bkgPeakMGauss21", Bmass, bkgGauss2_mean1, bkgGauss2_sigma1);
    RooGaussian f_bkgPeakMGauss22("f_bkgPeakMGauss22","f_bkgPeakMGauss22", Bmass, bkgGauss2_mean2, bkgGauss2_sigma2);
    RooAddPdf f_bkgPeakM1("f_bkgPeakM1","f_bkgPeakM1", RooArgList(f_bkgPeakMGauss11, f_bkgPeakMGauss12), bkgM_frac1);
    RooAddPdf f_bkgPeakM2("f_bkgPeakM2","f_bkgPeakM2", RooArgList(f_bkgPeakMGauss21, f_bkgPeakMGauss22), bkgM_frac2);
    RooAddPdf f_bkgPeakM12("f_bkgPeakM12","f_bkgPeakM12", RooArgList(f_bkgPeakM1,f_bkgPeakM2), bkgM_frac12);
    
    RooRealVar nbkgPeak("nbkgPeak","nbkgPeak",1E1,1,1E7);
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
        double val[3]={0,0,0};
        val[0] = bkgGauss1_mean1.getVal();val[1] = bkgGauss1_mean1.getError();
        writeParam(iBin, "bkgGauss1_mean1", val);
        val[0] = bkgGauss1_mean2.getVal();val[1] = bkgGauss1_mean2.getError();
        writeParam(iBin, "bkgGauss1_mean2", val);
        val[0] = bkgGauss1_sigma1.getVal();val[1] = bkgGauss1_sigma1.getError();
        writeParam(iBin, "bkgGauss1_sigma1", val);
        val[0] = bkgGauss1_sigma2.getVal();val[1] = bkgGauss1_sigma2.getError();
        writeParam(iBin, "bkgGauss1_sigma2", val);
        val[0] = bkgM_frac1.getVal();val[1] = bkgM_frac1.getError();
        writeParam(iBin, "bkgM_frac1", val);
        val[0] = bkgGauss2_mean1.getVal();val[1] = bkgGauss2_mean1.getError();
        writeParam(iBin, "bkgGauss2_mean1", val);
        val[0] = bkgGauss2_mean2.getVal();val[1] = bkgGauss2_mean2.getError();
        writeParam(iBin, "bkgGauss2_mean2", val);
        val[0] = bkgGauss2_sigma1.getVal();val[1] = bkgGauss2_sigma1.getError();
        writeParam(iBin, "bkgGauss2_sigma1", val);
        val[0] = bkgGauss2_sigma2.getVal();val[1] = bkgGauss2_sigma2.getError();
        writeParam(iBin, "bkgGauss2_sigma2", val);
        val[0] = bkgM_frac2.getVal();val[1] = bkgM_frac2.getError();
        writeParam(iBin, "bkgM_frac2", val);
        val[0] = bkgM_frac12.getVal();val[1] = bkgM_frac12.getError();
        writeParam(iBin, "bkgM_frac12", val);
        if (is7TeVCheck){
            switch (iBin) {
                case 2:
                    val[0]=470;val[1]=11;
                    break;
                case 4:
                    val[0]=155;val[1]=6.9;
                    break;
                case 6:
                    val[0]=6.8;val[1]=1.4;
                    break;
                default:
                    val[0] = nbkgPeak.getVal();val[1] = nbkgPeak.getError();
            }
        }else{
            switch (iBin) {
                case 2:
                    val[0]=1410;val[1]=33;
                    break;
                case 4:
                    val[0]=465;val[1]=20.7;
                    break;
                case 6:
                    val[0]=20.4;val[1]=4.2;
                    break;
                default:
                    val[0] = nbkgPeak.getVal();val[1] = nbkgPeak.getError();
            }
            val[0] = nbkgPeak.getVal();val[1] = nbkgPeak.getError();
        }
        writeParam(iBin, "nbkgPeak", val);
    }
}//}}}
void angular3D_2a_PkPl(int iBin, const char outfile[] = "angular3D_2a_PkPl", bool keepParam = false)
{//{{{
    // Gaussian constraint on yields and mass is needed.
    if (iBin ==0 || iBin%2 == 1){
        // Pm is flat(and the yield is 0) for bins other than 2,4,6
        if (keepParam){
            double val[3]={0,0,0};
            writeParam(iBin, "bkgPeakL_c1", val);
            writeParam(iBin, "bkgPeakL_c2", val);
            writeParam(iBin, "bkgPeakL_c3", val);
            writeParam(iBin, "bkgPeakL_c4", val);
            writeParam(iBin, "bkgPeakK_c1", val);
            writeParam(iBin, "bkgPeakK_c2", val);
            writeParam(iBin, "bkgPeakK_c3", val);
            writeParam(iBin, "bkgPeakK_c4", val);
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
    RooRealVar nbkgPeak("nbkgPeak","nbkgPeak",50,0.,1E4);
    RooExtendPdf f_bkgPeakA_ext("f_bkgPeakA_ext","f_bkgPeakA_ext",f_bkgPeakA,nbkgPeak);

    // Gaussian Constraint
    RooGaussian gaus_nbkgPeak("gaus_nbkgPeak","gaus_nbkgPeak",nbkgPeak,RooConst(readParam(iBin,"nbkgPeak ", 0)),RooConst(readParam(iBin, "nbkgPeak ", 1)));
    
    // Get data
    RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(CosThetaK,CosThetaL,Q2),q2range[iBin],0);
    RooFitResult *f_fitresult = f_bkgPeakA_ext.fitTo(*data,Save(kTRUE),Extended(),Minimizer("Minuit"),ExternalConstraints(gaus_nbkgPeak));

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
        double val[3] = {0,0,0};
        val[0] = bkgPeakL_c1.getVal();val[1] = bkgPeakL_c1.getError();
        writeParam(iBin, "bkgPeakL_c1", val);
        val[0] = bkgPeakL_c2.getVal();val[1] = bkgPeakL_c2.getError();
        writeParam(iBin, "bkgPeakL_c2", val);
        val[0] = bkgPeakL_c3.getVal();val[1] = bkgPeakL_c3.getError();
        writeParam(iBin, "bkgPeakL_c3", val);
        val[0] = bkgPeakL_c4.getVal();val[1] = bkgPeakL_c4.getError();
        writeParam(iBin, "bkgPeakL_c4", val);
        val[0] = bkgPeakK_c1.getVal();val[1] = bkgPeakK_c1.getError();
        writeParam(iBin, "bkgPeakK_c1", val);
        val[0] = bkgPeakK_c2.getVal();val[1] = bkgPeakK_c2.getError();
        writeParam(iBin, "bkgPeakK_c2", val);
        val[0] = bkgPeakK_c3.getVal();val[1] = bkgPeakK_c3.getError();
        writeParam(iBin, "bkgPeakK_c3", val);
        val[0] = bkgPeakK_c4.getVal();val[1] = bkgPeakK_c4.getError();
        writeParam(iBin, "bkgPeakK_c4", val);
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
    RooRealVar bkgCombL_c1("bkgCombL_c1","c1",0.,-2.5,2.5);
    RooRealVar bkgCombL_c2("bkgCombL_c2","c2",0.,-2.5,2.5);
    RooRealVar bkgCombL_c3("bkgCombL_c3","c3",0.,-2.5,2.5);
    RooRealVar bkgCombL_c4("bkgCombL_c4","c4",0.,-2.5,2.5);
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
    RooRealVar bkgCombK_c1("bkgCombK_c1","c1",0.,-2.5,2.5);
    RooRealVar bkgCombK_c2("bkgCombK_c2","c2",0.,-2.5,2.5);
    RooRealVar bkgCombK_c3("bkgCombK_c3","c3",0.,-5,5);
    RooRealVar bkgCombK_c4("bkgCombK_c4","c4",0.,-5,5);
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
        double val[3] = {0,0,0};
        val[0] = bkgCombL_c1.getVal();val[1] = bkgCombL_c1.getError();
        writeParam(iBin, "bkgCombL_c1", val);
        val[0] = bkgCombL_c2.getVal();val[1] = bkgCombL_c2.getError();
        writeParam(iBin, "bkgCombL_c2", val);
        val[0] = bkgCombL_c3.getVal();val[1] = bkgCombL_c3.getError();
        writeParam(iBin, "bkgCombL_c3", val);
        val[0] = bkgCombL_c4.getVal();val[1] = bkgCombL_c4.getError();
        writeParam(iBin, "bkgCombL_c4", val);
        val[0] = bkgCombK_c1.getVal();val[1] = bkgCombK_c1.getError();
        writeParam(iBin, "bkgCombK_c1", val);
        val[0] = bkgCombK_c2.getVal();val[1] = bkgCombK_c2.getError();
        writeParam(iBin, "bkgCombK_c2", val);
        val[0] = bkgCombK_c3.getVal();val[1] = bkgCombK_c3.getError();
        writeParam(iBin, "bkgCombK_c3", val);
        val[0] = bkgCombK_c4.getVal();val[1] = bkgCombK_c4.getError();
        writeParam(iBin, "bkgCombK_c4", val);
    }
}//}}}

std::vector<double> angular3D_bin(int iBin, const char outfile[] = "angular3D")
{//{{{
    // Remark: You must use RooFit!! It's better in unbinned fit.
    //         Extended ML fit is adopted by Mauro, just follow!!
    
    // Read data
    RooRealVar CosThetaK("CosThetaK", "cos#theta_{K}", -1., 1.);
    RooRealVar CosThetaL("CosThetaL", "cos#theta_{L}", -1., 1.);
    RooRealVar Bmass("Bmass","M_{K^{*}#Mu#Mu}",5.,5.56);
    RooRealVar Q2("Q2","q^{2}",0.5,20.);

    // Create parameters and PDFs
        // Signal double gaussian
    RooRealVar sigGauss_mean("sigGauss_mean","M_{K*#Mu#Mu}",5.28,5.26,5.30);
    RooRealVar sigGauss1_sigma("sigGauss1_sigma","#sigma_{1}",readParam(iBin,"sigGauss1_sigma ",0));
    sigGauss1_sigma.setError(readParam(iBin,"sigGauss1_sigma ",1));
    RooRealVar sigGauss2_sigma("sigGauss2_sigma","#sigma_{2}",readParam(iBin,"sigGauss2_sigma ",0));
    sigGauss2_sigma.setError(readParam(iBin,"sigGauss2_sigma ",1));
    RooRealVar sigM_frac("sigM_frac","sigM_frac",readParam(iBin,"sigM_frac ",0));
    sigM_frac.setError(readParam(iBin,"sigM_frac ",1));
        // Angular parameters
    RooRealVar afb("afb", "A_{FB}", 0., -1., 1.);
    RooRealVar fl("fl", "F_{L}", 0.8, 0., 1.);
    RooRealVar fs("fs","F_{S}",0.01,0.,1.);//Derive from B0ToKstarJpsi, Bin3
    RooRealVar as("as","A_{S}",-0.1,-1.,1.);//Derive from B0ToKstarJpsi, Bin3
    if (iBin != 3 && iBin != 5){
        // 2011 cross check
        fs.setVal(0.0129254);
        fs.setAsymError(-0.00898344,0.0101371);
        as.setVal(-0.0975919);
        as.setAsymError(-0.00490805,0.0049092);

        // read parameter from datacard
        //fs.setVal(readParam(3,"fs ",0));
        //fs.setAsymError(readParam(3,"fs ",1),readParam(3,"fs ",2));
        //as.setVal(readParam(3,"as ",0));
        //fs.setAsymError(readParam(3,"as ",1),readParam(3,"as ",2));
    }
        // Efficiency and acceptance
    double *arrAccPar, *arrAccParErr, *arrRecPar, *arrRecParErr;
    if (is7TeVCheck){
        arrAccPar       = arrAccPar2011[iBin];
        arrAccParErr    = arrAccParErr2011[iBin];
        arrRecPar       = arrRecPar2011[iBin];
        arrRecParErr    = arrRecParErr2011[iBin];
    }else{
        arrAccPar       = arrAccPar2012[iBin];
        arrAccParErr    = arrAccParErr2012[iBin];
        arrRecPar       = arrRecPar2012[iBin];
        arrRecParErr    = arrRecParErr2012[iBin];
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
    RooArgSet f_sigA_argset(CosThetaL,CosThetaK);
    f_sigA_argset.add(RooArgSet(fl,afb,fs,as));
    f_sigA_argset.add(RooArgSet(accK0L0,accK1L0,accK2L0,accK3L0));
    f_sigA_argset.add(RooArgSet(accK0L2,accK1L2,accK2L2,accK3L2));
    f_sigA_argset.add(RooArgSet(recK0L0,recK1L0,recK2L0,recK3L0));
    f_sigA_argset.add(RooArgSet(recK0L2,recK1L2,recK2L2,recK3L2));
    TString f_sigA_format;
    TString f_ang_format = "9/16*((2/3*fs+4/3*as*CosThetaK)*(1-CosThetaL**2)+(1-fs)*(2*fl*CosThetaK**2*(1-CosThetaL**2)+1/2*(1-fl)*(1-CosThetaK**2)*(1+CosThetaL**2)+4/3*afb*(1-CosThetaK**2)*CosThetaL))";
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
        f_acc_L2 = "(accK0L2+accK1L2*CosThetaK+accK2L2*(3*CosThetaK**2-1)/2+accK3L2*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**2";
        f_acc_L3 = "(accK0L3+accK1L3*CosThetaK+accK2L3*(3*CosThetaK**2-1)/2+accK3L3*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**3";
        f_acc_L4 = "(accK0L4+accK1L4*CosThetaK+accK2L4*(3*CosThetaK**2-1)/2+accK3L4*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**4";
        f_acc_L6 = "(accK0L6+accK1L6*CosThetaK+accK2L6*(3*CosThetaK**2-1)/2+accK3L6*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**6";
        f_rec_L0 = "(recK0L0+recK1L0*CosThetaK+recK2L0*(3*CosThetaK**2-1)/2+recK3L0*(5*CosThetaK**3-3*CosThetaK)/2)";
        f_rec_L2 = "(recK0L2+recK1L2*CosThetaK+recK2L2*(3*CosThetaK**2-1)/2+recK3L2*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**2";
        f_rec_L3 = "(recK0L3+recK1L3*CosThetaK+recK2L3*(3*CosThetaK**2-1)/2+recK3L3*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**3";
        f_rec_L4 = "(recK0L4+recK1L4*CosThetaK+recK2L4*(3*CosThetaK**2-1)/2+recK3L4*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**4";
        f_rec_L6 = "(recK0L6+recK1L6*CosThetaK+recK2L6*(3*CosThetaK**2-1)/2+recK3L6*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**6";
    }

    if (iBin == 0) {
        f_sigA_argset.add(RooArgSet(accK0L4,accK1L4,accK2L4,accK3L4));
        f_sigA_argset.add(RooArgSet(accK0L6,accK1L6,accK2L6,accK3L6));
        f_sigA_argset.add(RooArgSet(recK0L4,recK1L4,recK2L4,recK3L4));
        f_sigA_argset.add(RooArgSet(recK0L6,recK1L6,recK2L6,recK3L6));
        f_sigA_format = TString::Format("(%s+%s+%s+%s)*(%s+%s+%s+%s)*%s",f_acc_L0.Data(),f_acc_L2.Data(),f_acc_L4.Data(),f_acc_L6.Data(),f_rec_L0.Data(),f_rec_L2.Data(),f_rec_L4.Data(),f_rec_L6.Data(),f_ang_format.Data());
    }else if (iBin == 1) {
        f_sigA_argset.add(RooArgSet(accK0L4,accK1L4,accK2L4,accK3L4));
        f_sigA_argset.add(RooArgSet(recK0L4,recK1L4,recK2L4,recK3L4));
        f_sigA_format = TString::Format("(%s+%s+%s)*(%s+%s+%s)*%s",f_acc_L0.Data(),f_acc_L2.Data(),f_acc_L4.Data(),f_rec_L0.Data(),f_rec_L2.Data(),f_rec_L4.Data(),f_ang_format.Data());
    }else if (iBin > 1 && iBin < 6) {
        f_sigA_argset.add(RooArgSet(accK0L3,accK1L3,accK2L3,accK3L3));
        f_sigA_argset.add(RooArgSet(accK0L4,accK1L4,accK2L4,accK3L4));
        f_sigA_argset.add(RooArgSet(recK0L3,recK1L3,recK2L3,recK3L3));
        f_sigA_argset.add(RooArgSet(recK0L4,recK1L4,recK2L4,recK3L4));
        f_sigA_format = TString::Format("(%s+%s+%s+%s)*(%s+%s+%s+%s)*%s",f_acc_L0.Data(),f_acc_L2.Data(),f_acc_L3.Data(),f_acc_L4.Data(),f_rec_L0.Data(),f_rec_L2.Data(),f_rec_L3.Data(),f_rec_L4.Data(),f_ang_format.Data());
    }else{
        f_sigA_argset.add(RooArgSet(accK0L3,accK1L3,accK2L3,accK3L3));
        f_sigA_argset.add(RooArgSet(recK0L3,recK1L3,recK2L3,recK3L3));
        f_sigA_format = TString::Format("(%s+%s+%s)*(%s+%s+%s)*%s",f_acc_L0.Data(),f_acc_L2.Data(),f_acc_L3.Data(),f_rec_L0.Data(),f_rec_L2.Data(),f_rec_L3.Data(),f_ang_format.Data());
    }
        // angular map of signal
    RooGenericPdf f_sigA("f_sigA", f_sigA_format,f_sigA_argset);

    // Create signal distribution
        // mass distro of signal
    RooGaussian f_sigMGauss1("f_sigMGauss1","f_sigMGauss1", Bmass, sigGauss_mean, sigGauss1_sigma);//double gaussian with shared mean
    RooGaussian f_sigMGauss2("f_sigMGauss2","f_sigMGauss2", Bmass, sigGauss_mean, sigGauss2_sigma);//double gaussian with shared mean
    RooAddPdf f_sigM("f_sigM","f_sigM", RooArgList(f_sigMGauss1, f_sigMGauss2), sigM_frac);
    RooProdPdf f_sig("f_sig","f_sig",f_sigM,f_sigA);
    printf("INFO: f_sig prepared.\n");

    // Create combinatorial background distribution
    RooRealVar bkgCombM_c("bkgCombM_c","c1",0.,-20,1);
    RooRealVar offset("offset","offset",-5.);
    RooAddition Bmass_offset("Bmass_offset","Bmass_offset",RooArgList(Bmass,offset));
    RooExponential f_bkgCombM("f_bkgCombM","f_bkgCombM",Bmass_offset,bkgCombM_c);// exponential decay
    RooRealVar bkgCombL_c1("bkgCombL_c1","c1",readParam(iBin,"bkgCombL_c1",0),-2.5,2.5);
    RooRealVar bkgCombL_c2("bkgCombL_c2","c2",readParam(iBin,"bkgCombL_c2",0),-2.5,2.5);
    RooRealVar bkgCombL_c3("bkgCombL_c3","c3",readParam(iBin,"bkgCombL_c3",0),-2.5,2.5);
    RooRealVar bkgCombL_c4("bkgCombL_c4","c4",readParam(iBin,"bkgCombL_c4",0),-2.5,2.5);
    RooArgSet f_bkgCombL_argset;
    switch (iBin) {
        case 7:
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1));
            bkgCombL_c2.setVal(0.);
            bkgCombL_c3.setVal(0.);
            bkgCombL_c4.setVal(0.);
            bkgCombL_c2.setConstant(kTRUE);
            bkgCombL_c3.setConstant(kTRUE);
            bkgCombL_c4.setConstant(kTRUE);
            break;
        case 0:
        case 1:
        case 4:
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1,bkgCombL_c2));
            bkgCombL_c3.setVal(0.);
            bkgCombL_c4.setVal(0.);
            bkgCombL_c3.setConstant(kTRUE);
            bkgCombL_c4.setConstant(kTRUE);
            break;
        case 2:
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1,bkgCombL_c2,bkgCombL_c3));
            bkgCombL_c4.setVal(0.);
            bkgCombL_c4.setConstant(kTRUE);
            break;
        case 3:
        case 5:
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1,bkgCombL_c2,bkgCombL_c3,bkgCombL_c4));
            break;
        default:
            bkgCombL_c1.setVal(0.);
            bkgCombL_c2.setVal(0.);
            bkgCombL_c3.setVal(0.);
            bkgCombL_c4.setVal(0.);
            bkgCombL_c1.setConstant(kTRUE);
            bkgCombL_c2.setConstant(kTRUE);
            bkgCombL_c3.setConstant(kTRUE);
            bkgCombL_c4.setConstant(kTRUE);
            break;
    }
    RooPolynomial f_bkgCombL("f_bkgCombL","f_bkgCombL",CosThetaL,f_bkgCombL_argset);
    RooRealVar bkgCombK_c1("bkgCombK_c1","c1",readParam(iBin,"bkgCombK_c1",0),-2.5,2.5);
    RooRealVar bkgCombK_c2("bkgCombK_c2","c2",readParam(iBin,"bkgCombK_c2",0),-2.5,2.5);
    RooRealVar bkgCombK_c3("bkgCombK_c3","c3",readParam(iBin,"bkgCombK_c3",0),-2.5,2.5);
    RooRealVar bkgCombK_c4("bkgCombK_c4","c4",readParam(iBin,"bkgCombK_c4",0),-2.5,2.5);
    RooArgSet f_bkgCombK_argset;
    switch (iBin) {
        case 2:
            f_bkgCombK_argset.add(RooArgSet(bkgCombK_c1));
            bkgCombK_c2.setVal(0.);
            bkgCombK_c3.setVal(0.);
            bkgCombK_c4.setVal(0.);
            bkgCombK_c2.setConstant(kTRUE);
            bkgCombK_c3.setConstant(kTRUE);
            bkgCombK_c4.setConstant(kTRUE);
            break;
        case 0:
        case 1:
            f_bkgCombK_argset.add(RooArgSet(bkgCombK_c1,bkgCombK_c2));
            bkgCombK_c3.setVal(0.);
            bkgCombK_c4.setVal(0.);
            bkgCombK_c3.setConstant(kTRUE);
            bkgCombK_c4.setConstant(kTRUE);
            break;
        case 3:
        case 4:
        case 5:
            f_bkgCombK_argset.add(RooArgSet(bkgCombK_c1,bkgCombK_c2,bkgCombK_c3,bkgCombK_c4));
            break;
        default:
            bkgCombK_c1.setVal(0.);
            bkgCombK_c2.setVal(0.);
            bkgCombK_c3.setVal(0.);
            bkgCombK_c4.setVal(0.);
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
    RooRealVar bkgGauss1_mean1("bkgGauss1_mean1","M_{K*#Mu#Mu}",readParam(iBin,"bkgGauss1_mean1 ",0));
    bkgGauss1_mean1.setError(readParam(iBin,"bkgGauss1_mean1",1));
    RooRealVar bkgGauss1_mean2("bkgGauss1_mean2","M_{K*#Mu#Mu}",readParam(iBin,"bkgGauss1_mean2 ",0));
    bkgGauss1_mean2.setError(readParam(iBin,"bkgGauss1_mean2",1));
    RooRealVar bkgGauss1_sigma1("bkgGauss1_sigma1","#sigma_{11}",readParam(iBin,"bkgGauss1_sigma1 ",0));
    bkgGauss1_sigma1.setError(readParam(iBin,"bkgGauss1_sigma1",1));
    RooRealVar bkgGauss1_sigma2("bkgGauss1_sigma2","#sigma_{12}",readParam(iBin,"bkgGauss1_sigma2 ",0));
    bkgGauss1_sigma2.setError(readParam(iBin,"bkgGauss1_sigma2",1));
    RooRealVar bkgM_frac1("bkgM_frac1","bkgM_frac1",readParam(iBin,"bkgM_frac1",0));
    bkgM_frac1.setError(readParam(iBin,"bkgM_frac1",1));
    RooRealVar bkgGauss2_mean1("bkgGauss2_mean1","M_{K*#Mu#Mu}",readParam(iBin,"bkgGauss2_mean1 ",0));
    bkgGauss2_mean1.setError(readParam(iBin,"bkgGauss2_mean1 ",1));
    RooRealVar bkgGauss2_mean2("bkgGauss2_mean2","M_{K*#Mu#Mu}",readParam(iBin,"bkgGauss2_mean2 ",0));
    bkgGauss2_mean2.setError(readParam(iBin,"bkgGauss2_mean2 ",1));
    RooRealVar bkgGauss2_sigma1("bkgGauss2_sigma1","#sigma_{21}",readParam(iBin,"bkgGauss2_sigma1",0));
    bkgGauss2_sigma1.setError(readParam(iBin,"bkgGauss2_sigma1",1));
    RooRealVar bkgGauss2_sigma2("bkgGauss2_sigma2","#sigma_{22}",readParam(iBin,"bkgGauss2_sigma2",0));
    bkgGauss2_sigma2.setError(readParam(iBin,"bkgGauss2_sigma2",1));
    RooRealVar bkgM_frac2("bkgM_frac2","bkgM_frac2",readParam(iBin,"bkgM_frac2",0));
    bkgM_frac2.setError(readParam(iBin,"bkgM_frac2",1));
    RooRealVar bkgM_frac12("bkgM_frac12","bkgM_frac12",readParam(iBin,"bkgM_frac12",0));
    bkgM_frac12.setError(readParam(iBin,"bkgM_frac12",1));
    RooGaussian f_bkgPeakMGauss11("f_bkgPeakMGauss11","f_bkgPeakMGauss11", Bmass, bkgGauss1_mean1, bkgGauss1_sigma1);
    RooGaussian f_bkgPeakMGauss12("f_bkgPeakMGauss12","f_bkgPeakMGauss12", Bmass, bkgGauss1_mean2, bkgGauss1_sigma2);
    RooGaussian f_bkgPeakMGauss21("f_bkgPeakMGauss21","f_bkgPeakMGauss21", Bmass, bkgGauss2_mean1, bkgGauss2_sigma1);
    RooGaussian f_bkgPeakMGauss22("f_bkgPeakMGauss22","f_bkgPeakMGauss22", Bmass, bkgGauss2_mean2, bkgGauss2_sigma2);
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
        case 3:
        case 5:
            bkgPeakK_c1.setConstant(kFALSE);
            bkgPeakL_c1.setConstant(kFALSE);
            bkgPeakK_c2.setConstant(kFALSE);
            bkgPeakL_c2.setConstant(kFALSE);
            bkgPeakK_c3.setConstant(kFALSE);
            bkgPeakL_c3.setConstant(kFALSE);
            bkgPeakK_c4.setConstant(kFALSE);
            bkgPeakL_c4.setConstant(kFALSE);
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
    RooPolynomial f_bkgPeakL("f_bkgPeakL","f_bkgPeakL",CosThetaL,f_bkgPeakL_argset);
    RooPolynomial f_bkgPeakK("f_bkgPeakK","f_bkgPeakK",CosThetaK,f_bkgPeakK_argset);
    RooProdPdf f_bkgPeakA("f_bkgPeakA", "f_bckPeakA",f_bkgPeakK,f_bkgPeakL);
    RooProdPdf f_bkgPeak("f_bkgPeak", "f_bkgPeak",f_bkgPeakA,f_bkgPeakM12);
    printf("INFO: f_bkgPeak prepared.\n");

    // Observed spectrum = model*fullEfficiency
    RooRealVar nsig("nsig","nsig",10,0,5E3);
    RooRealVar nbkgComb("nbkgComb","nbkgComb",20,0,1E4);
    //RooRealVar nbkgPeak("nbkgPeak","nbkgPeak",readParam(iBin,"nbkgPeak ",0));
    //nbkgPeak.setError(readParam(iBin,"nbkgPeak ",1));
    RooRealVar nbkgPeak("nbkgPeak","nbkgPeak",50,0.,1E4);
    if (iBin == 0 || iBin %2 == 1){
        nbkgPeak.setMin(NULL,0.);
        nbkgPeak.setVal(0.);
        nbkgPeak.setConstant(kTRUE);
    }
    //RooAddPdf kernel("kernel","kernel",RooArgList(f_sig,f_bkgComb,f_bkgPeak),RooArgList(nsig,nbkgComb,nbkgPeak));
    RooAddPdf f("kernel","kernel",RooArgList(f_bkgComb,f_bkgPeak,f_sig),RooArgList(nbkgComb,nbkgPeak,nsig));// test without penalty term

    // Extra penalty term to confine As, Fs, Fl, Afb.
    //RooRealVar t_penalty("t_penalty","t",0.01);
    //RooGenericPdf f_penaltyAfb("f_penaltyAfb","(1-TMath::Erf((afb-0.75*(1-fl))/(1.5*t_penalty*(1-fl))))*(1-TMath::Erf((-afb-0.75*(1-fl))/(1.5*t_penalty*(1-fl))))",RooArgSet(afb,fl,t_penalty));
    //RooGenericPdf f_penaltyAs("f_penaltyAfb","(1-TMath::Erf((afb-2*(1-fl)/3)/(1.5*t_penalty*(1-fl))))*(1-TMath::Erf((-afb-0.75*(1-fl))/(1.5*t_penalty*(1-fl))))",RooArgSet(afb,fl,t_penalty));
    //RooProdPdf f_penalty("f_penalty","f_penalty",f_penaltyAfb,f_penaltyAs);
    //RooProdPdf f("f","f",f_model,f_penalty);
    printf("INFO: f_penalty NOT prepared.\n");

    // Gaussian constraints
    RooGaussian gaus_sigGauss1_sigma("gaus_sigGauss1_sigma","gaus_sigGauss1_sigma",sigGauss1_sigma,RooConst(readParam(iBin,"sigGauss1_sigma ",0)),RooConst(readParam(iBin,"sigGauss1_sigma",1)));
    RooGaussian gaus_sigGauss2_sigma("gaus_sigGauss2_sigma","gaus_sigGauss2_sigma",sigGauss2_sigma,RooConst(readParam(iBin,"sigGauss2_sigma ",0)),RooConst(readParam(iBin,"sigGauss2_sigma",1)));
    RooGaussian gaus_sigM_frac("gaus_sigM_frac","gaus_sigM_frac",sigM_frac,RooConst(readParam(iBin,"sigM_frac ",0)),RooConst(readParam(iBin,"sigM_frac ",1)));

    RooGaussian gaus_nbkgPeak("gaus_nbkgPeak","gaus_nbkgPeak",nbkgPeak,RooConst(readParam(iBin,"nbkgPeak ",0)),RooConst(readParam(iBin,"nbkgPeak ",1)));
    RooGaussian gaus_bkgGauss1_mean1("gaus_bkgGauss1_mean1","gaus_bkgGauss1_mean1",bkgGauss1_mean1,RooConst(readParam(iBin,"bkgGauss1_mean1 ",0)),RooConst(readParam(iBin,"bkgGauss1_mean1 ",1)));
    RooGaussian gaus_bkgGauss1_mean2("gaus_bkgGauss1_mean2","gaus_bkgGauss1_mean2",bkgGauss1_mean2,RooConst(readParam(iBin,"bkgGauss1_mean2 ",0)),RooConst(readParam(iBin,"bkgGauss1_mean2 ",1)));
    RooGaussian gaus_bkgGauss1_sigma1("gaus_bkgGauss1_sigma1","gaus_bkgGauss1_sigma1",bkgGauss1_sigma1,RooConst(readParam(iBin,"bkgGauss1_sigma1 ",0)),RooConst(readParam(iBin,"bkgGauss1_sigma1 ",1)));
    RooGaussian gaus_bkgGauss1_sigma2("gaus_bkgGauss1_sigma2","gaus_bkgGauss1_sigma2",bkgGauss1_sigma2,RooConst(readParam(iBin,"bkgGauss1_sigma2 ",0)),RooConst(readParam(iBin,"bkgGauss1_sigma2 ",1)));
    RooGaussian gaus_bkgM_frac1("gaus_bkgM_frac1","gaus_bkgM_frac1",bkgM_frac1,RooConst(readParam(iBin,"bkgM_frac1 ",0)),RooConst(readParam(iBin,"bkgM_frac1 ",1)));
    RooGaussian gaus_bkgGauss2_mean1("gaus_bkgGauss2_mean1","gaus_bkgGauss2_mean1",bkgGauss2_mean1,RooConst(readParam(iBin,"bkgGauss2_mean1 ",0)),RooConst(readParam(iBin,"bkgGauss2_mean1 ",1)));
    RooGaussian gaus_bkgGauss2_mean2("gaus_bkgGauss2_mean2","gaus_bkgGauss2_mean2",bkgGauss2_mean2,RooConst(readParam(iBin,"bkgGauss2_mean2 ",0)),RooConst(readParam(iBin,"bkgGauss2_mean2 ",1)));
    RooGaussian gaus_bkgGauss2_sigma1("gaus_bkgGauss2_sigma1","gaus_bkgGauss2_sigma1",bkgGauss2_sigma1,RooConst(readParam(iBin,"bkgGauss2_sigma1 ",0)),RooConst(readParam(iBin,"bkgGauss2_sigma1 ",1)));
    RooGaussian gaus_bkgGauss2_sigma2("gaus_bkgGauss2_sigma2","gaus_bkgGauss2_sigma2",bkgGauss2_sigma2,RooConst(readParam(iBin,"bkgGauss2_sigma2 ",0)),RooConst(readParam(iBin,"bkgGauss2_sigma2 ",1)));
    RooGaussian gaus_bkgM_frac2("gaus_bkgM_frac2","gaus_bkgM_frac2",bkgM_frac2,RooConst(readParam(iBin,"bkgM_frac2 ",0)),RooConst(readParam(iBin,"bkgM_frac2 ",1)));
    RooGaussian gaus_bkgM_frac12("gaus_bkgM_frac12","gaus_bkgM_frac12",bkgM_frac12,RooConst(readParam(iBin,"bkgM_frac12 ",0)),RooConst(readParam(iBin,"bkgM_frac12 ",1)));
    
    RooGaussian gaus_bkgPeakL_c1("gaus_bkgPeakL_c1","gaus_bkgPeakL_c1",bkgPeakL_c1,RooConst(readParam(iBin,"bkgPeakL_c1 ",0)),RooConst(readParam(iBin,"bkgPeakL_c1 ",1)));
    RooGaussian gaus_bkgPeakL_c2("gaus_bkgPeakL_c2","gaus_bkgPeakL_c2",bkgPeakL_c2,RooConst(readParam(iBin,"bkgPeakL_c2 ",0)),RooConst(readParam(iBin,"bkgPeakL_c2 ",1)));
    RooGaussian gaus_bkgPeakL_c3("gaus_bkgPeakL_c3","gaus_bkgPeakL_c3",bkgPeakL_c3,RooConst(readParam(iBin,"bkgPeakL_c3 ",0)),RooConst(readParam(iBin,"bkgPeakL_c3 ",1)));
    RooGaussian gaus_bkgPeakL_c4("gaus_bkgPeakL_c4","gaus_bkgPeakL_c4",bkgPeakL_c4,RooConst(readParam(iBin,"bkgPeakL_c4 ",0)),RooConst(readParam(iBin,"bkgPeakL_c4 ",1)));
    RooGaussian gaus_bkgPeakK_c1("gaus_bkgPeakK_c1","gaus_bkgPeakK_c1",bkgPeakK_c1,RooConst(readParam(iBin,"bkgPeakK_c1 ",0)),RooConst(readParam(iBin,"bkgPeakK_c1 ",1)));
    RooGaussian gaus_bkgPeakK_c2("gaus_bkgPeakK_c2","gaus_bkgPeakK_c2",bkgPeakK_c2,RooConst(readParam(iBin,"bkgPeakK_c2 ",0)),RooConst(readParam(iBin,"bkgPeakK_c2 ",1)));
    RooGaussian gaus_bkgPeakK_c3("gaus_bkgPeakK_c3","gaus_bkgPeakK_c3",bkgPeakK_c3,RooConst(readParam(iBin,"bkgPeakK_c3 ",0)),RooConst(readParam(iBin,"bkgPeakK_c3 ",1)));
    RooGaussian gaus_bkgPeakK_c4("gaus_bkgPeakK_c4","gaus_bkgPeakK_c4",bkgPeakK_c4,RooConst(readParam(iBin,"bkgPeakK_c4 ",0)),RooConst(readParam(iBin,"bkgPeakK_c4 ",1)));
    //RooBifurGauss gaus_fs("gaus_fs","gaus_fs",fs,RooConst(readParam(iBin,"fs ",0)),RooConst(readParam(iBin,"fs ",1)),RooConst(readParam(iBin,"fs ",2)));
    //RooBifurGauss gaus_as("gaus_as","gaus_as",as,RooConst(readParam(iBin,"as ",0)),RooConst(readParam(iBin,"as ",1)),RooConst(readParam(iBin,"as ",2)));
    RooBifurGauss gaus_fs("gaus_fs","gaus_fs",fs,RooConst(0.0129254),RooConst(0.00898344),RooConst(0.0101371));// 2011 result
    RooBifurGauss gaus_as("gaus_as","gaus_as",as,RooConst(-0.0975919),RooConst(0.00490805),RooConst(0.0049092));
    
    RooArgSet gausConstraints(gaus_sigGauss1_sigma,gaus_sigGauss2_sigma,gaus_sigM_frac);
    switch (iBin) {
        case 2:
            //1 double guassian ,4+4 deg. ploy
            gausConstraints.add(RooArgSet(gaus_bkgPeakL_c1,gaus_bkgPeakL_c2,gaus_bkgPeakL_c3,gaus_bkgPeakL_c4));
            gausConstraints.add(RooArgSet(gaus_bkgPeakK_c1,gaus_bkgPeakK_c2,gaus_bkgPeakK_c3,gaus_bkgPeakK_c4));
            gausConstraints.add(RooArgSet(gaus_bkgGauss1_mean1,gaus_bkgGauss1_mean2,gaus_bkgGauss1_sigma1,gaus_bkgGauss1_sigma2,gaus_bkgM_frac1));
            gausConstraints.add(gaus_nbkgPeak);
            break;
        case 4:
            //2 double guassian ,4+4 deg. ploy
            gausConstraints.add(RooArgSet(gaus_bkgPeakL_c1,gaus_bkgPeakL_c2,gaus_bkgPeakL_c3,gaus_bkgPeakL_c4));
            gausConstraints.add(RooArgSet(gaus_bkgPeakK_c1,gaus_bkgPeakK_c2,gaus_bkgPeakK_c3,gaus_bkgPeakK_c4));
            gausConstraints.add(RooArgSet(gaus_bkgGauss1_mean1,gaus_bkgGauss1_mean2,gaus_bkgGauss1_sigma1,gaus_bkgGauss1_sigma2,gaus_bkgM_frac1));
            gausConstraints.add(RooArgSet(gaus_bkgGauss2_mean1,gaus_bkgGauss2_mean2,gaus_bkgGauss2_sigma1,gaus_bkgGauss2_sigma2,gaus_bkgM_frac2));
            gausConstraints.add(gaus_bkgM_frac12);
            gausConstraints.add(gaus_nbkgPeak);
            break;
        case 6:
            //1 guassian ,2+2 deg. ploy
            gausConstraints.add(RooArgSet(gaus_bkgPeakL_c1,gaus_bkgPeakL_c2));
            gausConstraints.add(RooArgSet(gaus_bkgPeakK_c1,gaus_bkgPeakK_c2));
            gausConstraints.add(RooArgSet(gaus_bkgGauss2_mean1,gaus_bkgGauss2_sigma1));
            gausConstraints.add(gaus_nbkgPeak);
            break;
    }
    if (iBin == 3 || iBin == 5) gausConstraints.add(RooArgSet(gaus_fs,gaus_as));
    
    // Get data and apply unbinned fit
    RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Bmass,CosThetaK,CosThetaL,Q2),q2range[iBin],0);
    f.fitTo(*data,Extended(kTRUE),ExternalConstraints(gausConstraints),Minimizer("Minuit"));

    // Draw the frame on the canvas
    TCanvas* c = new TCanvas("c");
    RooPlot* framemass = Bmass.frame();
    data->plotOn(framemass,Binning(20));
    f.plotOn(framemass,LineColor(1));
    f.plotOn(framemass,Components(f_sig),LineColor(4),LineWidth(2));
    f.plotOn(framemass,Components(f_bkgComb),LineColor(2),LineWidth(2),LineStyle(2));
    f.plotOn(framemass,Components(f_bkgPeak),LineColor(6),LineWidth(2),LineStyle(2));

    framemass->SetTitle("");
    framemass->SetMinimum(0);
    framemass->Draw();
    
    TLatex *t1 = new TLatex();
    t1->SetNDC();
    double fixNDC = 0;
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    c->Print(TString::Format("./plots/%s_bin%d.pdf",outfile,iBin));

    // Draw projection to CosThetaK
    RooPlot* framecosk = CosThetaK.frame(); 
    data->plotOn(framecosk,Binning(20)); 
    f.plotOn(framecosk,LineColor(1)); 
    f.plotOn(framecosk,Components(f_sig),LineColor(4),LineWidth(2));
    f.plotOn(framecosk,Components(f_bkgComb),LineColor(2),LineWidth(2),LineStyle(2));
    f.plotOn(framecosk,Components(f_bkgPeak),LineColor(6),LineWidth(2),LineStyle(2));

    framecosk->SetTitle("");
    framecosk->SetMinimum(0);
    framecosk->Draw();

    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    c->Update();
    c->Print(TString::Format("./plots/%s_cosk_bin%d.pdf",outfile,iBin));

    // Draw projection to CosThetaL
    RooPlot* framecosl = CosThetaL.frame(); 
    data->plotOn(framecosl,Binning(20)); 
    f.plotOn(framecosl,LineColor(1)); 
    f.plotOn(framecosl,Components(f_sig),LineColor(4),LineWidth(2));
    f.plotOn(framecosl,Components(f_bkgComb),LineColor(2),LineWidth(2),LineStyle(2));
    f.plotOn(framecosl,Components(f_bkgPeak),LineColor(6),LineWidth(2),LineStyle(2));

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
    double val[3]={0,0,0};
    val[0] = fl.getVal();val[1] = fl.getError();
    writeParam(iBin, "fl", val);
    val[0] = afb.getVal();val[1] = afb.getError();
    writeParam(iBin, "afb",val);
    val[0] = fs.getVal();val[1] = fs.getErrorHi();val[2] = fs.getErrorLo();
    writeParam(iBin, "fs", val, 3);
    val[0] = as.getVal();val[1] = as.getErrorHi();val[2] = fs.getErrorLo();
    writeParam(iBin, "as", val, 3);

    std::vector<double> output;
    output.push_back(fl.getVal());
    output.push_back(fl.getError());
    output.push_back(afb.getVal());
    output.push_back(afb.getError());
    output.push_back(fs.getVal());
    output.push_back(fs.getError());
    output.push_back(as.getVal());
    output.push_back(as.getError());
    return output;
}//}}}

void angular(const char outfile[] = "angular", bool doFit = true)
{//{{{

    double x[8]={1.5,3.15,6.49,9.385,11.475,13.52,15.09,17.5};
    double xerr[8]={0.5,1.15,2.09,0.705,1.385,0.66,0.91,1.5};
    double yafb[8],yerrafb[8],yfl[8],yerrfl[8];

    if (doFit){
        angular3D_bin(3);
        angular3D_bin(5);

        angular3D_bin(0);
        angular3D_bin(1);
        angular3D_bin(2);
        angular3D_bin(4);
        angular3D_bin(6);
        angular3D_bin(7);
    }

    // Checkout input data
    for(int ibin = 0; ibin < 8; ibin++){
        yfl[ibin]       = readParam(ibin,"fl ",0);
        yerrfl[ibin]    = readParam(ibin,"fl ",1);
        yafb[ibin]      = readParam(ibin,"afb ",0);
        yerrafb[ibin]   = readParam(ibin,"afb ",1);
        printf("yafb[%d]=%6.4f +- %6.4f\n",ibin,yafb[ibin],yerrafb[ibin]);
        printf("yfl [%d]=%6.4f +- %6.4f\n",ibin,yfl[ibin],yerrfl[ibin]);
    }
    
    // Draw
    TCanvas *c = new TCanvas();
    TH2F *frame = new TH2F("frame","",18,1,19,10,-1,1);
    frame->SetStats(kFALSE);

    frame->SetXTitle("q^{2} [(GeV)^{2}]");
    frame->SetYTitle("F_{L}");
    frame->SetAxisRange(0,1,"Y");
    frame->Draw();
    TGraphAsymmErrors *g_fl  = new TGraphAsymmErrors(8,x,yfl,xerr,xerr,yerrfl,yerrfl);
    g_fl->SetFillColor(2);
    g_fl->SetFillStyle(3001);
    g_fl->Draw("P2");
    c->Print(TString::Format("./plots/%s_fl.pdf",outfile));
    c->Clear();

    frame->SetTitle("");
    frame->SetXTitle("q^{2} [(GeV)^{2}]");
    frame->SetYTitle("A_{FB}");
    frame->SetAxisRange(-1,1,"Y");
    frame->Draw();
    TGraphAsymmErrors *g_afb = new TGraphAsymmErrors(8,x,yafb,xerr,xerr,yerrafb,yerrafb);
    g_afb->SetFillColor(2);
    g_afb->SetFillStyle(3001);
    g_afb->Draw("P2");
    c->Print(TString::Format("./plots/%s_afb.pdf",outfile));
}//}}}

//_________________________________________________________________________________
//_________________________________________________________________________________
int main(int argc, char** argv) {
    // Tags
    is7TeVCheck = false;   
    
    if (argc <= 2) {
        printf("Usage       : ./fit Function infile\n");
        printf("Functions   :\n");
        printf("    bmass               Fit to mass spectrum using a double Gaussian signal and Chebyshev bkg.\n");
        printf("    fl_gen              Derive F_{L} from cosThetaK distribution at GEN level.\n");
        printf("    angular_gen         Derive F_{L} and A_{FB} from angular distribution at GEN level.\n");
        printf("    acceptance          Get 2D acceptance map from unfiltered signal GEN, |Mu pT| > 2.5 GeV, |Mu eta| < 2.3.\n");
        printf("    efficiency          Get 2D reconstruction efficiency map from signal simulation.\n");
        printf("    angular2D           Same as angular_gen, but fit to data with efficiency correction, bkg component is NOT considered.\n");
        printf("    angular3D_1a_Sm     Leading step1 to angular3D, determine signal shape from simulation.\n");
        printf("    angular3D_1b_YpPm   Leading step2 to angular3D, determine mass spectrum of peaking bkg from simulation.\n");
        printf("    angular3D_2a_PkPl   Leading step3 to angular3D, determine angular dist. of peaking bkg from simulation.\n");
        printf("    angular3D_prior     Leading step4 to angular3D, fit to data sideband to get initial values of combinatorial bkg.\n");
        printf("    angular3D           Derive F_{L} and A_{FB} by fitting to mass and angular distribution.\n");
        printf("Remark      :Outputs will be stored in ./plots, please keep the directory.\n");
        return 0;
    }

    // main
    TString func    = argv[1];
    TString infile  = argv[2];
  
    if (func == "bmass") {
        if (argc != 4){
            printf("./fit bmass infile binID\n");
            for (int i = 0; i < 10; i++) {
                printf("    Bin %d : %s\n",i,q2range[i]);
            }
            return 0;
        }
        int iBin = atoi(argv[3]);

        ch->Add(infile.Data());
        if (ch == NULL) gSystem->Exit(0);
        const char outfile[]="bmass";
        bmass(iBin, outfile); 
    }else if (func == "fl_gen"){
        ch->Add(infile.Data());
        if (ch == NULL) gSystem->Exit(0);
        const char outfile[]="fl_gen";
        fl_gen(outfile);
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
    }else if (func == "angular2D"){
        ch->Add(infile.Data());
        if (ch == NULL) gSystem->Exit(0);
        const char outfile[]="angular2D";
        for (int iBin = 0; iBin < 8; iBin++) {
            if (iBin == 3 || iBin == 5) continue;
            angular2D_bin(iBin);
        }
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
            fx(iBin,func,true);
        }
    }else if (func == "angular3D"){
        ch->Add(infile.Data());
        if (ch == NULL) gSystem->Exit(0);
        const char outfile[]="angular3D";
        angular(outfile, false);
        //angular3D_bin(3);
        //angular3D_bin(5);
    }else if (func == "test"){
        ch->Add(infile.Data());
        if (ch == NULL) gSystem->Exit(0);
        const char outfile[]="test";
        for (int i = 0; i < 8; i++) {
            // Test whatever you want!
            drawEfficiency(i);
        }
    }else{ 
        cerr << "No function available for: " << func.Data() << endl; 
    }
    printf("%lld entries processed.\n",ch->GetEntries());
    gSystem->Exit(0);

    return 0 ;
}
