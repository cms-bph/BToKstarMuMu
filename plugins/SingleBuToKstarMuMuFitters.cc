// vim: sw=4 ts=4 fdm=marker et:

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
#include <RooGaussian.h>
#include <RooChebychev.h> 
#include <RooAddPdf.h>
#include <RooDataSet.h>
#include <RooAbsData.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooGenericPdf.h> 
#include <RooPolynomial.h>
#include <RooProdPdf.h>
#include <RooDataHist.h>

#include "tools.h" 

using namespace std; 
using namespace RooFit;

TChain *ch=new TChain("tree");
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

TF2  *f2_fcn = NULL;
double model_2D(double x, double y,int nPar, double* par)
{//{{{
    for (int i = 0; i < nPar; i++) f2_fcn->SetParameter(i,par[i]);
    return f2_fcn->Eval(x,y);
}//}}}

TH2F *h2_fcn = NULL;
// nParameters, ???, return fcn value, parameter array, strategy
void fcn_binnedChi2_2D(int &npar, double *gin, double &f, double *par, int iflag)
{//{{{
    f=0;
    for (int i = 1; i <= h2_fcn->GetNbinsX(); i++) {
        for (int j = 1; j <= h2_fcn->GetNbinsY(); j++) {
            int gBin = h2_fcn->GetBin(i,j);
            double xx = h2_fcn->GetXaxis()->GetBinCenter(i);
            double yy = h2_fcn->GetYaxis()->GetBinCenter(j);
            double measure  = h2_fcn->GetBinContent(gBin);
            double error    = h2_fcn->GetBinError(gBin);
            double func     = model_2D(xx,yy,npar, par);
            double delta    = (measure-func)/error;
            f+=delta*delta;
        }
    }
    //printf("FCN in calls = %f\n",f);
}//}}}


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
    //f.plotOn(framecosk,Components(f_sig),LineColor(4),LineWidth(2));
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

std::vector<double> acceptance(int iBin) // acceptance
{//{{{
    std::vector<double> output;

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
    
    TF2 f2_model("f2_model","([0]+[1]*y+[2]*y**2+[3]*y**3)+([4]+[5]*y+[6]*y**2+[7]*y**3)*x**2+([8]+[9]*y+[10]*y**2+[11]*y**3)*x**3+([12]+[13]*y+[14]*y**2+[15]*y**3)*x**4+([16]+[17]*y+[18]*y**2+[19]*y**3)*x**6",-1.,1.,-1.,1.);
    f2_fcn = &f2_model;
    gMinuit->DefineParameter( 0, "k0l0", .01,    .1,    -1E-2, 1E-1);
    gMinuit->DefineParameter( 1, "k1l0",  0.,    .1,    -1E-2, 1E-2);
    gMinuit->DefineParameter( 2, "k2l0",  0.,    .1,    -1E-2, 1E-2);
    gMinuit->DefineParameter( 3, "k3l0",  0.,    .1,    -1E-2, 1E-2);
    gMinuit->DefineParameter( 4, "k0l2", .01,    .1,    -1E-2, 1E-1);
    gMinuit->DefineParameter( 5, "k1l2",  0.,    .1,    -1E-2, 1E-2);
    gMinuit->DefineParameter( 6, "k2l2",  0.,    .1,    -1E-2, 1E-2);
    gMinuit->DefineParameter( 7, "k3l2",  0.,    .1,    -1E-2, 1E-2);
    gMinuit->DefineParameter( 8, "k0l3",  0.,    .1,    -1E-2, 1E-2);
    gMinuit->DefineParameter( 9, "k1l3",  0.,    .1,    -1E-2, 1E-2);
    gMinuit->DefineParameter(10, "k2l3",  0.,    .1,    -1E-2, 1E-2);
    gMinuit->DefineParameter(11, "k3l3",  0.,    .1,    -1E-2, 1E-2);
    gMinuit->DefineParameter(12, "k0l4",1E-3,  1E-3,    -1   , 1   );
    gMinuit->DefineParameter(13, "k1l4",  0.,    .1,    -1E-1, 1E-1);
    gMinuit->DefineParameter(14, "k2l4",  0.,    .1,    -1E-1, 1E-1);
    gMinuit->DefineParameter(15, "k3l4",  0.,    .1,    -1E-1, 1E-1);
    gMinuit->DefineParameter(16, "k0l6",1E-3,  1E-3,    -1   , 1   );
    gMinuit->DefineParameter(17, "k1l6",  0.,    .1,    -1E-1, 1E-1);
    gMinuit->DefineParameter(18, "k2l6",  0.,    .1,    -1E-1, 1E-1);
    gMinuit->DefineParameter(19, "k3l6",  0.,    .1,    -1E-1, 1E-1);
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
    gMinuit->Command("MINI");
    gMinuit->Command("MINI");
    gMinuit->Command("MINI");
    gMinuit->Command("MINI");
    //gMinuit->Command("HESSE");
    gMinuit->Command("IMPROVE");
    gMinuit->Command("MINOS");

    double arrPar[nPar];
    double arrParErr[nPar];
    for (int iPar = 0; iPar < nPar; iPar++) gMinuit->GetParameter(iPar,arrPar[iPar],arrParErr[iPar]);
    for (int iPar = 0; iPar < nPar; iPar++) f2_model.SetParameter(iPar,arrPar[iPar]);

    // Calculate chi2
    double chi2Val=0;
    for (int i = 1; i <= 6; i++) {//thetaL
        for (int j = 1; j <= 5; j++) {//thetaK
            chi2Val += pow((f2_model.Eval(h2_acc.GetXaxis()->GetBinCenter(i),h2_acc.GetYaxis()->GetBinCenter(j))-h2_acc.GetBinContent(i,j))/h2_acc.GetBinError(i,j),2);
            printf("FitValue= (%9.6f - %9.6f)/%f at i=%2d(%.2f),j=%2d(%.2f)\n",f2_model.Eval(h2_acc.GetXaxis()->GetBinCenter(i),h2_acc.GetYaxis()->GetBinCenter(j)),h2_acc.GetBinContent(i,j),h2_acc.GetBinError(i,j),i,h2_acc.GetXaxis()->GetBinCenter(i),j,h2_acc.GetYaxis()->GetBinCenter(j));
        }
    }
    printf("Chi2=%f \n",chi2Val);
    
    // Prepare draw
    TCanvas canvas("canvas");
    TLatex *latex = new TLatex();
    
    // Draw efficiency
    h2_acc.SetStats(0);
    h2_acc.SetMaximum(.02);
    h2_acc.Draw("LEGO2 Z");
    latex->DrawLatexNDC(0.35,0.95,TString::Format("Acceptance_{RECO} in Bin%d",iBin));
    
    // Draw FitResult
    f2_model.SetTitle("");
    f2_model.SetMaximum(.02);
    f2_model.SetLineWidth(1);
    f2_model.Draw("SURF3 SAME ");
    canvas.Print(TString::Format("./plots/acceptance_2D_bin%d.pdf",iBin));
    

    //// Draw compare
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

    chi2Val = 0;
    TF1  *f_cosk = new TF1("f_cosk"
                          ,"([0]+[1]*x+[2]*x**2+[3]*x**3)+([4]+[5]*x+[6]*x**2+[7]*x**3)/3+([8]+[9]*x+[10]*x**2+[11]*x**3)*0+([12]+[13]*x+[14]*x**2+[15]*x**3)/5+([16]+[17]*x+[18]*x**2+[19]*x**3)/7"
                          ,-1.,1.);
    for (int iPar = 0; iPar < nPar; iPar++) f_cosk->SetParameter(iPar,arrPar[iPar]);
    f_cosk->Draw("SAME");
    for (int kBin = 1; kBin <= 5; kBin++) {
        chi2Val += pow((f_cosk->Eval(h_cosk->GetBinCenter(kBin))-h_cosk->GetBinContent(kBin))/h_cosk->GetBinError(kBin),2);
    }
    latex->DrawLatexNDC(0.01,0.95,TString::Format("#chi^{2} = %f",chi2Val));
    canvas.Update();
    canvas.Print(TString::Format("./plots/recoEff_cosK_bin%d.pdf",iBin));
    
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

    chi2Val = 0;
    TF1  *f_cosl = new TF1("f_cosl"
                          ,"([0]+[1]*x*0+[2]/3+[3]*0)+([4]+[5]*0+[6]/3+[7]*0)*x**2+([8]+[9]*0+[10]/3+[11]*0)*x**3+([12]+[13]*0+[14]/3+[15]*0)*x**4+([16]+[17]*0+[18]/3+[19]*0)*x**6"
                          ,-1.,1.);
    for (int iPar = 0; iPar < nPar; iPar++) f_cosl->SetParameter(iPar,arrPar[iPar]);
    f_cosl->Draw("SAME");
    for (int lBin = 1; lBin <= 6; lBin++) {
        chi2Val += pow((f_cosl->Eval(h_cosl->GetBinCenter(lBin))-h_cosl->GetBinContent(lBin))/h_cosl->GetBinError(lBin),2);
    }
    latex->DrawLatexNDC(0.01,0.95,TString::Format("#chi^{2} = %f",chi2Val));
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
    std::vector<double> output;
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

    // Using RooFit //{{{
    /*
    RooDataHist *data  = new RooDataHist("data","data",RooArgList(genCosThetaL,genCosThetaK),Import(h2_rec,kFALSE)) ; 
    
    // Define fit function
    RooRealVar k0l0("k0l0","k0l0",1,-1E2,1E2);// As the normalizer, forced to be 1.
    RooRealVar k1l0("k1l0","k1l0",0,-1E2,1E2);
    RooRealVar k2l0("k2l0","k2l0",0,-1E2,1E2);
    RooRealVar k3l0("k3l0","k3l0",0,-1E2,1E2);
    RooRealVar k0l2("k0l2","k0l2",1,-1E2,1E2);
    RooRealVar k1l2("k1l2","k1l2",0,-1E2,1E2);
    RooRealVar k2l2("k2l2","k2l2",0,-1E2,1E4);
    RooRealVar k3l2("k3l2","k3l2",0,-1E2,1E2);
    RooRealVar k0l3("k0l3","k0l3",0,-1E2,1E2);
    RooRealVar k1l3("k1l3","k1l3",0,-1E2,1E2);
    RooRealVar k2l3("k2l3","k2l3",0,-1E2,1E2);
    RooRealVar k3l3("k3l3","k3l3",0,-1E2,1E2);
    RooRealVar k0l4("k0l4","k0l4",1,-1E2,1E2);
    RooRealVar k1l4("k1l4","k1l4",0,-1E2,1E2);
    RooRealVar k2l4("k2l4","k2l4",0,-1E2,1E2);
    RooRealVar k3l4("k3l4","k3l4",0,-1E2,1E2);
    RooRealVar k0l6("k0l6","k0l6",1,-1E2,1E2);
    RooRealVar k1l6("k1l6","k1l6",0,-1E2,1E2);
    RooRealVar k2l6("k2l6","k2l6",0,-1E2,1E2);
    RooRealVar k3l6("k3l6","k3l6",0,-1E2,1E2);
    RooAddPdf *f_addition = 0;
    //RooGenericPdf *f_kl0 = new RooGenericPdf("f_kl0","f_kl0","(1+k1l0*genCosThetaK+k2l0*genCosThetaK**2+k3l0*genCosThetaK**3)/(2+2*k2l0/3)"                ,RooArgList(genCosThetaK,genCosThetaL,k1l0,k2l0,k3l0));
    RooGenericPdf *f_kl0 = new RooGenericPdf("f_kl0","f_kl0","(1+k1l0*genCosThetaK+k2l0*genCosThetaK**2+k3l0*genCosThetaK**3)/(2+2/3*k2l0)*0.5+(k0l3+k1l3*genCosThetaK+k2l3*genCosThetaK**2+k3l3*genCosThetaK**3)*genCosThetaL**3",RooArgList(genCosThetaK,genCosThetaL,k1l0,k2l0,k3l0,k0l3,k1l3,k2l3,k3l3));
    RooGenericPdf *f_kl2 = new RooGenericPdf("f_kl2","f_kl2","(1+k1l2*genCosThetaK+k2l2*genCosThetaK**2+k3l2*genCosThetaK**3)/(2+2/3*k2l2)*1.5*genCosThetaL**2",RooArgList(genCosThetaK,genCosThetaL,k1l2,k2l2,k3l2));
    //RooGenericPdf *f_kl3 = new RooGenericPdf("f_kl3","f_kl3","(1+k1l3*genCosThetaK+k2l3*genCosThetaK**2+k3l3*genCosThetaK**3)/(2+2/3*k2l2)*genCosThetaL**3",RooArgList(genCosThetaK,genCosThetaL,k1l3,k2l3,k3l3));
    RooGenericPdf *f_kl4 = new RooGenericPdf("f_kl4","f_kl4","(1+k1l4*genCosThetaK+k2l4*genCosThetaK**2+k3l4*genCosThetaK**3)/(2+2/3*k2l4)*2.5*genCosThetaL**4",RooArgList(genCosThetaK,genCosThetaL,k1l4,k2l4,k3l4));
    RooGenericPdf *f_kl6 = new RooGenericPdf("f_kl6","f_kl6","(1+k1l6*genCosThetaK+k2l6*genCosThetaK**2+k3l6*genCosThetaK**3)/(2+2/3*k2l6)*3.5*genCosThetaL**6",RooArgList(genCosThetaK,genCosThetaL,k1l6,k2l6,k3l6));
    if (iBin == 0){
        f_addition = new RooAddPdf("f","f",RooArgList(*f_kl0,*f_kl2,*f_kl4,*f_kl6),RooArgList(k0l0,k0l2,k0l4,k0l6));
    }else if (iBin == 1){
        f_addition = new RooAddPdf("f","f",RooArgList(*f_kl0,*f_kl2,*f_kl4),RooArgList(k0l0,k0l2,k0l4));
        k0l3.setVal(0.);
        k1l3.setVal(0.);
        k2l3.setVal(0.);
        k3l3.setVal(0.);
        k0l3.setConstant(kTRUE);
        k1l3.setConstant(kTRUE);
        k2l3.setConstant(kTRUE);
        k3l3.setConstant(kTRUE);
    }else if (iBin >=1 && iBin < 6){
        f_addition = new RooAddPdf("f","f",RooArgList(*f_kl0,*f_kl2,*f_kl4),RooArgList(k0l0,k0l2,k0l4));
    }else{
        f_addition = new RooAddPdf("f","f",RooArgList(*f_kl0,*f_kl2),RooArgList(k0l0,k0l2));
    }
    RooAddPdf *f = f_addition;
    
    // Fit reco efficiency
    printf("Start fitting...\n");
    //f->chi2FitTo(*data,Extended());
    //f->chi2FitTo(*data,Extended());
    RooFitResult *f_fitresult = f->chi2FitTo(*data,Save(kTRUE),Extended());
    printf("End fitting...\n");
    */ //}}}

    // Using pure TMinuit
    int nPar = 20;
    TMinuit *gMinuit = new TMinuit(nPar);
    h2_fcn = &h2_rec;
    gMinuit->SetFCN(fcn_binnedChi2_2D);
    
    TF2 f2_model("f2_model","([0]+[1]*y+[2]*y**2+[3]*y**3)+([4]+[5]*y+[6]*y**2+[7]*y**3)*x**2+([8]+[9]*y+[10]*y**2+[11]*y**3)*x**3+([12]+[13]*y+[14]*y**2+[15]*y**3)*x**4+([16]+[17]*y+[18]*y**2+[19]*y**3)*x**6",-1.,1.,-1.,1.);
    f2_fcn = &f2_model;
    gMinuit->DefineParameter( 0, "k0l0", .01,    .1,    -1E-2, 1E-1);
    gMinuit->DefineParameter( 1, "k1l0",  0.,    .1,    -1E-2, 1E-2);
    gMinuit->DefineParameter( 2, "k2l0",  0.,    .1,    -1E-2, 1E-2);
    gMinuit->DefineParameter( 3, "k3l0",  0.,    .1,    -1E-2, 1E-2);
    gMinuit->DefineParameter( 4, "k0l2", .01,    .1,    -1E-2, 1E-1);
    gMinuit->DefineParameter( 5, "k1l2",  0.,    .1,    -1E-2, 1E-2);
    gMinuit->DefineParameter( 6, "k2l2",  0.,    .1,    -1E-2, 1E-2);
    gMinuit->DefineParameter( 7, "k3l2",  0.,    .1,    -1E-2, 1E-2);
    gMinuit->DefineParameter( 8, "k0l3",  0.,    .1,    -1E-2, 1E-2);
    gMinuit->DefineParameter( 9, "k1l3",  0.,    .1,    -1E-2, 1E-2);
    gMinuit->DefineParameter(10, "k2l3",  0.,    .1,    -1E-2, 1E-2);
    gMinuit->DefineParameter(11, "k3l3",  0.,    .1,    -1E-2, 1E-2);
    gMinuit->DefineParameter(12, "k0l4",1E-3,  1E-3,    -1   , 1   );
    gMinuit->DefineParameter(13, "k1l4",  0.,    .1,    -1E-1, 1E-1);
    gMinuit->DefineParameter(14, "k2l4",  0.,    .1,    -1E-1, 1E-1);
    gMinuit->DefineParameter(15, "k3l4",  0.,    .1,    -1E-1, 1E-1);
    gMinuit->DefineParameter(16, "k0l6",1E-3,  1E-3,    -1   , 1   );
    gMinuit->DefineParameter(17, "k1l6",  0.,    .1,    -1E-1, 1E-1);
    gMinuit->DefineParameter(18, "k2l6",  0.,    .1,    -1E-1, 1E-1);
    gMinuit->DefineParameter(19, "k3l6",  0.,    .1,    -1E-1, 1E-1);
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
    gMinuit->Command("MINI");
    gMinuit->Command("MINI");
    gMinuit->Command("MINI");
    gMinuit->Command("MINI");
    //gMinuit->Command("HESSE");
    gMinuit->Command("IMPROVE");
    gMinuit->Command("MINOS");

    double arrPar[nPar];
    double arrParErr[nPar];
    for (int iPar = 0; iPar < nPar; iPar++) gMinuit->GetParameter(iPar,arrPar[iPar],arrParErr[iPar]);
    for (int iPar = 0; iPar < nPar; iPar++) f2_model.SetParameter(iPar,arrPar[iPar]);

    // Calculate chi2
    double chi2Val=0;
    for (int i = 1; i <= 6; i++) {//thetaL
        for (int j = 1; j <= 5; j++) {//thetaK
            chi2Val += pow((f2_model.Eval(h2_rec.GetXaxis()->GetBinCenter(i),h2_rec.GetYaxis()->GetBinCenter(j))-h2_rec.GetBinContent(i,j))/h2_rec.GetBinError(i,j),2);
            printf("FitValue= (%9.6f - %9.6f)/%f at i=%2d(%.2f),j=%2d(%.2f)\n",f2_model.Eval(h2_rec.GetXaxis()->GetBinCenter(i),h2_rec.GetYaxis()->GetBinCenter(j)),h2_rec.GetBinContent(i,j),h2_rec.GetBinError(i,j),i,h2_rec.GetXaxis()->GetBinCenter(i),j,h2_rec.GetYaxis()->GetBinCenter(j));
        }
    }
    printf("Chi2=%f \n",chi2Val);
    
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
    latex->DrawLatexNDC(0.3,0.95,TString::Format("#varepsilon_{RECO,measured} / #varepsilon_{RECO,fit} in Bin%d",iBin));
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

    chi2Val = 0;
    TF1  *f_cosk = new TF1("f_cosk"
                          ,"([0]+[1]*x+[2]*x**2+[3]*x**3)+([4]+[5]*x+[6]*x**2+[7]*x**3)/3+([8]+[9]*x+[10]*x**2+[11]*x**3)*0+([12]+[13]*x+[14]*x**2+[15]*x**3)/5+([16]+[17]*x+[18]*x**2+[19]*x**3)/7"
                          ,-1.,1.);
    for (int iPar = 0; iPar < nPar; iPar++) f_cosk->SetParameter(iPar,arrPar[iPar]);
    f_cosk->Draw("SAME");
    for (int kBin = 1; kBin <= 5; kBin++) {
        chi2Val += pow((f_cosk->Eval(h_cosk->GetBinCenter(kBin))-h_cosk->GetBinContent(kBin))/h_cosk->GetBinError(kBin),2);
    }
    latex->DrawLatexNDC(0.01,0.95,TString::Format("#chi^{2} = %f",chi2Val));
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

    chi2Val = 0;
    TF1  *f_cosl = new TF1("f_cosl"
                          ,"([0]+[1]*x*0+[2]/3+[3]*0)+([4]+[5]*0+[6]/3+[7]*0)*x**2+([8]+[9]*0+[10]/3+[11]*0)*x**3+([12]+[13]*0+[14]/3+[15]*0)*x**4+([16]+[17]*0+[18]/3+[19]*0)*x**6"
                          ,-1.,1.);
    for (int iPar = 0; iPar < nPar; iPar++) f_cosl->SetParameter(iPar,arrPar[iPar]);
    f_cosl->Draw("SAME");
    for (int lBin = 1; lBin <= 6; lBin++) {
        chi2Val += pow((f_cosl->Eval(h_cosl->GetBinCenter(lBin))-h_cosl->GetBinContent(lBin))/h_cosl->GetBinError(lBin),2);
    }
    latex->DrawLatexNDC(0.01,0.95,TString::Format("#chi^{2} = %f",chi2Val));
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
    for (int iPar = 0; iPar < nPar; iPar++){
        output.push_back(arrPar[iPar]);
        output.push_back(arrParErr[iPar]);
    }
    for (int i = 0; i < output.size(); i=i+2) {
        printf("%f +- %f\n",output[i],output[i+1]);
    }
    return output;
}//}}}

/*
std::vector<double> angular_bin(int iBin, const char outfile[] = "angular")
{//{{{

    // Fit efficiency
    RooRealVar 

    RooDataSet *

    // Fit to data
    RooRealVar CosThetaK("CosThetaK", "cos#theta_{K}", -1, 1);
    RooRealVar CosThetaL("CosThetaL", "cos#theta_{L}", -1, 1);
    RooRealVar Q2("Q2","q^{2}",0.5,20.);
    RooRealVar fl("fl", "F_{L}", 0.99, -0.5, 1.5);
    RooRealVar afb("afb", "A_{FB}", -1., -10., 10.);
    RooRealVar fs("fs","F_{S}",0.013);//B0ToKstarJpsi
    RooRealVar as("as","A_{S}",-0.098);//B0ToKstarJpsi
    fs.setConstant(kTRUE);
    as.setConstant(kTRUE);

    RooGenericPdf f_bkg("f_bkg", "(fs+2*as*CosThetaK)*(1-CosThetaL*CosThetaL)*3/8", RooArgSet(CosThetaK,CosThetaL,fs,as));
    RooGenericPdf f_sig("f_sig", "(1-fs)*9/16*(2*fl*CosThetaK*CosThetaK*(1-CosThetaL*CosThetaL)+0.5*(1-fl)*(1-CosThetaK*CosThetaK)*(1+CosThetaL*CosThetaL)+4*afb*(1-CosThetaK*CosThetaK)*CosThetaL/3)", RooArgSet(CosThetaK,CosThetaL,fl,afb,fs,as));
    RooRealVar nsig("nsig","nsig",100,1E2,1E8);
    RooRealVar nbkg("nbkg","nbkg",10,10,1E8);
    RooAddPdf f("f","f",RooArgList(f_sig,f_bkg),RooArgList(nsig,nbkg));

    RooGenericPdf f_eff("f_eff", "", RooArgSet());
    
    RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(CosThetaK,CosThetaL,Q2),q2range[iBin],0);
    RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE));

    RooPlot* framecosk = CosThetaK.frame(); 
    data->plotOn(framecosk); 
    f.plotOn(framecosk); 

    // Draw the frame on the canvas
    TCanvas* c = new TCanvas("c");
    framecosk->SetTitle("");
    //framecosk->SetMaximum(0);
    framecosk->SetMinimum(0);
    framecosk->Draw();

    TLatex *t1 = new TLatex();
    t1->SetNDC();
    t1->DrawLatex(.50,.86,TString::Format("%s",q2range[iBin]));
    t1->DrawLatex(.50,.80,TString::Format("F_{L}=%5.3f#pm%5.3f",fl.getVal(),fl.getError()));
    t1->DrawLatex(.50,.74,TString::Format("A_{FB}=%5.3f#pm%5.3f",afb.getVal(),afb.getError()));
    if ( f_fitresult->status() == 0){
        //t1->DrawLatex(.50,.68,TString::Format("Fit status: %s","GOOD"));
        t1->DrawLatex(.50,.68,TString::Format("Fit status: %d",f_fitresult->covQual()));
    }else{
        t1->DrawLatex(.50,.68,TString::Format("Fit status: %s","BAD"));
        //t1->DrawLatex(.50,.68,TString::Format("Fit status: %d",f_fitresult->covQual()));
    }
    c->Print(TString::Format("./plots/%s_bin%d.pdf",outfile,iBin));

    std::vector<double> output;
    output.push_back(fl.getVal());
    output.push_back(fl.getError());
    output.push_back(afb.getVal());
    output.push_back(afb.getError());
    return output;

}//}}}
*/
//_________________________________________________________________________________


int main(int argc, char** argv) {
    if (argc <= 2) {
        printf("Need at least 2 arguments.\n");
        printf("./fit [bmass, fl, angular_gen, acceptance, efficiency, angular] infile\n");
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
        //acceptance(0);
    }else if (func == "efficiency") {
        ch->Add(infile.Data());
        if (ch == NULL) gSystem->Exit(0);
        for (int iBin = 0; iBin < 8; iBin++) {
            efficiency(iBin);
        }
        //efficiency(0);
    }else if (func == "angular"){
        ch->Add(infile.Data());
        if (ch == NULL) gSystem->Exit(0);
        const char outfile[]="angular";
        printf("This function is under developing...\n");
        //angular_bin(8,outfile);
        //angular(outfile);
    }else{ 
        cerr << "No function available for: " << func.Data() << endl; 
    }
    printf("%lld entries processed.\n",ch->GetEntries());
    gSystem->Exit(0);

    return 0 ;
}
