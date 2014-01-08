// vim: sw=4 ts=4 fdm=marker et:

// -----------------------------------------------
//       Author: Xin Shi <Xin.Shi@cern.ch> 
//       Created:   [2013-08-15 Thu 14:54] 
// -----------------------------------------------
#include <sstream>
#include <TSystem.h>
#include <TH1.h>
#include <TH2F.h>
#include <TFile.h>
#include <TPad.h> 
#include <TCanvas.h> 
#include <TChain.h> 
#include <TPaveText.h>
#include <TLatex.h>
#include <TString.h>
#include <TGraphAsymmErrors.h>

#include <RooRealVar.h>
#include <RooGaussian.h>
#include <RooChebychev.h> 
#include <RooAddPdf.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooGenericPdf.h> 
#include <RooPolynomial.h>
#include <RooProdPdf.h>
#include <RooDataHist.h>

#include "tools.h" 

using namespace std; 
using namespace RooFit ;

void bmass(TString datatype, TString label, TString cut, TString outfile)
{//{{{
  bool test = false; 

  // Importing a  TTree into a RooDataSet with cuts 
  // --------------------------------------------------------------------------
  TChain* ch = add_chain(datatype, label, cut); 
  if (ch == NULL) gSystem->Exit(0);

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
  paveText->SetBorderSize(0);
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

}//}}}

std::vector<double> fl_bin(int iBin, const char outfile[] = "fl")
{//{{{
  // From fomula (2) in LHCb 2012 PRL108, 181806(2012)
  // integrated over theta_l and phi: 
  // 
  // 1/Gamma * d^2 Gamma/d cos(theta_K) dq^2 = 3/2 * F_L cos^2(theta_K)
  // + 3/4(1-F_L)(1-cos^2theta_K)
  // 
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

  bool test = false;

  TChain *treein = new TChain("tree");
  //treein->Add("BuToKstarMuMu_test.root");
  treein->Add("BuToKstar_merged.root");

  RooRealVar genCosThetaK("genCosThetaK", "cos#theta_{K}", -1, 1);
  RooRealVar Q2("Q2","q^{2}",0.5,20.);
  RooRealVar fl("fl", "F_{L}", 0.5, -0.5, 1.5);

  RooGenericPdf f("f", "1.5*fl*genCosThetaK*genCosThetaK+0.75*(1-fl)*(1-genCosThetaK*genCosThetaK)", RooArgSet(genCosThetaK,fl));
  RooDataSet* data;
  
  if (test){
      fl.setVal(0.5);
      data = f.generate(RooArgSet(genCosThetaK,Q2), 10000);
  }else{
      data = new RooDataSet("data","data",treein,RooArgSet(genCosThetaK,Q2),q2range[iBin],0);
  }
  
  f.fitTo(*data); 
  //f.fitTo(*data,Extended()); 

  RooPlot* framecosk = genCosThetaK.frame(); 
  data->plotOn(framecosk); 
  f.plotOn(framecosk); 

  // Draw the frame on the canvas
  TCanvas* c = new TCanvas("c"); 
  framecosk->SetTitle("");
  framecosk->Draw();

  TLatex *t1 = new TLatex();
  t1->SetNDC();
  t1->DrawLatex(.40,.85,TString::Format("%s",q2range[iBin]));
  t1->DrawLatex(.40,.79,TString::Format("F_{L}=%5.3f#pm%5.3f",fl.getVal(),fl.getError()));

  c->Print(TString::Format("./plots/%s_bin%d.pdf",outfile,iBin));

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
}//}}}

std::vector<double> angular_gen_bin(int iBin, const char outfile[] = "angular_gen")
{//{{{
    char q2range[10][32] = {"genQ2 < 2.00 && genQ2 > 1.00",
                            "genQ2 < 4.30 && genQ2 > 2.00",
                            "genQ2 < 8.68 && genQ2 > 4.30",
                            "genQ2 <10.09 && genQ2 > 8.68",
                            "genQ2 <12.86 && genQ2 >10.09",
                            "genQ2 <14.18 && genQ2 >12.86",
                            "genQ2 <16.00 && genQ2 >14.18",
                            "genQ2 <19.00 && genQ2 >16.00",
                            "genQ2 <19.00 && genQ2 > 1.00",
                            "genQ2 < 6.00 && genQ2 > 1.00"};

    TChain *treein = new TChain("tree");
    //treein->Add("BuToKstarMuMu_PtEtaFilter_7TeV_1E6.root");
    treein->Add("BuToKstarMuMu_NoFilter_staticB_1E6.root");
    //treein->Add("BuBarToKstarMuMu_NoFilter_1E6.root");

    RooRealVar genCosThetaK("genCosThetaK", "cos#theta_{K}", -1., 1.);
    RooRealVar genCosThetaL("genCosThetaL", "cos#theta_{L}", -1., 1.);
    RooRealVar genQ2("genQ2","q^{2}",0.5,20.);
    RooRealVar fl("fl", "F_{L}", 0.8, -0.2, 1.2);
    RooRealVar afb("afb", "A_{FB}", 0., -1., 1.);
    RooRealVar fs("fs","F_{S}",0.);//B0ToKstarJpsi
    RooRealVar as("as","A_{S}",0.);//B0ToKstarJpsi
    fs.setConstant(kTRUE);
    as.setConstant(kTRUE);

    RooRealVar nsig("nsig","nsig",1E4,1E2,1E8);
    RooRealVar nbkg("nbkg","nbkg",10,0.1,1E3);
    
    RooGenericPdf f_sig("f_sig", "9/16*((2/3*fs+4/3*as*genCosThetaK)*(1-genCosThetaL*genCosThetaL)+(1-fs)*(2*fl*genCosThetaK*genCosThetaK*(1-genCosThetaL*genCosThetaL)+1/2*(1-fl)*(1-genCosThetaK*genCosThetaK)*(1+genCosThetaL*genCosThetaL)+4/3*afb*(1-genCosThetaK*genCosThetaK)*genCosThetaL))", RooArgSet(genCosThetaK,genCosThetaL,fl,afb,fs,as));
    RooGenericPdf f_bkg("f_bkg", "1",RooArgSet());
    //nbkg.setConstant(kTRUE);
    RooAddPdf f("f","f",RooArgList(f_sig,f_bkg),RooArgList(nsig,nbkg));
    RooDataSet *data = new RooDataSet("data","data",treein,RooArgSet(genCosThetaK,genCosThetaL,genQ2),q2range[iBin],0);

    RooFitResult *f_fitresult = f.fitTo(*data,Extended(),Save(kTRUE),Minimizer("Minuit2"));

    // Draw the frame on the canvas
    TCanvas* c = new TCanvas("c");
    RooPlot* framecosk = genCosThetaK.frame(); 
    data->plotOn(framecosk,Binning(100)); 
    f.plotOn(framecosk); 
    //f.plotOn(framecosk,Components(f_sig),LineColor(4),LineWidth(2));
    f.plotOn(framecosk,Components(f_bkg),LineStyle(2),LineColor(8),LineWidth(2));

    framecosk->SetTitle("");
    //framecosk->SetMaximum(0);
    framecosk->SetMinimum(0);
    framecosk->Draw();

    TLatex *t1 = new TLatex();
    t1->SetNDC();
    double fixNDC = -0.5;
    if (iBin < 5) fixNDC = 0.;
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
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
    //framecosl->SetMaximum(0);
    framecosl->SetMinimum(0);
    framecosl->Draw();

    if (iBin < 5) fixNDC = 0.;
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    c->Update();
    c->Print(TString::Format("./plots/%s_cosl_gen_bin%d.pdf",outfile,iBin));

    // Make 2-D plot
    TH1 *h1 = data->createHistogram("genCosThetaL,genCosThetaK", 100, 100);
    h1->Draw("LEGO");
    c->Update();
    c->Print(TString::Format("./plots/%s_2d_gen_bin%d.pdf",outfile,iBin));

    std::vector<double> output;
    output.push_back(fl.getVal());
    output.push_back(fl.getError());
    output.push_back(afb.getVal());
    output.push_back(afb.getError());
    return output;

}//}}}

void angular_gen(const char outfile[] = "angular")
{//{{{

    TCanvas *c = new TCanvas();
    TH2F *frame = new TH2F("frame","",18,1,19,10,0,1.5);
    frame->SetStats(kFALSE);
    frame->SetXTitle("q^{2} [(GeV)^{2}]");
    frame->SetYTitle("F_{L}");
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
    g_fl->Draw("P*");
    c->Print(TString::Format("./plots/%s_fl.pdf",outfile));

    c->Clear();
    frame->SetYTitle("A_{FB}");
    frame->Draw();
    //TGraphAsymmErrors *g_afb = new TGraphAsymmErrors(8,x,yafb,xerr,xerr,yerrafb,yerrafb);
    //g_afb->Draw("ap");
    //c->Print(TString::Format("./plots/%s_afb.pdf",outfile));
}//}}}

std::vector<double> efficiency(int iBin) // acceptance times reconstruction efficiency
{//{{{
    double q2rangedn[10] = {1.00 , 2.00 , 4.30 , 8.68  , 10.09 , 12.86 , 14.18 , 16.00 ,  1.00 , 1.00};
    double q2rangeup[10] = {2.00 , 4.30 , 8.68 , 10.09 , 12.86 , 14.18 , 16.00 , 19.00 , 19.00 , 6.00};
    std::vector<double> output;

    double BMass = 0;
    double genQ2 = 0;
    double gCosThetaK = 0;
    double gCosThetaL = 0;

    TChain *treein = new TChain("tree");
    treein->Add("BuToKstarMuMu_1M.root");
    treein->SetBranchStatus("*",0);
    treein->SetBranchStatus("genQ2"         , 1);
    treein->SetBranchStatus("genCosTheta*"  , 1);
    treein->SetBranchStatus("BMass"         , 1);
    treein->SetBranchAddress("genQ2"        , &genQ2);
    treein->SetBranchAddress("genCosThetaK" , &gCosThetaK);
    treein->SetBranchAddress("genCosThetaL" , &gCosThetaL);
    treein->SetBranchAddress("BMass"        , &BMass);

    // Fill histograms
    TH2F h2_ngen("h2_ngen" ,"h2_ngen" ,10,-1.,1.,10,-1.,1.); 
    TH2F h2_nacc("h2_nacc" ,"h2_nacc" ,10,-1.,1.,10,-1.,1.); 
    TH2F h2_nreco("h2_nreco","h2_nreco",10,-1.,1.,10,-1.,1.);
    for (int entry = 0; entry < treein->GetEntries(); entry++) {
        treein->GetEntry(entry);
        if (genQ2 > q2rangeup[iBin] || genQ2 < q2rangedn[iBin]) continue;
        h2_ngen.Fill(gCosThetaK,gCosThetaL);
        
        if (BMass != 0) h2_nreco.Fill(gCosThetaK,gCosThetaL);
    }
    
    // Calculate efficiency
    TH2F h2_eff("h2_eff","h2_eff",10,-1.,1.,10,-1.,1.);
    for (int i = 1; i <= 10; i++) {
        for (int j = 1; j <= 10; j++) {
            if (h2_ngen.GetBinContent(i,j) == 0) {
                printf("ERROR: h2_gen_bin%d(%d,%d)=0\n",iBin,i,j);
                h2_ngen.SetBinContent(i,j,0.01);
            }
            h2_eff.SetBinContent(i,j,h2_nreco.GetBinContent(i,j)/h2_ngen.GetBinContent(i,j));
        }
    }

    // Draw efficiency
    TCanvas canvas("canvas");
    h2_eff.Draw("LEGO");
    canvas.Print(TString::Format("./plots/efficiency_bin%d.pdf",iBin));

    // Define fit function
    RooRealVar genCosThetaL("genCosThetaL","genCosThetaL",-1,1);
    RooRealVar genCosThetaK("genCosThetaK","genCosThetaK",-1,1);
    RooRealVar c2("c2","c2",1,-1E2,1E2);
    RooRealVar c1("c1","c1",1,-1E2,1E2);
    RooRealVar c0("c0","c0",1,-1E2,1E2);
    RooRealVar b2("b2","b2",1,-1E2,1E2);
    RooRealVar b1("b1","b1",1,-1E2,1E2);
    RooRealVar b0("b0","b0",1,-1E2,1E2);
    RooRealVar a2("a2","a2",1,-1E2,1E2);
    RooRealVar a1("a1","a1",1,-1E2,1E2);
    RooPolynomial *f_amp = new RooPolynomial("f_amp","f_amp",genCosThetaK,RooArgList(a2,a1));
    RooGenericPdf *f_exp = new RooGenericPdf("f_exp","f_exp","exp(-0.5*((genCosThetaL-c2*genCosThetaK**2-c1*genCosThetaK-c0)/(b2*genCosThetaK**2+b1*genCosThetaK+b0))**2)",RooArgList(genCosThetaL,genCosThetaK,c2,c1,c0,b2,b1,b0));
    RooProdPdf *f = new RooProdPdf("f","f",RooArgList(*f_amp,*f_exp));
    
    // Fit efficiency
    RooDataHist *data = new RooDataHist("data","data",RooArgList(genCosThetaK,genCosThetaL),&h2_eff);
    RooFitResult *f_fitresult = f->fitTo(*data,Save(kTRUE));
    output.push_back(c2.getVal());
    output.push_back(c2.getError());
    output.push_back(c1.getVal());
    output.push_back(c1.getError());
    output.push_back(c0.getVal());
    output.push_back(c0.getError());
    output.push_back(b2.getVal());
    output.push_back(b2.getError());
    output.push_back(b1.getVal());
    output.push_back(b1.getError());
    output.push_back(b0.getVal());
    output.push_back(b0.getError());
    output.push_back(a2.getVal());
    output.push_back(a2.getError());
    output.push_back(a1.getVal());
    output.push_back(a1.getError());
    
    // Draw FitResult
    RooPlot* framecosk = genCosThetaK.frame(); 
    data->plotOn(framecosk); 
    f->plotOn(framecosk); 
    framecosk->Draw();
    canvas.Update();
    canvas.Print(TString::Format("./plots/fitEfficiency_cosK_bin%d.pdf",iBin));
    
    RooPlot* framecosl = genCosThetaL.frame(); 
    data->plotOn(framecosl); 
    f->plotOn(framecosl); 
    framecosl->Draw();
    canvas.Update();
    canvas.Print(TString::Format("./plots/fitEfficiency_cosL_bin%d.pdf",iBin));

    // Clear
    treein = 0;
    delete treein;

    return output;
}//}}}

/*
std::vector<double> angular_bin(int iBin, const char outfile[] = "angular")
{//{{{
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

    TChain *treein = new TChain("tree");
    //treein->Add("BuToKstarMuMu_test.root");
    treein->Add("BuToKstar_merged.root");

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
    
    RooDataSet *data = new RooDataSet("data","data",treein,RooArgSet(CosThetaK,CosThetaL,Q2),q2range[iBin],0);
    RooFitResult *f_fitresult = f.fitTo(*data,Extended(),Save(kTRUE));

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


int main(int argc, char** argv) {
    if (argc == 1) {
        printf("argv[0]=%s\n",argv[0]);
        return 0;
    }

    TString func    = argv[1];
    TString infile  = argv[2];

    if (func == "help" || func == "-h"){
        printf("./fit [bmass, angular_gen,angular,fl] infile\n");
    }else if (func == "bmass") {
        TString datatype = argv[3]; 
        TString label    = argv[4]; 
        TString cut      = argv[5]; 
        TString outfile  = argv[6]; 

        bmass(datatype, label, cut, "bmass"); 
    }else if (func == "fl"){
        fl("fl");
    }else if (func == "angular_gen"){
        const char outfile[]="angular_gen";
        //angular_bin_gen(8,outfile);
        angular_gen(outfile);
    }else if (func == "angular"){
        const char outfile[]="angular";
        //angular_bin(8,outfile);
        //angular(outfile);
    }else if (func == "efficiency") {
        for (int iBin = 0; iBin < 8; iBin++) {
            efficiency(iBin);
        }
    }else{ 
        cerr << "No function available for: " << func.Data() << endl; 
    }
    gSystem->Exit(0);

    return 0 ;
}
