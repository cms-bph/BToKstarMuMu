#define SingleBdToKstarMuMuSelector_xCheck2011_cxx
// The class definition in SingleBdToKstarMuMuSelector_xCheck2011.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("SingleBdToKstarMuMuSelector_xCheck2011.C")
// Root > T->Process("SingleBdToKstarMuMuSelector_xCheck2011.C","some options")
// Root > T->Process("SingleBdToKstarMuMuSelector_xCheck2011.C+")
//

#include <iostream>
#include <sstream>
#include <map>
#include "SingleBdToKstarMuMuSelector_xCheck2011.h"
#include <TH2.h>
#include <TStyle.h>
#include <TProof.h>
#include <TLorentzVector.h>

// Global Constants
const double KSTAR_MASS = 0.89166; // GeV
const double KSTAR_WIDTH = 0.0508; // GeV 
const double MUON_MASS = 0.10565837;
const double KAON_MASS = 0.493677;
const double PION_MASS = 0.13957018;
const double KSHORT_MASS = 0.497614;

// user defined variables
TDatime t_begin_ , t_now_ ;
int n_processed_, n_selected_; 

TTree *tree_; 

// Branch variables for new tree
int    Nb          = 0;
double Mumumass    = 0;
double Mumumasserr = 0;
double Kstarmass   = 0;
double Trkpt       = 0;
double Trkdcasigbs = 0;

double Kshortpt    = 0;
double Pimpt       = 0;
double Pippt       = 0;

double Bmass       = 0;
double Bpt         = 0;
double Bphi        = 0;
int    Bchg        = 0;
double Bvtxcl      = 0;
double Blxysig     = 0;
double Bcosalphabs = 0;
double Bctau       = 0;

double Q2 = 0;
double CosThetaL   = 999;
double CosThetaK   = 999;

// Branches NOT included in 2012 analysis
int TrkChg = 999;
double EvtWeight = 1;

// Branches for Generator level information
int     genBChg      = 999;
double  genBPt       = 0;
double  genBEta      = 0;
double  genBPhi      = 0;
double  genBVtxX     = 0;
double  genBVtxY     = 0;
double  genBVtxZ     = 0;
double  genMupPt     = 0;
double  genMupEta    = 0;
double  genMupPhi    = 0;
double  genMumPt     = 0;
double  genMumEta    = 0;
double  genMumPhi    = 0;
double  genKstPt     = 0;
double  genKstEta    = 0;
double  genKstPhi    = 0;
int     genTkChg     = 999;
double  genTkPt      = 0;
double  genTkEta     = 0;
double  genTkPhi     = 0;
double  genKPt       = 0;//Could be K+-/Kshort
double  genKEta      = 0;
double  genKPhi      = 0;
double  genKVtxX     = 0;
double  genKVtxY     = 0;
double  genKVtxZ     = 0;
double  genPipPt     = 0;//pion pair decayed from Kshort
double  genPipEta    = 0;
double  genPipPhi    = 0;
double  genPimPt     = 0;
double  genPimEta    = 0;
double  genPimPhi    = 0;
double  genQ2        = 0;
double  genCosThetaL = 999;
double  genCosThetaK = 999;

void ClearEvent()
{//{{{
    Nb          = 0;
    Mumumass    = 0;
    Mumumasserr = 0;
    Kstarmass   = 0;
    Trkpt       = 0;
    Trkdcasigbs = 0;
    
    Kshortpt    = 0;
    Pimpt       = 0;
    Pippt       = 0;
    
    Bmass       = 0;
    Bpt         = 0;
    Bphi        = 0;
    Bchg        = 0;
    Bvtxcl      = 0;
    Blxysig     = 0;
    Bcosalphabs = 0;
    Bctau       = 0;

    Q2          = 0;
    CosThetaL   = 999;
    CosThetaK   = 999;

        // 2011 special
    TrkChg      = 999;
    EvtWeight   = 1;

        // mc
    genBChg = 999;
    genBPt = 0;
    genBEta= 0;
    genBPhi= 0;
    genBVtxX = 0;
    genBVtxY = 0;
    genBVtxZ = 0;
    genMupPt = 0;
    genMupEta= 0;
    genMupPhi= 0;
    genMumPt = 0;
    genMumEta= 0;
    genMumPhi= 0;
    genKstPt = 0;
    genKstEta= 0;
    genKstPhi= 0;
    genTkChg = 999;
    genTkPt = 0;
    genTkEta= 0;
    genTkPhi= 0;
    genKPt = 0;//Could be K+-/Kshort
    genKEta= 0;
    genKPhi= 0;
    genKVtxX = 0;
    genKVtxY = 0;
    genKVtxZ = 0;
    genPipPt = 0;//pion pair decayed from Kshort
    genPipEta= 0;
    genPipPhi= 0;
    genPimPt = 0;
    genPimEta= 0;
    genPimPhi= 0;
    genQ2 = 0;
    genCosThetaL = 999;
    genCosThetaK = 999;
}//}}}

void str_replace(std::string& str, const std::string& oldStr, const std::string& newStr)
{//{{{
  size_t pos = 0;
  while((pos = str.find(oldStr, pos)) != std::string::npos)
  {
    str.replace(pos, oldStr.length(), newStr);
    pos += newStr.length();
  }
}//}}}

string get_option_value(string option, string name)
{//{{{
  vector<string> args;
  istringstream f(option);
  string s;    
  while (getline(f, s, ';')) {
    args.push_back(s);
  }
  
  string value; 
  for(vector<string>::iterator it = args.begin(); it != args.end(); ++it) {
    value = *it; 
    unsigned found = value.find(name);
    if (found == 0) {
      str_replace(value, name+"=", ""); 
      break; 
    }
  }
  return value; 
}//}}}

void SingleBdToKstarMuMuSelector_xCheck2011::Begin(TTree * /*tree*/)
{//{{{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).
  t_begin_.Set(); 
  printf("\n ---------- Begin Job ---------- \n");
  t_begin_.Print();

  n_processed_ = 0;
  n_selected_ = 0;
}//}}}

void SingleBdToKstarMuMuSelector_xCheck2011::SlaveBegin(TTree * /*tree*/)
{//{{{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

    string option = GetOption();
    tree_ = new TTree("tree", "tree"); 
    
    tree_->Branch("Mumumass"    , &Mumumass    , "Mumumass/D");
    tree_->Branch("Mumumasserr" , &Mumumasserr , "Mumumasserr/D");
    tree_->Branch("Kstarmass"   , &Kstarmass   , "Kstarmass/D");
    tree_->Branch("Trkpt"       , &Trkpt       , "Trkpt/D");
    tree_->Branch("Trkdcasigbs" , &Trkdcasigbs , "Trkdcasigbs/D");

    tree_->Branch("Kshortpt"    , &Kshortpt    , "Kshortpt/D");
    tree_->Branch("Pimpt"       , &Pimpt       , "Pimpt/D");
    tree_->Branch("Pippt"       , &Pippt       , "Pippt/D");

    tree_->Branch("Bmass"       , &Bmass       , "Bmass/D");
    tree_->Branch("Bpt"         , &Bpt         , "Bpt/D");
    tree_->Branch("Bphi"         , &Bphi         , "Bphi/D");
    tree_->Branch("Bchg"        , &Bchg        , "Bchg/I");

    tree_->Branch("Bvtxcl"      , &Bvtxcl      , "Bvtxcl/D");
    tree_->Branch("Blxysig"     , &Blxysig     , "Blxysig/D");
    tree_->Branch("Bcosalphabs" , &Bcosalphabs , "Bcosalphabs/D");
    tree_->Branch("Bctau"       , &Bctau       , "Bctau/D");
    
    tree_->Branch("Q2"          , &Q2          , "Q2/D");
    tree_->Branch("CosThetaL"   , &CosThetaL   , "CosThetaL/D");
    tree_->Branch("CosThetaK"   , &CosThetaK   , "CosThetaK/D");

    tree_->Branch("TrkChg"      , &TrkChg       , "TrkChg/I");
    tree_->Branch("EvtWeight"   , &EvtWeight    , "EvtWeight/D");

    string datatype = get_option_value(option, "datatype"); 
    std::map<string,int> maptype;
    maptype.insert(std::pair<string,int>("data",1));
    maptype.insert(std::pair<string,int>("mc.lite",2));
    maptype.insert(std::pair<string,int>("mc",999));
    switch (maptype[datatype]) {
        case 1:
            break;
        case 2:
            tree_->Branch("genBChg"      , &genBChg      , "genBChg/I");
            tree_->Branch("genBPhi"      , &genBPhi      , "genBPhi/D");
            tree_->Branch("genMupPt"     , &genMupPt     , "genMupPt/D");
            tree_->Branch("genMupEta"    , &genMupEta    , "genMupEta/D");
            tree_->Branch("genMumPt"     , &genMumPt     , "genMumPt/D");
            tree_->Branch("genMumEta"    , &genMumEta    , "genMumEta/D");
            tree_->Branch("genQ2"        , &genQ2        , "genQ2/D");
            tree_->Branch("genCosThetaL" , &genCosThetaL , "genCosThetaL/D");
            tree_->Branch("genCosThetaK" , &genCosThetaK , "genCosThetaK/D");
            break;
        case 999:
            tree_->Branch("genBChg"      , &genBChg      , "genBChg/I");
            tree_->Branch("genBPt"       , &genBPt       , "genBPt/D");
            tree_->Branch("genBEta"      , &genBEta      , "genBEta/D");
            tree_->Branch("genBPhi"      , &genBPhi      , "genBPhi/D");
            tree_->Branch("genBVtxX"     , &genBVtxX     , "genBVtxX/D");
            tree_->Branch("genBVtxY"     , &genBVtxY     , "genBVtxY/D");
            tree_->Branch("genBVtxZ"     , &genBVtxZ     , "genBVtxZ/D");
            tree_->Branch("genMupPt"     , &genMupPt     , "genMupPt/D");
            tree_->Branch("genMupEta"    , &genMupEta    , "genMupEta/D");
            tree_->Branch("genMupPhi"    , &genMupPhi    , "genMupPhi/D");
            tree_->Branch("genMumPt"     , &genMumPt     , "genMumPt/D");
            tree_->Branch("genMumEta"    , &genMumEta    , "genMumEta/D");
            tree_->Branch("genMumPhi"    , &genMumPhi    , "genMumPhi/D");
            tree_->Branch("genKstPt"     , &genKstPt     , "genKstPt/D");
            tree_->Branch("genKstEta"    , &genKstEta    , "genKstEta/D");
            tree_->Branch("genKstPhi"    , &genKstPhi    , "genKstPhi/D");
            tree_->Branch("genTkChg"     , &genTkChg     , "genTkChg/I");
            tree_->Branch("genTkPt"      , &genTkPt      , "genTkPt/D");
            tree_->Branch("genTkEta"     , &genTkEta     , "genTkEta/D");
            tree_->Branch("genTkPhi"     , &genTkPhi     , "genTkPhi/D");
            tree_->Branch("genKPt"       , &genKPt       , "genKPt/D");
            tree_->Branch("genKEta"      , &genKEta      , "genKEta/D");
            tree_->Branch("genKPhi"      , &genKPhi      , "genKPhi/D");
            tree_->Branch("genKVtxX"     , &genKVtxX     , "genKVtxX/D");
            tree_->Branch("genKVtxY"     , &genKVtxY     , "genKVtxY/D");
            tree_->Branch("genKVtxZ"     , &genKVtxZ     , "genKVtxZ/D");
            tree_->Branch("genPipPt"     , &genPipPt     , "genPipPt/D");
            tree_->Branch("genPipEta"    , &genPipEta    , "genPipEta/D");
            tree_->Branch("genPipPhi"    , &genPipPhi    , "genPipPhi/D");
            tree_->Branch("genPimPt"     , &genPimPt     , "genPimPt/D");
            tree_->Branch("genPimEta"    , &genPimEta    , "genPimEta/D");
            tree_->Branch("genPimPhi"    , &genPimPhi    , "genPimPhi/D");
            tree_->Branch("genQ2"        , &genQ2        , "genQ2/D");
            tree_->Branch("genCosThetaL" , &genCosThetaL , "genCosThetaL/D");
            tree_->Branch("genCosThetaK" , &genCosThetaK , "genCosThetaK/D");
            break;
        default:
            printf("No compatible datatype found. Please check use following types...\n\t\t[");
            for (std::map<string,int>::iterator iType = maptype.begin(); iType != maptype.end(); iType++){
                if (iType->second != 0) printf("%s,",iType->first.c_str());
            }
            printf("]\n");
            break;
        }

            fOutput->AddAll(gDirectory->GetList()); 

    }//}}}

Bool_t SingleBdToKstarMuMuSelector_xCheck2011::Process(Long64_t entry)
{//{{{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either SingleBdToKstarMuMuSelector_xCheck2011::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

    ClearEvent();

    string option = GetOption();
    string datatype = get_option_value(option, "datatype"); 
    string cut = get_option_value(option, "cut"); 

    GetEntry(entry); 
    n_processed_ += 1; 
    Nb = nB; 
    //printf("Enter entry %lld with %d B candidate(s).\n",entry,Nb);

    if (datatype != "data") SaveGen();

    int i = SelectB(cut); 
    if ( i != -1 ) {
        n_selected_ += 1; 
        SaveEvent(i);     
    }

    tree_->Fill();	   
    return kTRUE;

}//}}}

void SingleBdToKstarMuMuSelector_xCheck2011::SlaveTerminate()
{//{{{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}//}}}

void SingleBdToKstarMuMuSelector_xCheck2011::Terminate()
{//{{{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
   //
    string option = GetOption();
    TString outfile = get_option_value(option, "ofile"); 

    TFile file(outfile.Data(), "recreate"); 
    fOutput->Write();

    t_now_.Set(); 
    printf(" \n ---------- End Job ---------- \n" ) ;
    t_now_.Print();  
    printf(" processed: %i \n selected: %i \n \
            duration: %i sec \n rate: %g evts/sec\n",
            n_processed_, n_selected_, 
            t_now_.Convert() - t_begin_.Convert(), 
            float(n_processed_)/(t_now_.Convert()-t_begin_.Convert()) );

}//}}}

int SingleBdToKstarMuMuSelector_xCheck2011::SelectB(string cut)
{//{{{
    // To be fixed

    int best_idx = -1; 
    double best_bvtxcl = 0.0; 

    if (cut == "singleCand") {
        best_idx = 0;
    }else if (cut == "singleCandSig"){
        if ( (mumuMass->at(0) > 3.096916+3*mumuMassE->at(0) || mumuMass->at(0) < 3.096916-5*mumuMassE->at(0)) && fabs(mumuMass->at(0)-3.686109) > 3*mumuMassE->at(0) ) best_idx = 0;
    }else if (cut == "singleCandBkg"){
        if ( (mumuMass->at(0) < 3.096916+3*mumuMassE->at(0) && mumuMass->at(0) > 3.096916-5*mumuMassE->at(0)) || fabs(mumuMass->at(0)-3.686109) < 3*mumuMassE->at(0) ) best_idx = 0;
    }else if (cut == "cut0"){
        for (int i = 0; i < nB; i++) {

            if ( ! HasGoodDimuon(i) ) continue; 

            if (bVtxCL->at(i) > best_bvtxcl) {
                best_bvtxcl = bVtxCL->at(i); 
                best_idx = i; 
            }
        }
    }else{
        printf("ERROR: Unknown cut %s\n",cut.c_str());
    }
    return best_idx;
}//}}}

bool SingleBdToKstarMuMuSelector_xCheck2011::HasGoodDimuon(int i)
{//{{{
    // To be fixed
    return true; 
}//}}}

void SingleBdToKstarMuMuSelector_xCheck2011::SaveEvent(int i)
{//{{{
    TLorentzVector B_4vec, Kst_4vec, Mup_4vec, Mum_4vec, Tk_4vec, K_4vec, Pip_4vec, Pim_4vec, buff1, buff2, buff3;
    B_4vec.SetPtEtaPhiM(B0pT,B0Eta,B0Phi,B0MassArb);
    Mup_4vec.SetXYZM(mupPx->at(i),mupPy->at(i),mupPz->at(i),MUON_MASS);
    Mum_4vec.SetXYZM(mumPx->at(i),mumPy->at(i),mumPz->at(i),MUON_MASS);
    Kst_4vec.SetXYZM(kstPx->at(i),kstPy->at(i),kstPz->at(i),kstMass->at(i));
    // In 2011 analysis, both mass hypothesis are computed.
    if (B0notB0bar){
        K_4vec.SetXYZM(kstTrkpPx->at(i),kstTrkpPy->at(i),kstTrkpPz->at(i),KAON_MASS);//B0->K^+ + Pi^-
        Tk_4vec.SetXYZM(kstTrkmPx->at(i),kstTrkmPy->at(i),kstTrkmPz->at(i),PION_MASS);
        TrkChg = -1;
        Trkdcasigbs = fabs(kstTrkmDCABS->at(i)/kstTrkmDCABSE->at(i));
    }else{
        Tk_4vec.SetXYZM(kstTrkpPx->at(i),kstTrkpPy->at(i),kstTrkpPz->at(i),PION_MASS);
        K_4vec.SetXYZM(kstTrkmPx->at(i),kstTrkmPy->at(i),kstTrkmPz->at(i),KAON_MASS);
        TrkChg = 1;
        Trkdcasigbs = fabs(kstTrkpDCABS->at(i)/kstTrkpDCABSE->at(i));
    }
    Pip_4vec.SetXYZM(0,0,0,0);
    Pim_4vec.SetXYZM(0,0,0,0);

    Bmass = B0MassArb; 
    Bchg = 0; 
    Bvtxcl = bVtxCL->at(i); 
    Blxysig = (bLBS->at(i)/bLBSE->at(i));
    Bcosalphabs = bCosAlphaBS->at(i); 
    Bctau = bctauPVBS->at(i); 

    //Kshortpt = sqrt( (kspx->at(i))*(kspx->at(i)) + (kspy->at(i))*(kspy->at(i)) );
    //Pimpt = sqrt( (pimpx->at(i))*(pimpx->at(i)) + (pimpy->at(i))*(pimpy->at(i)) );
    //Pippt = sqrt( (pippx->at(i))*(pippx->at(i)) + (pippy->at(i))*(pippy->at(i)) );

    Bpt = B_4vec.Pt(); 
    Bphi = B_4vec.Phi();

    Mumumass = mumuMass->at(i); 
    Mumumasserr = mumuMassE->at(i); 
    Kstarmass = kstMass->at(i); 
    Trkpt = Tk_4vec.Pt(); 

    Q2 = pow(mumuMass->at(i),2);
    CosThetaL = CosThetaMuArb;
    CosThetaK = CosThetaKArb;
    EvtWeight = evWeight;
}//}}}

void SingleBdToKstarMuMuSelector_xCheck2011::SaveGen()
{//{{{
    TLorentzVector genB_4vec, genKst_4vec, genMup_4vec, genMum_4vec, genTk_4vec, genK_4vec, genPip_4vec, genPim_4vec, buff1, buff2, buff3;
    genB_4vec.SetXYZM(genB0Px,genB0Py,genB0Pz,genB0Mass);
    genMup_4vec.SetXYZM(genMupPx,genMupPy,genMupPz,MUON_MASS);
    genMum_4vec.SetXYZM(genMumPx,genMumPy,genMumPz,MUON_MASS);
    genKst_4vec.SetXYZM(genKstPx,genKstPy,genKstPz,KSTAR_MASS);
    // In 2011 analysis, both mass hypothesis are computed.
    genTk_4vec.SetXYZM(genKstTrkmPx,genKstTrkmPy,genKstTrkmPz,PION_MASS);
    genK_4vec.SetXYZM(genKstTrkpPx,genKstTrkpPy,genKstTrkpPz,KAON_MASS);//K*0 -> K+ Pi-
    buff1 = genTk_4vec + genK_4vec;
    genK_4vec.SetXYZM(genKstTrkmPx,genKstTrkmPy,genKstTrkmPz,KAON_MASS);
    genTk_4vec.SetXYZM(genKstTrkpPx,genKstTrkpPy,genKstTrkpPz,PION_MASS);
    buff2 = genTk_4vec + genK_4vec;
    if ( fabs(buff1.Mag()-KSTAR_MASS) < fabs(buff2.Mag()-KSTAR_MASS) ){
        genTk_4vec.SetXYZM(genKstTrkmPx,genKstTrkmPy,genKstTrkmPz,PION_MASS);
        genK_4vec.SetXYZM(genKstTrkpPx,genKstTrkpPy,genKstTrkpPz,KAON_MASS);
        genTkChg = -1;
    }else{
        genK_4vec.SetXYZM(genKstTrkmPx,genKstTrkmPy,genKstTrkmPz,KAON_MASS);
        genTk_4vec.SetXYZM(genKstTrkpPx,genKstTrkpPy,genKstTrkpPz,PION_MASS);
        genTkChg = 1;
    }
    genPip_4vec.SetXYZM(0,0,0,0);
    genPim_4vec.SetXYZM(0,0,0,0);

    genBChg      = 0;
    genBPt       = genB_4vec.Pt();
    genBEta      = genB_4vec.Eta();
    genBPhi      = genB_4vec.Phi();
    genBVtxX     = genPriVtxX;
    genBVtxY     = genPriVtxY;
    genBVtxZ     = genPriVtxZ;
    genMupPt     = genMup_4vec.Pt();
    genMupEta    = genMup_4vec.Eta();
    genMupPhi    = genMup_4vec.Phi();
    genMumPt     = genMum_4vec.Pt();
    genMumEta    = genMum_4vec.Eta();
    genMumPhi    = genMum_4vec.Phi();
    genKstPt     = genKst_4vec.Pt();
    genKstEta    = genKst_4vec.Eta();
    genKstPhi    = genKst_4vec.Phi();
    genTkPt      = genTk_4vec.Pt();
    genTkEta     = genTk_4vec.Eta();
    genTkPhi     = genTk_4vec.Phi();
    genKPt       = genK_4vec.Pt();
    genKEta      = genK_4vec.Eta();
    genKPhi      = genK_4vec.Phi();
    genKVtxX     = 0;
    genKVtxY     = 0;
    genKVtxZ     = 0;
    genPipPt     = 0;
    genPipEta    = 0;
    genPipPhi    = 0;
    genPimPt     = 0;
    genPimEta    = 0;
    genPimPhi    = 0;
    genQ2        = (genMup_4vec+genMum_4vec).Mag2();

    buff1 = genB_4vec;
    buff2 = genMup_4vec+genMum_4vec;
    buff1.Boost(-buff2.X()/buff2.T(),-buff2.Y()/buff2.T(),-buff2.Z()/buff2.T());
    if ( genTkChg < 0){
        buff3 = genMum_4vec;//Take mu- to avoid extra minus sign.
    }else{
        buff3 = genMup_4vec;
    }
    buff3.Boost(-buff2.X()/buff2.T(),-buff2.Y()/buff2.T(),-buff2.Z()/buff2.T());
    genCosThetaL = buff1.Vect().Dot(buff3.Vect())/buff1.Vect().Mag()/buff3.Vect().Mag();

    buff1 = genB_4vec;
    buff2 = genKst_4vec;
    buff1.Boost(-buff2.X()/buff2.T(),-buff2.Y()/buff2.T(),-buff2.Z()/buff2.T());
    buff3 = genTk_4vec;//Take pion to avoid extra minus.
    buff3.Boost(-buff2.X()/buff2.T(),-buff2.Y()/buff2.T(),-buff2.Z()/buff2.T());
    genCosThetaK = buff1.Vect().Dot(buff3.Vect())/buff1.Vect().Mag()/buff3.Vect().Mag();
}//}}}

#ifndef __CINT__ 
#include <algorithm>

char* get_option(char ** begin, char ** end, const std::string & option)
{//{{{
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)  return *itr;
    return 0;
}//}}}

bool option_exists(char** begin, char** end, const std::string& option)
{//{{{
    return std::find(begin, end, option) != end;
}//}}}

void print_usage()
{//{{{
    cerr << "Usage: ./sel2011 datatype cut infile outfile [-n] [-s] [-j] [-h]\n"
        << "   datatype: data, mc, mc.lite.\n"
        << "   cut     : singleCand, singleCandSig, singleCandbkg, cut0.\n"
        << "Options: \n" 
        << "  -h \t\tPrint this info.\n"
        << "  -n \t\tNumber of entries.\n" 
        << "  -s \t\tStarting run number.\n"
        << "  -j \t\tNumber of workers.\n" 
        << endl; 
}//}}}

int main(int argc, char** argv) {
    if ( (argc < 3) or option_exists(argv, argv+argc, "-h") ){
        print_usage() ;  
        return -1; 
    }

    TString datatype = argv[1]; 
    TString cut = argv[2]; 
    TString infile = argv[3]; 
    TString outfile = argv[4]; 

    Printf("datatype: '%s'", datatype.Data());
    Printf("cut: '%s'", cut.Data());
    Printf("input file: '%s'", infile.Data());
    Printf("output file: '%s'", outfile.Data());

    TChain *ch = new TChain("B0SingleCand/B0KstMuMuSingleCandNTuple"); 
    ch->Add(infile.Data()); 

    char * j = get_option(argv, argv+argc, "-j");
    if (j) {
        TProof::Open(Form("workers=%s", j));
        ch->SetProof(); 
    }

    Long64_t nentries = 1000000000; 
    char * n = get_option(argv, argv+argc, "-n");  
    if (n) nentries = atoi(n);
    
    
    int     iStart = 0;
    char *s = get_option(argv, argv+argc, "-s");
    if (s) {
        iStart = atoi(s);
        if (iStart > ch->GetEntries()){
            printf("ERROR: Number of entries is %lld.\n",ch->GetEntries());
            return -1;
        }
    }

    TString option; 
    option.Form("datatype=%s;cut=%s;ofile=%s_s%d.root", datatype.Data(), cut.Data(), outfile.Data(), iStart); 
    
    // It's not allowed to run with fat trees!
    if (datatype == "mc" && (!(s) || !(n))){
        printf("WARNING: You must specify #entries(-n) and start run(-s) for datatype '%s'.\n",datatype.Data());
        return -1;
    }

    ch->Process("SingleBdToKstarMuMuSelector_xCheck2011.cc+", option, nentries, iStart); 

    gSystem->Exit(0);

    return 0 ;
}

#endif
