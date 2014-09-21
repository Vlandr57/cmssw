// c++ classes
#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>

//ROOT
#include "TCanvas.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "TH2F.h"
#include "TString.h"
#include "TROOT.h"
#include "TStyle.h"

#include "TH1F.h"
#include "TF1.h"
#include "TFile.h"
#include "TPaveText.h"
#include "TTree.h"

#include "TMath.h"
#include <cmath>

Double_t twoGammaProfile(Double_t *x,Double_t *par) {
  double twoGamma=0;

  twoGamma = par[5] * par[0] * pow(x[0]*par[2],(par[1]-1.0)) 
           * exp(-x[0]*par[2]) / (ROOT::Math::tgamma(par[1])) 
           + par[5] * (1.-par[0]) * pow(x[0]*par[4],(par[3]-1.0))
           * exp(-x[0]*par[4]) / (ROOT::Math::tgamma(par[3]));


/*
  twoGamma = par[0] * pow(x[0]*par[2],(par[1]-1.0))
           * exp(-x[0]*par[2]) / (ROOT::Math::tgamma(par[1]))
           + (1.-par[0]) * pow(x[0]*par[4],(par[3]-1.0))
           * exp(-x[0]*par[4]) / (ROOT::Math::tgamma(par[3]));
*/

  return twoGamma;
}

void PionFit() {

   TFile *input(0);
   TString fname ="example_pi50Gev_new.root";
   input = TFile::Open( fname ); 

   double energyGen=0;
   double phiGen=0;
   double etaGen=0;
   double startingPoint=0;
   double interactPoint=0;
   double ePi0first=0;
   double ePi0tot=0;

   std::vector<double> *myx=0;
   std::vector<double> *myy=0;
   std::vector<double> *myz=0;
   std::vector<double> *myenergy=0;
   std::vector<double> *mytype=0;

//------------------------------------------------------
// mytype (1.0-8.0) determines where the deposited energy 
// is located in Ecal or Hcal (active or dead material):
//
// 1.0 - energy in HB scintillator  (barrel Hcal)
// 2.0 - energy in HB dead material (barrel Hcal)
//
// 3.0 - energy in HE scintillator  (endcap Hcal)
// 4.0 - energy in HE dead material (endcap Hcal)
//
// 5.0 - energy in EE active material (endcap Ecal)
// 6.0 - energy in EE dead material   (endcap Ecal)
//
// 7.0 - energy in EB active material (barrel Ecal) 
// 8.0 - energy in EB dead material   (barrel Ecal)
//----------------------------------------------------

   TTree* t1 = (TTree*)input->Get("demo/Event");
   TTree* t2 = (TTree*)input->Get("demo/Hits");

   t1->SetBranchAddress("energyGen",&energyGen);
   t1->SetBranchAddress("phiGen",&phiGen);
   t1->SetBranchAddress("etaGen",&etaGen);
   t1->SetBranchAddress("startingPoint",&startingPoint);
   t1->SetBranchAddress("interactPoint",&interactPoint);
   t1->SetBranchAddress("ePi0first",&ePi0first);
   t1->SetBranchAddress("ePi0tot",&ePi0tot);

   t2->SetBranchAddress("myx",&myx);
   t2->SetBranchAddress("myy",&myy);
   t2->SetBranchAddress("myz",&myz);
   t2->SetBranchAddress("myenergy",&myenergy);
   t2->SetBranchAddress("mytype",&mytype);

   int n = t1->GetEntries();
   cout << " Number of events =  " << n << endl;

   const int nHist(n);
   TH1D *h1;
   TH1D *hev[nHist];
//   TH1D *hev[100];
   char name[7];
   double alpha[100];
   double err_alpha[100];

   int    nBins = 24;
//   int    nBins = 36;
//   int    nBins = 48;
   double maxLength = 120.0;
   double z_depth;
   double e_dep;

   TH1D *h1 = new TH1D("h1","alpha " , 100, 0.0, 10.0);
   TH1D *h2 = new TH1D("h2","beta "  , 100, 0.0, 0.5);
   TH1D *h3 = new TH1D("h3","lambda" , 100, 0.0, 140.0);
   TH1D *h4 = new TH1D("h4","theta"  , 100, 0.0, 5.0);
   TH1D *h5 = new TH1D("h5","fsig"   , 100, 0.0, 1.0);

   TH1D *h6 = new TH1D("h6","alfa vs n" ,   100, 0.0, 100.0);
   TH1D *h7 = new TH1D("h7","beta vs n" ,   100, 0.0, 100.0);
   TH1D *h8 = new TH1D("h8","lambda vs n" , 100, 0.0, 100.0);
   TH1D *h9 = new TH1D("h9","theta vs n" ,  100, 0.0, 100.0);

   TH2D *h10 = new TH2D("h10","alfa vs beta" ,  100,0.0,10.0, 100,0.0,0.5);
   TH2D *h11 = new TH2D("h11","lambda vs theta",100,0.0,140.0,100,0.0,5.0);

//   for (int iev=0; iev<n; ++iev) {
   for (int iev=10; iev<20; ++iev) {

     t1->GetEntry(iev);
     t2->GetEntry(iev);

     snprintf(name,7,"#hev%d",iev);
     hev[iev] = new TH1D(name,"Shower depth ", nBins, 0.0, maxLength);
     hev[iev]->Sumw2();

     e_dep = 0;
     for (int ih=0; ih<myenergy->size(); ++ih) {
        z_depth = ( (*myz)[ih] - interactPoint)/10.; // in [cm] 
        if( ((*mytype)[ih]>=3 || (*mytype)[ih]<=4) && z_depth<=120.0 ) 
            e_dep +=(*myenergy)[ih];
     }

     for (int ih=0; ih<myenergy->size(); ++ih) {
        z_depth = ( (*myz)[ih] - interactPoint)/10.; // in [cm]
        if( e_dep>0.0 && ((*mytype)[ih]>=3 || (*mytype)[ih]<=4) ) 
        hev[iev]->Fill(z_depth,(*myenergy)[ih]/((maxLength/nBins)*e_dep));
//        hev[iev]->Fill(z_depth,(*myenergy)[ih]/(e_dep));
     }

// find maximum bin content
//----------------------------
   double maxCont = 0.0;
   for (Int_t bin=1; bin<hev[iev]->GetSize()-1; ++bin)
   {
      if( hev[iev]->GetBinContent(bin) > maxCont )
                    maxCont = hev[iev]->GetBinContent(bin); 
   }

//    cout << " Integral = " << hev[iev]->Integral() << endl;
//    cout << " maxCont = " << maxCont << endl;

    if(hev[iev]->Integral() < 0.05) continue;

//  TF1 *func = new TF1("twoGammaProfile", twoGammaProfile, 0, maxLength,5);
  TF1 *func = new TF1("twoGammaProfile", twoGammaProfile, 0, maxLength,6);

//  func->SetParameters(0.5, 3.0, 0.1, 20.0, 0.5);
  func->SetParameters(0.5, 3.0, 0.2, 20.0, 0.5, 0.5);
//  func->SetParNames ("c","alfa1","beta1","alfa2","beta2");
  func->SetParNames ("c","alfa1","beta1","alfa2","beta2","d");
  hev[iev]->Fit("twoGammaProfile");

  double par0 = func->GetParameter(0);
  double par1 = func->GetParameter(1);
  double par2 = func->GetParameter(2);
  double par3 = func->GetParameter(3);
  double par4 = func->GetParameter(4);
  double par5 = func->GetParameter(5);

//  double par0 = func->GetParError(0);
//  double par0 = func->GetChisquare();
//  double par0 = func->GetNDF();

  alpha[iev] = par1;
//   cout << " par0 = " << par0 << endl;

  h10->Fill(par1,par2);
  h11->Fill(par3,par4);
/*
  form1 = new TFormula("form1","([2]*pow(x*[1],[0]-1.)*exp(-x*[1])/(ROOT::Math::tgamma([0])))");
  gamma1 = new TF1("gamma1","form1",0,120);
  gamma2 = new TF1("gamma2","form1",0,120);
  gamma1->SetParameters(par1,par2,par0);
  gamma2->SetParameters(par3,par4,par5);
*/

  TCanvas *c1 = new TCanvas();
  c1->Divide(1,1);
  gStyle->SetOptFit(kTRUE);  


//  h1->GetYaxis()->SetRangeUser(0,1.0);
//  hev[iev]->SetAxisRange(0,0.06,"Y");
  hev[iev]->SetAxisRange(0,1.15*maxCont,"Y");
//  func->Draw();
//  hev[iev]->Draw("same");
  hev[iev]->Draw();

//  gamma1->SetLineColor(kYellow);
//  gamma1->Draw("same");
//  gamma2->SetLineColor(kBlue);
//  gamma2->Draw("same");
  func->Draw("same");


//  }

/*
  TLatex *t5 = new TLatex();
  t5->SetNDC();
  t5->SetTextAlign(22);
  t5->SetTextFont(63);
  t5->SetTextSizePixels(22);

  char name1[50];
  snprintf(name1,50,"#alpha_{h} = %2.4f #pm %2.4f",par1,par2);
  t5->DrawLatex(0.7,0.8,name1);
*/
}

  TCanvas *c2 = new TCanvas();
  c2->Divide(2,1);

  c2->cd(1);
  h10->SetOption("box");;
  h10->Draw();

  c2->cd(2);
  h11->SetOption("box");;
  h11->Draw();

//
  TCanvas *c3 = new TCanvas();
  c3->Divide(3,3);

  TLatex *t5 = new TLatex();
  t5->SetNDC();
  t5->SetTextAlign(22);  
  t5->SetTextFont(63);
//  t5->SetTextSizePixels(22);
  t5->SetTextSizePixels(10);

  char name1[50];
  for (int ip=11; ip<20; ++ip) { 
    c3->cd(ip-10);
    hev[ip]->Draw();
//    snprintf(name1,50,"#alpha_{h} = %2.4f #pm %2.4f",alpha[ip],alpha[ip]);
//    t5->DrawLatex(0.65,0.75,name1);
  }

}







