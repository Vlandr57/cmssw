// c++ classes
#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>

//RooFit
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooGamma.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooHist.h"
#include "RooMCStudy.h"
#include "RooFitResult.h"
#include "TMatrixDSym.h"

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

using namespace RooFit;

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

   const int nEvents(t1->GetEntries());
   cout << " Number of events =  " << nEvents << endl;

   TH1D *h1;
   TH1D *hev[nEvents];
   char name[7];

   int    nBins = 24;
//   int    nBins = 48;
   double maxLength = 120.0;
   double z_depth;
   double e_dep;

   TH1D *h1 = new TH1D("h1","alpha " , 100, 0.0, 10.0);
   TH1D *h2 = new TH1D("h2","beta "  , 100, 0.0, 0.5);
   TH1D *h3 = new TH1D("h3","lambda" , 100, 0.0, 140.0);
   TH1D *h4 = new TH1D("h4","theta"  , 100, 0.0, 5.0);
   TH1D *h5 = new TH1D("h5","fsig"   , 100, 0.0, 1.0);

   TH1D *h6 = new TH1D("h6","alfa vs #event" ,  nEvents, 0.0, nEvents);
   TH1D *h7 = new TH1D("h7","beta vs #event" ,  nEvents, 0.0, nEvents);
   TH1D *h8 = new TH1D("h6","lambda vs #event", nEvents, 0.0, nEvents);
   TH1D *h9 = new TH1D("h7","theta vs #event",  nEvents, 0.0, nEvents);

   TH2D *h10 = new TH2D("h10","alfa vs beta" ,  100,0.0,10.0, 100,0.0,0.5);
   TH2D *h11 = new TH2D("h11","lambda vs theta",100,0.0,140.0,100,0.0,5.0);

   for (int iev=0; iev<nEvents; ++iev) {
//   for (int iev=0; iev<5; ++iev) {

     t1->GetEntry(iev);
     t2->GetEntry(iev);

     snprintf(name,7,"#hev%d",iev);
     hev[iev] = new TH1D(name,"Shower depth ", nBins, 0.0, maxLength);
     hev[iev]->Sumw2();

//     cout << " Number of hits =  " << myenergy->size() << endl;

     e_dep = 0;
     for (int ih=0; ih<myenergy->size(); ++ih) 
       if((*mytype)[ih]>=3 || (*mytype)[ih]<=4) e_dep += (*myenergy)[ih]; // [MeV]

//     cout << " Total energy =  " << e_dep << endl;

     for (int ih=0; ih<myenergy->size(); ++ih) {
        z_depth = ( (*myz)[ih] - interactPoint)/10.; // in [cm]
        if( e_dep>0.0 && ((*mytype)[ih]>=3 || (*mytype)[ih]<=4) )
        hev[iev]->Fill(z_depth,(*myenergy)[ih]/((maxLength/nBins)*e_dep));
//        hev[iev]->Fill(z_depth,(*myenergy)[ih]/(e_dep));
     }

// Declare observable depth
//--------------------------
   RooRealVar depth("depth","depth(cm)",0.,maxLength);

// Create a binned dataset that imports contents of TH1 and associates 
// its contents to observable 'x' -------------------------------------
//-------------------------------
//   RooDataHist data("dh","dh",depth,Import(*h1));
   RooDataHist data("data","dataset with depth",RooArgList(depth),hev[iev]);

// --- pdf ---

  RooRealVar alpha("alpha","alpha",2.5,0.,10.);
  alpha.removeMax();                            // set infinite range
  RooRealVar beta("beta","beta",0.1,0,10.);
  beta.removeMax();                             // set infinite range

  RooRealVar lambda("lambda","lambda",25.0,0.,160.);
  lambda.removeMax();                             // set infinite range
  RooRealVar theta("theta","theta",0.50,0.,10.);
  theta.removeMax();                              // set infinite range

  RooGenericPdf gamma_one("gamma_one","gamma_one pdf",
  "(pow((depth*beta),alpha-1.)*exp(-depth*beta)) / (ROOT::Math::tgamma(alpha))",
                               RooArgList(depth,alpha,beta));

  RooGenericPdf gamma_two("gamma_two","gamma_two pdf",
  "(pow((depth*theta),lambda-1.)*exp(-depth*theta)) / (ROOT::Math::tgamma(lambda))",
                               RooArgList(depth,lambda,theta));

  RooRealVar fsig("fsig","signal fraction",0.70,0.,1.);

  RooAddPdf model("model","model",RooArgList(gamma_one,gamma_two),fsig);

// Define "signal" range in x as [-3,3]
//-------------------------------------
//  depth.setRange("signal",0.,maxLength) ;
//  model.fitTo(data,Save(kTRUE),Range("signal"));
//  model.fitTo(data,Range("signal"));

//  model.fitTo(data,SumW2Error(kTRUE));
  model.fitTo(data);
//  model.fitTo(data,Save());
//  model.fitTo(data, Extended(kTRUE));

// set.writeToFile(“config.txt”);
// set.readFromFile(“config.txt”);

//  alpha.Print();
//  beta.Print();

  h1->Fill(alpha.getVal());
  h2->Fill(beta.getVal());
  h3->Fill(lambda.getVal());
  h4->Fill(theta.getVal());
  h5->Fill(fsig.getVal());

  h6->Fill(float(iev),lambda.getVal());
  h7->Fill(float(iev),theta.getVal());

  h10->Fill(alpha.getVal(),beta.getVal());
  h11->Fill(lambda.getVal(),theta.getVal());

// Plot data and PDF overlaid
//----------------------------

  RooPlot* xframe = depth.frame();
  data.plotOn(xframe);
  model.plotOn(xframe);
  model.plotOn(xframe,Components(gamma_one),LineStyle(kDashed),LineColor(kGreen),
               Name("gamma_one"));
  model.plotOn(xframe,Components(gamma_two),LineStyle(kDashed),LineColor(kRed),
               Name("gamma_one"));

//  model.paramOn(xframe,Layout(0.55,0.80,0.80));
// sum.paramOn(xframe, Layout(0.6), Format("NEU", AutoPrecision(1)),
// Parameters(RooArgList(bwmean, bwsigma, cbmean, cbsigma, bwsig,
// cbsig, b)));
//  data.statOn(xframe);
//  data.statOn(xframe,Layout(0.55,0.99,0.8));
//  model.plotOn(xframe,LineStyle(kDashed),LineColor(kRed));
//  gamma_one.plotOn(xframe,LineStyle(kDashed),LineColor(kRed));
//  gamma_two.plotOn(xframe,LineStyle(kDashed),LineColor(kGreen));
  xframe->Draw();
    
//}


//  TCanvas *c1 = new TCanvas();
//  gStyle->SetOptFit(kTRUE);
//  c1->Divide(1,1);
//  gStyle->SetOptStat(11111111);
//  c1->cd(1); 
//  hev[iev]->Draw();

}

  gStyle->SetOptStat(10111111);
  TCanvas *c2 = new TCanvas();
  c2->Divide(3,3);

  c2->cd(1);
  h1->Draw();

  c2->cd(2);
  h2->Draw();

  c2->cd(3);
  h3->Draw();

  c2->cd(4);
  h4->Draw();

  c2->cd(5);
  h5->Draw();

  c2->cd(6);
//  h10->SetOption("lego");
  h10->SetOption("box");
  h10->Draw();

  c2->cd(7);
//  h11->SetOption("lego");
  h11->SetOption("box");	
  h11->Draw();

}







