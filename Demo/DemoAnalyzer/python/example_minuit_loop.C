{
#include "TMath.h"
#include <cmath>

// Version of TTree ROOT analize 
// and using MINUIT for fitting
// =============================

gROOT->Reset();
//gStyle->SetOptFit(0001); 
//gStyle->SetOptFit(1111); 
//gStyle->SetOptFit(kTRUE); 

// Open input file
//================

   TFile *f = new TFile("example_pi50Gev_new.root");
   TTree* t1 = (TTree*)f->Get("demo/Event");
   TTree* t2 = (TTree*)f->Get("demo/Hits");

// Set variables
//==============

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

   const int n = t1->GetEntries();
   cout << " Number of events =  " << n << endl;

// Histos block
//=============

   TH1D *hev[n];
   char name[7];

   int    nBins = 24;
//   int    nBins = 36;
   double maxLength = 120.0;
   double z_depth;
   double e_dep;
   int nevShow = 0;
   int maxShow = 20; // display fit results only for the first maxShow events

   TH1D *h1 = new TH1D("h1","alpha " , 100, 0.0, 12.0);
   TH1D *h2 = new TH1D("h2","beta "  ,  50, 0.0,  1.0);
   TH1D *h3 = new TH1D("h3","lambda" , 100, 0.0, 140.0);
   TH1D *h4 = new TH1D("h4","theta"  ,  50, 0.0, 2.5);
   TH1D *h5 = new TH1D("h5","fsig"   ,  50, 0.0, 1.0);
   TH1D *h15= new TH1D("h15","pi0"   ,  50, 0.0, 1.0);

   TH1D *h6 = new TH1D("h6","alfa vs n" ,   100, 0.0, 100.0);
   TH1D *h7 = new TH1D("h7","beta vs n" ,   100, 0.0, 100.0);
   TH1D *h8 = new TH1D("h8","lambda vs n" , 100, 0.0, 100.0);
   TH1D *h9 = new TH1D("h9","theta vs n" ,  100, 0.0, 100.0);

   TH2D *h10 = new TH2D("h10","alfa vs beta" ,  100,0.0,12.0, 100,0.0,1.0);
   TH2D *h11 = new TH2D("h11","lambda vs theta",100,0.0,140.0, 50,0.0,2.5);

// Start loop over events
//========================

//   for (int iev=0; iev<n; ++iev) {
   for (int iev=0; iev<20; ++iev) {

     t1->GetEntry(iev);
     t2->GetEntry(iev);

     snprintf(name,7,"#hev%d",iev);
     hev[iev] = new TH1D(name,"Longitudinal shower profile",nBins,0.0,maxLength);
     hev[iev]->Sumw2();

// this example for the deposited energy in HE (mytype=3,4)
//--------------------------------------------------------
     e_dep = 0;
     for (int ih=0; ih<myenergy->size(); ++ih) {
        z_depth = ( (*myz)[ih]-interactPoint)/10.;        // in [cm] 
        if( ((*mytype)[ih]>=3 || (*mytype)[ih]<=4) && 
               z_depth>=0.0 && z_depth<=120.0) e_dep +=(*myenergy)[ih];
     }

     for (int ih=0; ih<myenergy->size(); ++ih) {
        z_depth = ( (*myz)[ih]-interactPoint)/10.; // in [cm]
        if( e_dep>0.0 && ((*mytype)[ih]>=3 || (*mytype)[ih]<=4) && 
                                    z_depth>=0.0 && z_depth<=120.0) 
        if( e_dep>0.0 && ((*mytype)[ih]>=3 || (*mytype)[ih]<=4) )
        hev[iev]->Fill(z_depth,(*myenergy)[ih]/((maxLength/nBins)*e_dep));
//        hev[iev]->Fill(z_depth,(*myenergy)[ih]/(e_dep));
     }


//=====================
//===  Start MINUT  ===
//=====================

// Set Minuti with 6 parameters
// -----------------------------
   TMinuit *minimizer = new TMinuit(10); 

// Select verbose level:
// default :          (58 lines in this test)
//      -1 : minimum  ( 4 lines in this test)
//       0 : low      (31 lines)
//       1 : medium   (61 lines)
//       2 : high     (89 lines)
//       3 : maximum (199 lines in this test)
//---------------------------------------------
   minimizer->SetPrintLevel();

// Registrate the function for minimazation
// -----------------------------------------
   minimizer->SetFCN(myFCN);

   Double_t arglist[10];
   Int_t ierflg = 0;
   arglist[0] = 1;

   minimizer->mnexcm("SET ERR", arglist ,1,ierflg);

// Set starting values and step sizes for parameters 
// --------------------------------------------------

//   Double_t vstart[6] = { 2.50, 0.2, 0.50, 25.0, 0.50, 1.0};
   Double_t vstart[6] = { 2.50, 0.2, 0.50, 25.0, 0.50, 0.2};
   Double_t   step[6] = { 0.10, 0.10, 0.10, 0.1, 0.10, 0.1};

//====================================// Set parameters:
// arg1 - parameter number
// arg2 - parameter name
// arg3 - first guess at parametr value
// arg4 - estimated distance to minimum
// arg5, arg6 - ignore for now
// arg7 - error flag
// -----------------
   minimizer->mnparm(0, "p0", vstart[0], step[0], 0, 0, ierflg);
   minimizer->mnparm(1, "p1", vstart[1], step[1], 0, 0, ierflg);
   minimizer->mnparm(2, "p2", vstart[2], step[2], 0, 0, ierflg);
   minimizer->mnparm(3, "p3", vstart[3], step[3], 0, 0, ierflg);
   minimizer->mnparm(4, "p4", vstart[4], step[4], 0, 0, ierflg);
   minimizer->mnparm(5, "p5", vstart[5], step[5], 0, 0, ierflg);

// Now ready for minimization step
// --------------------------------
//  arglist[0] = 1500;
  arglist[0] = 2000;       // number of itteration
  arglist[1] = 1.;

//  minimizer->mnexcm("MIGRAD", arglist ,2,ierflg);
  minimizer->mnexcm("MIGRAD", arglist ,1,ierflg);

//  minimizer->mnexcm("HESSE", arglist ,2,ierflg);
//  minimizer->mnexcm("SIMPLEX", arglist ,1,ierflg);
//  minimizer->mnexcm("SCAN", arglist ,1,ierflg);
//  minimizer->mnexcm("COMBINED", arglist ,1,ierflg);

// Print results
//==============

  double fParP0, fParP1, fParP2, fParP3, fParP4, fParP5;
  double fErrP0, fErrP1, fErrP2, fErrP3, fErrP4, fErrP5;

  minimizer->GetParameter(0,fParP0,fErrP0);
  minimizer->GetParameter(1,fParP1,fErrP1);
  minimizer->GetParameter(2,fParP2,fErrP2);
  minimizer->GetParameter(3,fParP3,fErrP3);
  minimizer->GetParameter(4,fParP4,fErrP4);
  minimizer->GetParameter(5,fParP5,fErrP5);

  cout << "\nPrint results from minuit\n" << endl;
  cout << " fParP0 = " << fParP0 << " +- " << fErrP0 << "\n" << endl;
  cout << " fParP1 = " << fParP1 << " +- " << fErrP1 << "\n" << endl;
  cout << " fParP2 = " << fParP2 << " +- " << fErrP2 << "\n" << endl;
  cout << " fParP3 = " << fParP3 << " +- " << fErrP3 << "\n" << endl;
  cout << " fParP4 = " << fParP4 << " +- " << fErrP4 << "\n" << endl;
  cout << " fParP5 = " << fParP5 << " +- " << fErrP5 << "\n" << endl;

//  cout << " int = " << hev[iev]->Integral(0.0,120.0) << endl;

  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  minimizer->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

// Skip events with bad fit
//===========================
  if(icstat != 3) continue;
  if( fParP0<0.0 || fParP1<0.0 || fParP3<0.0 || fParP4<0.0 ||
                                   fParP2<0.0 || fParP2>1.0 ) continue;


// Remove mismatching of 1-st and 2-nd gamma-functions
//----------------------------------------------------

//  cout << " 1 = " << fParP0/fParP1
//       << " 2 = " << fParP3/fParP4 << endl;

  if( (fParP0/fParP1 > fParP3/fParP4) ) {

//  if( (fParP0/fParP1>fParP3/fParP4 && fParP0>fParP3) || 
//      (fParP0/fParP1<fParP3/fParP4 && fParP0 > 20.0) ) {

    double t_fParP0 = fParP0;
    double t_fParP1 = fParP1;
    double t_fParP3 = fParP3;
    double t_fParP4 = fParP4;
    fParP0 = t_fParP3;
    fParP1 = t_fParP4;
    fParP3 = t_fParP0;
    fParP4 = t_fParP1;
    fParP2 = 1.0-fParP2;
  }


/*
  cout << "\n";
  cout << " Minimum chi square = " << amin << "\n";
  cout << " Estimated vert. distance to min. = " << edm << "\n";
  cout << " Number of variable parameters = " << nvpar << "\n";
  cout << " Highest number of parameters defined by user = " << nparx << "\n";
  cout << " Status of covariance matrix = " << icstat << "\n";

// minimizer->mnemat(correlation_matrix, number_of_parameters);
 Double_t emat[6][6];
 minimizer->mnemat(&emat[0][0], 6);
// printCorrelationMatrix (minimizer,6);
*/

// Prepare for display the fit functions
//======================================

  form1 = new TFormula("form1","[3]*([2]*pow(x*[1],[0]-1.)*exp(-x*[1])/(ROOT::Math::tgamma([0])))");
  form2 = new TFormula("form2","[3]*((1.-[2])*pow(x*[1],[0]-1.)*exp(-x*[1])/(ROOT::Math::tgamma([0])))");
  form3 = new TFormula("form3","[5]*([2]*pow(x*[1],[0]-1.)*exp(-x*[1])/(ROOT::Math::tgamma([0]))+(1.-[2])*pow(x*[4],[3]-1.)*exp(-x*[4])/(ROOT::Math::tgamma([3])))");

  gamma1 = new TF1("gamma1","form1",0,120);
  gamma2 = new TF1("gamma2","form2",0,120);
  sGamma = new TF1("sGamma","form3",0,120);

  gamma1->SetParameters(fParP0,fParP1,fParP2,fParP5);
  gamma2->SetParameters(fParP3,fParP4,fParP2,fParP5);
  sGamma->SetParameters(fParP0,fParP1,fParP2,fParP3,fParP4,fParP5);

//  cout << " gamma1 = " << gamma1->Integral(0.0,120.0) << endl;
//  cout << " gamma2 = " << gamma2->Integral(0.0,120.0) << endl;
//  cout << " sGamma = " << sGamma->Integral(0.0,120.0) << endl;
//}

//  double newPi0 = (gamma1->Integral(0.0,120.0))/(sGamma->Integral(0.0,120.0));   
//  fParP2 = (gamma1->Integral(0.0,120.0))/(sGamma->Integral(0.0,120.0));   
//  fParP2 = (gamma1->Integral(0.0,140.0))/(sGamma->Integral(0.0,140.0));   

// Fill histograms
//================

   h1->Fill(fParP0);
   h2->Fill(fParP1);
   h3->Fill(fParP3);
   h4->Fill(fParP4);
   h5->Fill(1.0-fParP2);
//   h5->Fill(1.0-newPi0);
   h15->Fill(ePi0tot);

   h6->Fill(float(iev),fParP0);
   h7->Fill(float(iev),fParP1);
   h8->Fill(float(iev),fParP3);
   h9->Fill(float(iev),fParP4);

   h10->Fill(fParP0,fParP1);
   h11->Fill(fParP3,fParP4);

//===============================

//  gStyle->SetOptFit(kTRUE);
//  gStyle->SetOptFit(1111);
//  gStyle->SetOptFit(1);
//   gStyle->SetOptStat(0);
  gStyle->SetOptStat(1);


// Display only first fit results
//===============================
  
  if( nevShow < maxShow ) {
    ++nevShow;  
    Double_t xl1=.05, yl1=0.75, xl2=xl1+.3, yl2=yl1+.125;
    TLegend *leg0 = new TLegend(xl1,yl1,xl2,yl2);

    TCanvas *c1 = new TCanvas();
    c1->Divide(1,1);
    c1->cd(1);

// find maximum bin content
//----------------------------
    double maxCont = 0.0;
    for (Int_t bin=1; bin<hev[iev]->GetSize()-1; ++bin) {
      if( hev[iev]->GetBinContent(bin) > maxCont ) 
                   maxCont = hev[iev]->GetBinContent(bin); 
    }

    TString xaxis="Depth in [cm]";
    hev[iev]->GetXaxis()->SetTitle(xaxis);
    hev[iev]->GetXaxis()->SetLabelSize(0.04);
    hev[iev]->GetXaxis()->SetTitleSize(0.04);

    hev[iev]->SetMarkerColor(kBlue);  
    hev[iev]->SetMarkerStyle(23);
    hev[iev]->SetMarkerSize(1.0);
    hev[iev]->SetAxisRange(0.0,1.15*maxCont,"Y");
    hev[iev]->Draw();

//  gamma1->SetLineColor(kYellow);
    gamma1->SetLineColor(3);
    gamma1->Draw("same");  
    gamma2->SetLineColor(kRed);
    gamma2->Draw("same");  
    sGamma->SetLineColor(kBlue);
    sGamma->SetLineStyle(9);
    sGamma->Draw("same");  

    TLatex *t5 = new TLatex();
    t5->SetNDC();
    t5->SetTextAlign(22);
    t5->SetTextFont(63);
    t5->SetTextSizePixels(22);

    char name1[50],name2[50],name3[50],name4[50],name5[50];
    snprintf(name1,50,"#alpha_{h} = %2.4f #pm %2.4f",fParP0,fErrP0);
    snprintf(name2,50,"#beta_{h} = %2.4f #pm %2.4f",fParP1,fErrP0);

    snprintf(name3,50,"#alpha_{e} = %2.4f #pm %2.4f",fParP3,fErrP3);
    snprintf(name4,50,"#beta_{e} = %2.4f #pm %2.4f",fParP4,fErrP4);

    snprintf(name5,50,"#pi^{0}-fraction = %2.4f #pm %2.4f",1.-fParP2,fErrP2);

    t5->DrawLatex(0.7,0.8,name1);
    t5->DrawLatex(0.7,0.7,name2);
    t5->DrawLatex(0.7,0.6,name3);
    t5->DrawLatex(0.7,0.5,name4);
    t5->DrawLatex(0.65,0.4,name5);

  } // end if condition

 } // end loop over events

  gStyle->SetOptStat(11111111);
//  gStyle->SetOptStat(kTRUE);
//  gStyle->SetOptStat(0);
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
   h15->Draw();

   c2->cd(7);
   h10->SetOption("box");;
   h10->Draw();

   c2->cd(8);
   h11->SetOption("box");;
   h11->Draw();

}


