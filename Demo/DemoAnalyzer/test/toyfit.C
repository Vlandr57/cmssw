// c++ classes
#include <iostream>
#include <vector>

//RooFit
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooGamma.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooHist.h"
#include "RooMCStudy.h"
#include "RooFitResult.h"
#include "TMatrixDSym.h"
//#include "RooCmdArg.h"

//ROOT
#include "TCanvas.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "TH2F.h"
#include "TString.h"
#include "TROOT.h"
#include "TStyle.h"


using namespace std;
using namespace RooFit;

void set_style();

int main(){


  RooRealVar zero("zero","",0,0,150);

  //Number of Events and particles per shower
  unsigned int events = 10;
  unsigned int particles = 250;

  //input parameters for the distribution
  double alpha_double = 1.836;
  double beta_double = 0.209;
  double lambda_double = 3.686;
  double theta_double = 2.066;
  RooRealVar fsig("fsig","signal fraction",0.75,0.,1.) ;
  //RooRealVar fsig("fsig","signal fraction",0.3) ;

  double beta_gflash_factor = 1.49;
  double theta_gflash_factor = 16.42;

  beta_double = beta_double/beta_gflash_factor;
  theta_double = theta_double/theta_gflash_factor;
  
  // --- Observable ---
  RooRealVar depth("depth","depth (cm)",0,0,150) ;
  
  // --- Gamma pdf ---
  RooRealVar alpha("alpha","alpha",alpha_double,0.,100.);
  alpha.removeMax(); // set infinite range
  RooRealVar beta("beta","beta",1./beta_double,0,100.) ;
  beta.removeMax(); // set infinite range
  //RooRealVar mu("mu","mu",1.49); // no ranges means variable is fixed in fit
  RooRealVar mu("mu","mu",0); // no ranges means variable is fixed in fit

  RooRealVar lambda("lambda","lambda",lambda_double,0.,100.);
  lambda.removeMax(); // set infinite range
  RooRealVar theta("theta","theta",1./theta_double,0.,100.) ;
  theta.removeMax(); // set infinite range
  RooRealVar nu("nu","nu",16.42); // no ranges means variable is fixed in fit

  RooGamma gamma_one("gamma_one","gamma_one pdf",depth,alpha,beta,mu) ;
  RooGamma gamma_two("gamma_two","gamma_two pdf",depth,lambda,theta,mu) ;

  RooAddPdf model("model","model",RooArgList(gamma_one,gamma_two),fsig) ;
  

  // Do the MonteCarlo study with RooMCStudy
  RooMCStudy* study_model = new RooMCStudy(model,depth,Silence(),FitOptions(Save(kTRUE)));
  study_model->generateAndFit(events,particles);


  //generate data for one histogram
  RooDataSet * data = model.generate(depth,particles);
  //RooDataSet * data = gamma_one.generate(depth,1000);
  //fsig.setVal(0.);
  model.fitTo(*data);
  //gamma_two.fitTo(*data);

  //get the correlation Matrix for the first generated and fitted histograms
  TMatrixDSym correlationMatrix = study_model->fitResult(0)->correlationMatrix();

  correlationMatrix.Print();

  //prepare to plot pull histograms
  
  RooPlot* pull_alpha = study_model->plotPull(alpha,Bins(20),FitGauss(kTRUE)) ;
  RooPlot* pull_beta = study_model->plotPull(beta,Bins(40),FitGauss(kTRUE)) ;
  RooPlot* pull_lambda = study_model->plotPull(lambda,Bins(40),FitGauss(kTRUE)) ;
  RooPlot* pull_theta = study_model->plotPull(theta,Bins(40),FitGauss(kTRUE)) ;
  RooPlot* pull_c = study_model->plotPull(fsig,Bins(40),FitGauss(kTRUE)) ;
  

  //prepare to plot one fit
  
  RooPlot * onegammeFrame = depth.frame(mu.getVal(),150,100);
  data->plotOn(onegammeFrame);
  gamma_two.plotOn(onegammeFrame);
  

  RooPlot * depthFrame = depth.frame(mu.getVal(),150,100);
  zero.plotOn(depthFrame,LineStyle(kDashed),LineColor(kBlack));
  data->plotOn(depthFrame);
  model.plotOn(depthFrame,Components(gamma_one),LineStyle(kDashed),LineColor(kGreen),Name("gamma_one"));
  model.plotOn(depthFrame,Components(gamma_two),LineStyle(kDashed),LineColor(kRed),Name("gamma_two"));
  model.plotOn(depthFrame,Name("model"));
  
  cout << alpha.getVal() << " " << 1/beta.getVal()*beta_gflash_factor << std::endl;
  cout << lambda.getVal() << " " << 1/theta.getVal()*theta_gflash_factor << std::endl;
  cout << fsig.getVal() << " "<<depthFrame->chiSquare()<< std::endl;


  RooHist* hpull = depthFrame->pullHist();
  hpull->SetLineColor(0);
  hpull->SetMarkerColor(kRed);
  hpull->SetMarkerStyle(21);

  depthFrame->addPlotable(hpull);
  depthFrame->SetMinimum(-5);
  

  RooPlot * pullhistFrame = depth.frame(Title("Pull Distribution"));
  pullhistFrame->addPlotable(hpull);
  pullhistFrame->SetMaximum(5);
  pullhistFrame->SetMinimum(-5);
  
  //prepare to dump all the plots in the ps file
  TCanvas * can = new TCanvas("can", "can", 600, 600);
  can->cd();
  set_style();
  
  TLegend * legend = new TLegend(0.55,0.65,.9,0.9);
  legend->SetTextFont(72);
  legend->SetTextSize(0.04);
  legend->SetFillColor(kWhite);
  legend->AddEntry(depthFrame->findObject("gamma_one"),"gamma one","L");
  legend->AddEntry(depthFrame->findObject("gamma_two"),"gamma two","L");
  legend->AddEntry(depthFrame->findObject("model")," gamma","L");
  
  can->Print("test.ps[");

  //fill 2D-Histograms for the variables of the diagram and the correlation
  TString variables[] = {"fsig","alpha","beta","lambda","theta"};


  //vector<TH2F *> correllation_hists;

  for(unsigned int i =0; i< 4444; ++i){
    for(unsigned int p =1+i; p < 5;++p){	
     double a_max = 50;
     if(i==0) a_max = 1;
     
     TH2F* corre_hist = new TH2F(variables[i]+variables[p],variables[i]+' '+variables[p],150,0,a_max,150,0,50);
     TH2F* corre_hist_corVal = new TH2F(variables[i]+variables[p]+"_corVal",variables[i]+' '+variables[p],150,0,a_max,150,0,50);
     for(unsigned int m= 0; m <events; m++){
       corre_hist->Fill(((RooRealVar*)study_model->fitResult(m)->floatParsFinal().find(variables[i]))->getVal(),((RooRealVar*)study_model->fitResult(m)->floatParsFinal().find(variables[p]))->getVal());
       corre_hist_corVal->Fill(((RooRealVar*)study_model->fitResult(m)->floatParsFinal().find(variables[i]))->getVal(),((RooRealVar*)study_model->fitResult(m)->floatParsFinal().find(variables[p]))->getVal(),study_model->fitResult(m)->correlation(variables[i],variables[p]));
     }
     corre_hist->Draw("box");
     can->Print("test.ps");
     corre_hist_corVal->Draw("colz");
     can->Print("test.ps");
     //correllation_hists.push_back(corre_hist);
    }
  }	

  //pull distibrutions are drawn
  pull_alpha->Draw();
  can->Print("test.ps");

  pull_beta->Draw();
  can->Print("test.ps");

  pull_lambda->Draw();
  can->Print("test.ps");

  pull_theta->Draw();
  can->Print("test.ps");

  pull_c->Draw();
  
  can->Print("test.ps");


  // one data sample and its fit is drawn
  depthFrame->Draw();
  legend->Draw();
  can->Print("test.ps");

  //pullhistFrame->Draw();
  //onegammeFrame->Draw();
  //can->Print("test.ps");

  can->Print("test.ps]");
 

  return 0;
}

void set_style()
{
  // general appearance and style
  gROOT->SetStyle("Plain");
  
  //gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(kWhite);
  gStyle->SetPadGridX(false);
  gStyle->SetPadGridY(false);
  gStyle->SetGridColor(0);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);
    
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(1);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameFillStyle(0);
  gStyle->SetFrameLineColor(1);
  gStyle->SetFrameLineStyle(1);
  gStyle->SetFrameLineWidth(1);
    
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
    
  gStyle->SetAxisColor(1, "XYZ");
  gStyle->SetStripDecimals(kTRUE);
  gStyle->SetTickLength(0.03, "XYZ");
  gStyle->SetNdivisions(510, "XYZ");
  
  gStyle->UseCurrentStyle();

}
