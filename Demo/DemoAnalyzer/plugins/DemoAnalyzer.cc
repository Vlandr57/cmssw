
// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/CaloHit.h"

#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"

#include "SimDataFormats/HcalTestBeam/interface/PHcalTB06Info.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"


//
// class declaration
//

class DemoAnalyzer : public edm::EDAnalyzer {
   public:
      explicit DemoAnalyzer(const edm::ParameterSet&);
      ~DemoAnalyzer();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

// ----------member data --------------------

      TTree*    tree_hit;
      TTree*    tree_evt;

      std::vector<double> myx;
      std::vector<double> myy;
      std::vector<double> myz;
      std::vector<double> myenergy;
      std::vector<double> mytype;

      double energyGen, phiGen, etaGen;
      double startingPoint, interactPoint;
      double ePi0first, ePi0tot;

      G4RotationMatrix*  beamline_RM;
      G4ThreeVector      xyzPosition, pvXYZ;
      G4ThreeVector      uvwPosition, pvUVW;
      G4ThreeVector      hitPosition, hitXYZ;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//

DemoAnalyzer::DemoAnalyzer(const edm::ParameterSet& iConfig) 
{
}

DemoAnalyzer::~DemoAnalyzer()
{
}

// -- method called once each job just before starting event loop  ---

void DemoAnalyzer::beginJob()
{

// Histos initilization
//----------------------     
  edm::Service<TFileService> fs;

  tree_evt = fs->make<TTree>("Event", "Event info");
  tree_evt->Branch("energyGen",&energyGen,"energyGen/D");
  tree_evt->Branch("phiGen",&phiGen,"phiGen/D");
  tree_evt->Branch("etaGen",&etaGen,"phiGen/D");
  tree_evt->Branch("startingPoint",&startingPoint,"startingPoint/D");
  tree_evt->Branch("interactPoint",&interactPoint,"interactPoint/D");
  tree_evt->Branch("ePi0first",&ePi0first,"ePi0first/D");
  tree_evt->Branch("ePi0tot",&ePi0tot,"ePi0tot/D");

  tree_hit = fs->make<TTree>("Hits", "Hits info");
  tree_hit->Branch("myx","std::vector<double>",&myx);
  tree_hit->Branch("myy","std::vector<double>",&myy);
  tree_hit->Branch("myz","std::vector<double>",&myz);
  tree_hit->Branch("myenergy","std::vector<double>",&myenergy);
  tree_hit->Branch("mytype","std::vector<double>",&mytype);

}

//
// member functions
//
// ------------ method called to for each event  ------------

void 
DemoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
//   using namespace reco;

// Put sero all variables
//-------------------------

   energyGen = phiGen = etaGen = 0.0;
   startingPoint = interactPoint = 0.0;
   ePi0first = ePi0tot = 0.0; 

// startingPoint - z-coordinate of the starting point in Hcal or Ecal  
//                 along particle direction
// interactPoint - z-coordinate of the first inelastic interaction point 
//                 in Hcal or Ecal along particle direction
// ePi0first     - energy of all pi-zero in the first inelastic interaction 
//
// ePi0tot       - total energy of all pi-zero in this event 
//
// xyzPosition   - starting point in Hcal or Ecal in global coordinate 
// pvXYZ         - starting point in Hcal or Ecal in local  coordinate
//
// uvwPosition   - first inelastic interaction point in Hcal or Ecal 
//                 in global coordinate
// pvUVW         - first inelastic interaction point in Hcal or Ecal  
//                 in local  coordinate
// hitPosition   - hit position in global coordinate system 
// hitXYZ        - hit position in local coordinate system

// Get information from PHcalTB06Info-branch
//------------------------------------------
   Handle<PHcalTB06Info> testhits;
   iEvent.getByLabel("g4SimHits",testhits);
   bool foundH = iEvent.getByLabel("g4SimHits",testhits);

   if (foundH) {

     PHcalTB06Info info(*testhits.product());
     std::vector<PHcalTB06Info::Hit> calohits = testhits->simHits();   

     double beamE    = info.initE();
     double beamThet = 2*atan(exp(-info.eta()));
     double beamPhi  = (info.phi() < 0.0) ? info.phi()+twopi : info.phi();

     beamline_RM = new G4RotationMatrix;
     beamline_RM->rotateZ(-beamPhi);
     beamline_RM->rotateY(-beamThet);

/*
     std::cout << std::endl;

     std::cout << " PHcalTB06Info testhits is found "
               << std::endl;

     std::cout << " Beam ID = " << info.partID()
               << " energy = "  << info.initE()
               << " eta = " << info.eta()
               << " phi = " << info.phi()  
               << std::endl;

     std::cout << " beamPhi = " << beamPhi
               << " beamThet = " << beamThet << std::endl;
     std::cout << " Rotation matrix = " << *beamline_RM << std::endl;

     std::cout << " Hcal first point = " 
               << " xP = " << info.vtxPrimX() 
               << " yP = " << info.vtxPrimY() 
               << " zP = " << info.vtxPrimZ() << std::endl;

     std::cout << " 1-st interaction pi0 energy = " << info.vtxEpi0f() 
               << " Total pi0 energy = " << info.vtxEpi0t()
               << std::endl; 
*/

     if ( beamE > 0.0 ) {
        ePi0first = info.vtxEpi0f()/beamE;
        ePi0tot   = info.vtxEpi0t()/beamE;
     }

     xyzPosition = G4ThreeVector(info.vtxPrimX(),info.vtxPrimY(),info.vtxPrimZ());
     pvXYZ       = (*beamline_RM)*(xyzPosition);

/*
     std::cout << " Local coord = " << info.vtxType()
               << " xU = " << pvXYZ.x()
               << " yV = " << pvXYZ.y()
               << " zW = " << pvXYZ.z() << std::endl;
*/
//     xfirst = pvXYZ.x();
//     yfirst = pvXYZ.y();
//     zfirst = pvXYZ.z();
     startingPoint = pvXYZ.z();
      
     uvwPosition = G4ThreeVector(info.vtxPrimU(),info.vtxPrimV(),info.vtxPrimW());
     pvUVW      = (*beamline_RM)*(uvwPosition);

/*
     std::cout << " Nuclear interaction point = "
               << " xP = " << info.vtxPrimU()
               << " yP = " << info.vtxPrimV()
               << " zP = " << info.vtxPrimW() << std::endl;

     std::cout << " Local coordinate = " 
               << " xU = " << pvUVW.x() 
               << " yV = " << pvUVW.y()
               << " zW = " << pvUVW.z() << std::endl;
*/

//     xinel = pvUVW.x();
//     yinel = pvUVW.y();
//     zinel = pvUVW.z();
     interactPoint = pvUVW.z();

// save
//-----------

     energyGen = beamE;
     phiGen    = beamPhi;
     etaGen    = beamThet;

     tree_evt->Fill();

     std::cout << " PHcaloTB06Info simHits = " << calohits.size() << std::endl;

// loop over hits and fill energy each hits:
//-------------------------------------------

     for(unsigned int i=0; i<calohits.size();++i) {

       hitPosition = G4ThreeVector(calohits[i].x,calohits[i].y,calohits[i].z);
       hitXYZ      = (*beamline_RM)*(hitPosition);
        
       myx.push_back(hitXYZ.x());
       myy.push_back(hitXYZ.y());
       myz.push_back(hitXYZ.z());
       myenergy.push_back(calohits[i].e);
       mytype.push_back(calohits[i].t);

//       myx.push_back(calohits[i].z);
//       myy.push_back(calohits[i].y);
//       myz.push_back(calohits[i].x);

/*
       std::cout << " hit = " << i
                 << " xx = " << hitPosition.x()
                 << " yy = " << hitPosition.y()
                 << " zz = " << hitPosition.z()
                 << std::endl;

       std::cout << " Local coordinate = "
                 << " xx = " << hitXYZ.x() 
                 << " yy = " << hitXYZ.y()
                 << " zz = " << hitXYZ.z() 
                 << std::endl;

       std::cout << " id  = " << calohits[i].id 
                 << " eta = " << calohits[i].eta
                 << " phi = " << calohits[i].phi 
                 << " e   = " << calohits[i].e 
                 << " l   = " << calohits[i].lay 
                 << std::endl;
*/
     }      
  }
  else {
     std::cout << " PHcaloTB06Info hcalhits is not found "
               << std::endl;
  }
  tree_hit->Fill();  

  myx.erase(myx.begin(), myx.end());
  myy.erase(myy.begin(), myy.end());
  myz.erase(myz.begin(), myz.end());
  myenergy.erase(myenergy.begin(), myenergy.end());
  mytype.erase(mytype.begin(), mytype.end());

// Get information for standard HCAL sim hits
//-------------------------------------------
  Handle<PCaloHitContainer> hcalhits;
  iEvent.getByLabel("g4SimHits","HcalHits",hcalhits);
  bool foundB = iEvent.getByLabel("g4SimHits","HcalHits",hcalhits);

  if (foundB) {

//    std::cout << " PCaloHitContainer hcalhits is found " << std::endl;

    PCaloHitContainer::const_iterator it=hcalhits.product()->begin();
    PCaloHitContainer::const_iterator itend=hcalhits.product()->end();

    for( ; it!=itend; ++it ) {
  
      DetId MyDetId = it->id();
/*
      std::cout << " Hcal MyDetId = " << it->id()
           << " det = "    << MyDetId.det()
           << " subdet = " << MyDetId.subdetId()
           << std::endl;
*/
      if(MyDetId.det()!=1 || MyDetId.subdetId()!=0) continue;
//      if(MyDetId.det()!=4 || MyDetId.subdetId()!=1) continue;
//      if(MyDetId.det()!=4 || MyDetId.subdetId()!=2) continue;

/*
      std::cout << " Hcal prints id = " << it->id() 
           << " subdet = " << myDetId.subdetId()
           << " st_phi = " << myDetId.iphi()
           << " st_eta = " << myDetId.ieta()
//           << " depth = " << myDetId.depth()
           << " energy = " << it->energy()
           << std::endl;
*/
    }
  }

// Get information for standard ECAL (EB) sim hits
//------------------------------------------------
  
  Handle<PCaloHitContainer> ecalhits;   
  iEvent.getByLabel("g4SimHits","EcalHitsEB",ecalhits);
//  bool foundE = iEvent.getByLabel("g4SimHits","EcalHitsEE",ecalhits); 
  bool foundE = iEvent.getByLabel("g4SimHits","EcalHitsEB",ecalhits); 
 
  if (foundE) {
    PCaloHitContainer::const_iterator jt=ecalhits.product()->begin();
    PCaloHitContainer::const_iterator jtend=ecalhits.product()->end();

    for( ; jt!=jtend; ++jt ) {
  
      DetId Eydetid = jt->id();
/*
     std::cout << " Ecal prints Det = " << Eydetid.det()
               << " subdet = " << Eydetid.subdetId()
               << std::endl;
*/
      if( Eydetid.det()!=3 || Eydetid.subdetId()!=1 ) continue;      // EB
//      if( Eydetid.det()!=3 || Eydetid.subdetId()!=2 ) continue;    // EE

//      EEDetId eeId (jt->id());
      EBDetId ebId (jt->id());

/*
     std::cout << " EB subdet = " << ebId.subdetId()
               << " ieta = " << ebId.ieta()
               << " iphi = " << ebId.iphi()
               << std::endl;

     std::cout << " ix = " << eeId.ix() << " iy = " << eeId.iy()
               << " ic = " << eeId.ic() << " isc = " << eeId.isc()
               << " iq = " << eeId.iquadrant()
               << std::endl;
*/

    }
  }

}

// --- method called once each job just after ending the event loop  ------

void DemoAnalyzer::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(DemoAnalyzer);

