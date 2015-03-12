///////////////////////////////////////////////////////////////////////////////
// File: FastHFShowerLibrary.cc
// Description: Shower library for Very forward hadron calorimeter
///////////////////////////////////////////////////////////////////////////////

#include "FastSimulation/ShowerDevelopment/interface/FastHFShowerLibrary.h"
#include "FastSimulation/Event/interface/FSimEvent.h"
#include "FastSimulation/Event/interface/FSimTrack.h"
#include "SimG4CMS/Calo/interface/HFFibreFiducial.h"

#include "SimDataFormats/CaloHit/interface/HFShowerLibraryEventInfo.h"
#include "DetectorDescription/Core/interface/DDFilter.h"
#include "DetectorDescription/Core/interface/DDFilteredView.h"
#include "DetectorDescription/Core/interface/DDValue.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "Randomize.hh"
#include "CLHEP/Units/GlobalSystemOfUnits.h"
#include "CLHEP/Units/GlobalPhysicalConstants.h"

// STL headers 
#include <vector>
#include <iostream>

//#define DebugLog

FastHFShowerLibrary::FastHFShowerLibrary(edm::ParameterSet const & p) : fastCalo(p)
{
  edm::ParameterSet m_HS   = p.getParameter<edm::ParameterSet>("HFShowerLibrary");
  applyFidCut              = m_HS.getParameter<bool>("ApplyFiducialCut");
  backProb                 = m_HS.getParameter<double>("BackProbability");

  edm::ParameterSet m_HF  = p.getParameter<edm::ParameterSet>("HFShower");
  probMax                 = m_HF.getParameter<double>("ProbMax");
}

FastHFShowerLibrary::~FastHFShowerLibrary() 
{
 if(hfshower)         delete hfshower;
 if(numberingScheme)  delete numberingScheme;
 if(numberingFromDDD) delete numberingFromDDD;
}

void const FastHFShowerLibrary::initHFShowerLibrary(const edm::EventSetup& iSetup) {

  edm::LogInfo("FastCalorimetry") << "initHFShowerLibrary::initialization"; 

  edm::ESTransientHandle<DDCompactView> cpv;
  iSetup.get<IdealGeometryRecord>().get(cpv);

  std::string name = "HcalHits";
  hfshower         = new HFShowerLibrary(name,*cpv,fastCalo);
  numberingFromDDD = new HcalNumberingFromDDD(name, *cpv);  
  numberingScheme  = new HcalNumberingScheme();

  initRun();  
}

void FastHFShowerLibrary::initRun() {

  geantinoPDG = 0; gammaPDG = 22;
  emPDG   = 11; epPDG    = -11; nuePDG   = 12; anuePDG   = -12;
  numuPDG = 14; anumuPDG = -14; nutauPDG = 16; anutauPDG = -16;
  pi0PDG = 111; etaPDG   = 221;

#ifdef DebugLog
  edm::LogInfo("FastCalorimetry") << "HFShowerLibrary: Particle codes for e- = " 
			   << emPDG << ", e+ = " << epPDG << ", gamma = " 
			   << gammaPDG << ", pi0 = " << pi0PDG << ", eta = " 
			   << etaPDG << ", geantino = " << geantinoPDG 
			   << "\n        nu_e = " << nuePDG << ", nu_mu = " 
			   << numuPDG << ", nu_tau = " << nutauPDG 
			   << ", anti_nu_e = " << anuePDG << ", anti_nu_mu = " 
			   << anumuPDG << ", anti_nu_tau = " << anutauPDG;
#endif
}

void FastHFShowerLibrary::recoHFShowerLibrary(const FSimTrack& myTrack) {

#ifdef DebugLog
  edm::LogInfo("FastCalorimetry") << "FastHFShowerLibrary: recoHFShowerLibrary ";
#endif 

  if(!myTrack.onVFcal()) {
#ifdef DebugLog
    edm::LogInfo("FastCalorimetry") << "FastHFShowerLibrary: we should not be here ";
#endif
  }

  hitMap.clear();
  double eGen  = 1000.*myTrack.vfcalEntrance().e();                // energy in [MeV]
  double delZv = (myTrack.vfcalEntrance().vertex().Z()>0.0) ? 50.0 : -50.0;
  G4ThreeVector vertex( 10.*myTrack.vfcalEntrance().vertex().X(),
                        10.*myTrack.vfcalEntrance().vertex().Y(),
                        10.*myTrack.vfcalEntrance().vertex().Z()+delZv); // in [mm]

  G4ThreeVector direction(myTrack.vfcalEntrance().Vect().X(),
                          myTrack.vfcalEntrance().Vect().Y(),
                          myTrack.vfcalEntrance().Vect().Z());

  bool ok;
  double weight = 1.0;                     // rad. damage 
  int parCode   = myTrack.type();

  std::vector<FastHFShowerLibrary::Hit> hits =
              getHits(vertex, direction, parCode, eGen, ok, weight, false);

  for (unsigned int i=0; i<hits.size(); ++i) {
    G4ThreeVector pos = hits[i].position;
    int depth         = hits[i].depth;
    double time       = hits[i].time;

    if (isItinFidVolume (pos)) {     
      int det = 5;
      int lay = 1;
      uint32_t id = 0;
      HcalNumberingFromDDD::HcalID tmp = numberingFromDDD->unitID(det, pos, depth, lay);
      id = numberingScheme->getUnitID(tmp);

      CaloHitID current_id(id,time,myTrack.id());
      std::map<CaloHitID,float>::iterator cellitr;
      cellitr = hitMap.find(current_id);
      if(cellitr==hitMap.end()) {
         hitMap.insert(std::pair<CaloHitID,float>(current_id,1.0));
      } else {
         cellitr->second += 1.0;
      }
    }  // end of isItinFidVolume check 

  } // end loop over hits

}

bool FastHFShowerLibrary::isItinFidVolume (G4ThreeVector& hitPoint) {
  bool flag = true;
  if (applyFidCut) {
    int npmt = HFFibreFiducial::PMTNumber(hitPoint);
#ifdef DebugLog
    edm::LogInfo("FastCalorimetry") << "HFShowerLibrary::isItinFidVolume:#PMT= " 
                                    << npmt << " for hit point " << hitPoint;
#endif
    if (npmt <= 0) flag = false;
  }
#ifdef DebugLog
    edm::LogInfo("FastCalorimetry") << "HFShowerLibrary::isItinFidVolume: point " 
                                    << hitPoint << " return flag " << flag;
#endif
  return flag;
}


std::vector<FastHFShowerLibrary::Hit> FastHFShowerLibrary::getHits(const G4ThreeVector & hitPoint,
                                  const G4ThreeVector & momDir, int parCode, double pin, 
                                  bool & ok, double weight, bool onlyLong) {

  std::vector<FastHFShowerLibrary::Hit> hit;
  ok = false;
  if (parCode == pi0PDG || parCode == etaPDG || parCode == nuePDG ||
      parCode == numuPDG || parCode == nutauPDG || parCode == anuePDG ||
      parCode == anumuPDG || parCode == anutauPDG || parCode == geantinoPDG) 
    return hit;
  ok = true;

  double tSlice = 0.1*hitPoint.mag()/29.98;
  double pz     = momDir.z(); 
  double zint   = hitPoint.z(); 

// if particle moves from interaction point or "backwards (halo)
  int backward = 0;
  if (pz * zint < 0.) backward = 1;
  
  double sphi   = sin(momDir.phi());
  double cphi   = cos(momDir.phi());
  double ctheta = cos(momDir.theta());
  double stheta = sin(momDir.theta());

#ifdef DebugLog
  edm::LogInfo("FastCalorimetry") << "HFShowerLibrary: getHits " << parCode
			          << " of energy " << pin/GeV << " GeV"
			          << "  dir.orts " << momDir.x() << ", " <<momDir.y() 
			          << ", " << momDir.z() << "  Pos x,y,z = " 
			          << hitPoint.x() << "," << hitPoint.y() << "," 
			          << hitPoint.z() << "," 
			          << " sphi,cphi,stheta,ctheta  = " << sphi 
			          << ","  << cphi << ", " << stheta << "," << ctheta; 
#endif    

  if (parCode == emPDG || parCode == epPDG || parCode == gammaPDG ) {
    if (pin<hfshower->pmom[hfshower->nMomBin-1]) {
      hfshower->interpolate(0, pin);
    } else {
      hfshower->extrapolate(0, pin);
    }
  } else {
    if (pin<hfshower->pmom[hfshower->nMomBin-1]) {
      hfshower->interpolate(1, pin);
    } else {
      hfshower->extrapolate(1, pin);
    }
  }

  int nHit = 0;
  FastHFShowerLibrary::Hit oneHit;
  for (int i = 0; i < hfshower->npe; i++) {
    double zv = std::abs(hfshower->pe[i].z()); // abs local z  

#ifdef DebugLog
    edm::LogInfo("FastCalorimetry") << "HFShowerLibrary: Hit " << i << " " << pe[i] << " zv " << zv;
#endif

    if (zv <= hfshower->gpar[1] && hfshower->pe[i].lambda() > 0 && 
       (hfshower->pe[i].z() >= 0 || (zv > hfshower->gpar[0] && (!onlyLong)))) {
      int depth = 1;
      if (onlyLong) {
      } else if (backward == 0) {    // fully valid only for "front" particles
	if (hfshower->pe[i].z() < 0) depth = 2; // with "front"-simulated shower lib.
      } else {                       // for "backward" particles - almost equal
	double r = G4UniformRand();  // share between L and S fibers
        if (r > 0.5) depth = 2;
      } 

      // Updated coordinate transformation from local
      //  back to global using two Euler angles: phi and theta
      double pex = hfshower->pe[i].x();
      double pey = hfshower->pe[i].y();

      double xx = pex*ctheta*cphi - pey*sphi + zv*stheta*cphi; 
      double yy = pex*ctheta*sphi + pey*cphi + zv*stheta*sphi;
      double zz = -pex*stheta + zv*ctheta;

      G4ThreeVector pos  = hitPoint + G4ThreeVector(xx,yy,zz);
      zv = std::abs(pos.z()) - hfshower->gpar[4] - 0.5*hfshower->gpar[1];
      G4ThreeVector lpos = G4ThreeVector(pos.x(),pos.y(),zv);

      zv = hfshower->getFibre()->zShift(lpos,depth,0);     // distance to PMT !

      double r  = pos.perp();
      double p  = hfshower->getFibre()->attLength(hfshower->pe[i].lambda());
      double fi = pos.phi();
      if (fi < 0) fi += twopi;
      int    isect = int(fi/hfshower->dphi) + 1;
      isect        = (isect + 1) / 2;
      double dfi   = ((isect*2-1)*hfshower->dphi - fi);
      if (dfi < 0) dfi = -dfi;
      double dfir  = r * sin(dfi);

#ifdef DebugLog
      edm::LogInfo("FastCalorimetry") << "HFShowerLibrary: Position shift " << xx 
   			              << ", " << yy << ", "  << zz << ": " << pos 
			              << " R " << r << " Phi " << fi << " Section " 
			              << isect << " R*Dfi " << dfir << " Dist " << zv;
#endif

      zz           = std::abs(pos.z());
      double r1    = G4UniformRand();
      double r2    = G4UniformRand();
      double r3    = -9999.;
      if (backward)     r3    = G4UniformRand();
      if (!applyFidCut) dfir += hfshower->gpar[5];

#ifdef DebugLog
      edm::LogInfo("FastCalorimetry") << "HFShowerLibrary: rLimits " << hfshower->rInside(r)
			       << " attenuation " << r1 <<":" << exp(-p*zv) 
			       << " r2 " << r2 << " r3 " << r3 << " rDfi "  
			       << hfshower->gpar[5] << " zz " 
			       << zz << " zLim " << hfshower->gpar[4] << ":" 
			       << hfshower->gpar[4]+hfshower->gpar[1] << "\n"
			       << "  rInside(r) :" << hfshower->rInside(r) 
			       << "  r1 <= exp(-p*zv) :" <<  (r1 <= exp(-p*zv))
			       << "  r2 <= probMax :"    <<  (r2 <= probMax*weight)
			       << "  r3 <= backProb :"   <<  (r3 <= backProb) 
			       << "  dfir > gpar[5] :"   <<  (dfir > hfshower->gpar[5])
			       << "  zz >= gpar[4] :"    <<  (zz >= hfshower->gpar[4])
			       << "  zz <= gpar[4]+gpar[1] :" 
			       << (zz <= hfshower->gpar[4]+hfshower->gpar[1]);   
#endif

      if (hfshower->rInside(r) && r1 <= exp(-p*zv) && r2 <= probMax*weight && 
	  dfir > hfshower->gpar[5] && zz >= hfshower->gpar[4] && 
          zz <= hfshower->gpar[4]+hfshower->gpar[1] && r3 <= backProb && 
          (depth != 2 || zz >= hfshower->gpar[4]+hfshower->gpar[0])) {

	oneHit.position = pos;
	oneHit.depth    = depth;
	oneHit.time     = (tSlice+(hfshower->pe[i].t())+
                           (hfshower->getFibre()->tShift(lpos,depth,1)));
	hit.push_back(oneHit);
#ifdef DebugLog
	edm::LogInfo("FastCalorimetry") << "HFShowerLibrary: Final Hit " << nHit 
				        <<" position " << (hit[nHit].position) 
				        << " Depth " << (hit[nHit].depth) <<" Time " 
				        << tSlice << ":" << hfshower->pe[i].t() << ":" 
				        << hfshower->getFibre()->tShift(lpos,depth,1) << ":" 
				        << (hit[nHit].time);
#endif
	nHit++;
      }
#ifdef DebugLog
      else  LogDebug("FastCalorimetry") << "HFShowerLibrary: REJECTED !!!";
#endif

      if (onlyLong && zz >= hfshower->gpar[4]+hfshower->gpar[0] && 
          zz <= hfshower->gpar[4]+hfshower->gpar[1]) {
	r1    = G4UniformRand();
	r2    = G4UniformRand();
	if (hfshower->rInside(r) && r1 <= exp(-p*zv) && r2 <= probMax && 
            dfir > hfshower->gpar[5]){
	  oneHit.position = pos;
	  oneHit.depth    = 2;
	  oneHit.time     = (tSlice+(hfshower->pe[i].t())+(hfshower->getFibre()->tShift(lpos,2,1)));
	  hit.push_back(oneHit);
#ifdef DebugLog
	  edm::LogInfo("FastCalorimetry") << "HFShowerLibrary: Final Hit " << nHit 
				   << " position " << (hit[nHit].position) 
				   << " Depth " << (hit[nHit].depth) <<" Time "
				   << tSlice << ":" << hfshower->pe[i].t() << ":" 
				   << hfshower->getFibre()->tShift(lpos,2,1) << ":" 
				   << (hit[nHit].time);
#endif
	  nHit++;
	}
      }
    }
  }

#ifdef DebugLog
  edm::LogInfo("FastCalorimetry") << "HFShowerLibrary: Total Hits " << nHit
			          << " out of " << npe << " PE";
#endif

  if (nHit > hfshower->npe && !onlyLong)
    edm::LogWarning("FastCalorimetry") << "HFShowerLibrary: Hit buffer " << hfshower->npe 
				       << " smaller than " << nHit << " Hits";
  return hit;

}


