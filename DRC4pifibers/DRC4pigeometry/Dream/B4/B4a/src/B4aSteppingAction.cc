//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: B4aSteppingAction.cc 68058 2013-03-13 14:47:43Z gcosmo $
// 
/// \file B4aSteppingAction.cc
/// \brief Implementation of the B4aSteppingAction class

#include "B4aSteppingAction.hh"
#include "B4aEventAction.hh"
#include "B4DetectorConstruction.hh"
#include "G4Material.hh"
#include "G4UnitsTable.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"
#include <stdlib.h>

#include "G4OpBoundaryProcess.hh"

#include <chrono>
#include <random>

#include "TLorentzVector.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4aSteppingAction::B4aSteppingAction(
                      const B4DetectorConstruction* detectorConstruction,
                      B4aEventAction* eventAction)
  : G4UserSteppingAction(),
    fDetConstruction(detectorConstruction),
    fEventAction(eventAction)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4aSteppingAction::~B4aSteppingAction()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4aSteppingAction::UserSteppingAction(const G4Step* step)
{
	
  // get volume of the current pre-step
  G4VPhysicalVolume* PreStepVolume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  G4double energydeposited = step->GetTotalEnergyDeposit();
  G4double steplength = step->GetStepLength();
  
  //define Birk's constant
  double k_B = 0.126; 
  G4double saturatedenergydeposited = 0.;

  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);

  std::poisson_distribution<int> cher_distribution(1*0.2723);

  std::poisson_distribution<int> scin_distribution(1.);
  
  //Add energy deposited in towers (copper only)  
  if (PreStepVolume->GetName() != "World"){
    fEventAction->Addenergy(energydeposited);
      if (PreStepVolume->GetLogicalVolume()->GetMaterial()->GetName() == "Copper"){
        G4double copynumbertower = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(); 
    	G4double copynumberslice = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1); 
    	if (copynumbertower > 0){          //in barrel right or endcap right
     	  fEventAction->AddVectorR(energydeposited,copynumbertower, copynumberslice);
    	}
  	if (copynumbertower < 0){          //im in barrel left or endcap left
  	  fEventAction->AddVectorL(energydeposited, copynumbertower, copynumberslice);
        }
     }
  }
  
  //compute leak
  if (PreStepVolume->GetName() == "leakageabsorber"){
    auto name = step->GetTrack()->GetDefinition()->GetParticleName();
      if (name=="nu_mu" || name=="nu_e" || name=="anti_nu_e" || name=="anti_nu_mu"){
      	fEventAction->Addneutrinoleakage(step->GetTrack()->GetKineticEnergy());
      	step->GetTrack()->SetTrackStatus(fStopAndKill);
      }
      else{
      fEventAction->Addleakage(step->GetTrack()->GetKineticEnergy());
      step->GetTrack()->SetTrackStatus(fStopAndKill);
      }
  };

  //compute em fraction
  if (PreStepVolume->GetName() != "World" ) {
    if (step->GetTrack()->GetDefinition()->GetParticleName() == "e-" || step->GetTrack()->GetDefinition()->GetParticleName() == "e+"){
      fEventAction->Addem(energydeposited);
    }
  }
 
 //primary particle energy
   if ( step->GetTrack()->GetTrackID() == 1 && step->GetTrack()->GetCurrentStepNumber() == 1){
    fEventAction->SavePrimaryParticle(step->GetTrack()->GetDefinition()->GetParticleName());
    fEventAction->SavePrimaryEnergy(step->GetTrack()->GetVertexKineticEnergy());
  }
		
  //compute and save all informations about scintillating and Cherenkov fibers
  std::string Fiber;
  std::string S_fiber = "fiberCoreScint";
  std::string C_fiber = "fiberCoreChere";
  Fiber = PreStepVolume->GetName(); //name of current step fiber
  
  if ( strstr(Fiber.c_str(),S_fiber.c_str())){          //in a scintillating fiber
    //Function to add up energy depoisted in scintillating fibers:
    //- as signal saturated by Birk's law in VectorSignals
    //- as regular energy deposition in all scintillating fibers in EnergyScin
    if(step->GetTrack()->GetDefinition()->GetPDGCharge() != 0.){
        if (steplength != 0){saturatedenergydeposited = (energydeposited/steplength) / ( 1+k_B*(energydeposited/steplength) ) * steplength;}
    }
    else{saturatedenergydeposited = 0.;}
	  
    fEventAction->AddScin(energydeposited);          //energy deposited in scin fibers (not saturated)
    int copynumbertower = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(2); 
    int copynumberslice = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(3); 
    G4int Sfibercopynumber = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1);
	  
    std::string LengthFibr =  step->GetPreStepPoint()->GetTouchableHandle()->GetVolume(1)->GetName(); 
    int S_fiber_ID = 0;

    scin_distribution = std::poisson_distribution<int> (saturatedenergydeposited*12.5); //smear light emission according to Poissonian distribution
    int s_signal = scin_distribution(generator);                                        //S signal
	
    if (copynumbertower > 0){                                                           //in barrel right or endcap right
      fEventAction->AddVectorScinEnergyR(s_signal,copynumbertower, copynumberslice);    //store signal in vector of towers
      fEventAction->AddVectorR(energydeposited, copynumbertower, copynumberslice);      //store energy deposited in vector of towers
       //unique Fiber ID needed: 168750000 is the max of Sfibercopynumber
      S_fiber_ID = Sfibercopynumber+(168750000*copynumberslice);
    }  
	  
    if (copynumbertower < 0){                                                           //in barrel left or endcap left
      fEventAction->AddVectorScinEnergyL(s_signal, copynumbertower, copynumberslice);
      fEventAction->AddVectorL(energydeposited, copynumbertower, copynumberslice);
      //unique Fiber ID needed: 168750000 is the max of Sfibercopynumber
      S_fiber_ID = Sfibercopynumber-(168750000*copynumberslice);
    }
    vector<float> test = {0.0,0.1,0.2};    
    //Fill vector of Fibers
    if (saturatedenergydeposited>0.){
      //local to global transformation
      G4TouchableHandle theTouchable = step->GetPreStepPoint()->GetTouchableHandle();
      G4ThreeVector origin(0.,0.,0.);
      G4ThreeVector zdir(0.,0.,1.);
      G4ThreeVector vectPos = theTouchable->GetHistory()->
      GetTopTransform().Inverse().TransformPoint(origin);
      G4ThreeVector direction = theTouchable->GetHistory()->
      GetTopTransform().Inverse().TransformAxis(zdir);
      G4double lengthfiber = atof(LengthFibr.c_str());
      G4ThreeVector Halffibervect = direction*lengthfiber/2;
      // Fibre tip position
      G4ThreeVector vectPostip = vectPos-Halffibervect;
      // SiPM position	G4ThreeVector SiPMvecPos = vectPos+Halffibervect;
      G4ThreeVector SiPMvecPos = vectPos+Halffibervect;
      if (s_signal>0.0){
        //fEventAction->WriteFiber_Info(S_fiber_ID,s_signal,1,vectPostip,copynumberslice,copynumbertower);//1 == S 0 == C
	//create vector with photons time of arrival
	//calculate distance from SiPM (mm)
	double distance = sqrt((SiPMvecPos[0]-step->GetTrack()->GetPosition().getX())*(SiPMvecPos[0]-step->GetTrack()->GetPosition().getX())+(SiPMvecPos[1]-step->GetTrack()->GetPosition().getY())*(SiPMvecPos[1]-step->GetTrack()->GetPosition().getY())+(SiPMvecPos[2]-step->GetTrack()->GetPosition().getZ())*(SiPMvecPos[2]-step->GetTrack()->GetPosition().getZ()));
	//calculate time of arrival 
	const float speed_s_fiber = 299.792458/1.59;     //mm/ns
	double time = distance/speed_s_fiber;
	vector<double> phtimes (s_signal, time);          //create vector of n identical times of arrival
	fEventAction->appendfiber(S_fiber_ID, 1, copynumberslice, copynumbertower, s_signal, vectPostip, phtimes);
	// Extract info for z time
	//std::ofstream TimeFile;
	//TimeFile.open("Time.txt", std::ios_base::app);
	//TimeFile<<"Scin "<< std::fixed << std::setprecision(3) <<S_fiber_ID<<" "<<vectPostip.getX()<<" "<<vectPostip.getY()<<" "<<vectPostip.getZ()<<" "<<s_signal<<" "<<sqrt((SiPMvecPos[0]-step->GetTrack()->GetPosition().getX())*(SiPMvecPos[0]-step->GetTrack()->GetPosition().getX())+(SiPMvecPos[1]-step->GetTrack()->GetPosition().getY())*(SiPMvecPos[1]-step->GetTrack()->GetPosition().getY())+(SiPMvecPos[2]-step->GetTrack()->GetPosition().getZ())*(SiPMvecPos[2]-step->GetTrack()->GetPosition().getZ()))<<" "<<step->GetTrack()->GetGlobalTime()<<G4endl;
	//TimeFile.close();
      } 
    } 
  }

  if ( strstr(Fiber.c_str(),C_fiber.c_str())){          //in a Cherenkov fiber
    fEventAction->AddCher(step->GetTotalEnergyDeposit());
  
    G4double copynumbertower = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(2); 
    G4double copynumberslice = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(3);

    if (copynumbertower>0) {fEventAction->AddVectorR(energydeposited, copynumbertower, copynumberslice);}
    if (copynumbertower<0) {fEventAction->AddVectorL(energydeposited, copynumbertower, copynumberslice);}
  }

  //part for cherenkov photons
  G4OpBoundaryProcessStatus theStatus = Undefined;

  G4ProcessManager* OpManager =
                    G4OpticalPhoton::OpticalPhoton()->GetProcessManager();

  if (OpManager) {
    G4int MAXofPostStepLoops =
    OpManager->GetPostStepProcessVector()->entries();
    G4ProcessVector* fPostStepDoItVector =
    OpManager->GetPostStepProcessVector(typeDoIt);

    for ( G4int i=0; i<MAXofPostStepLoops; i++) {
      G4VProcess* fCurrentProcess = (*fPostStepDoItVector)[i];
      fOpProcess = dynamic_cast<G4OpBoundaryProcess*>(fCurrentProcess);
      if (fOpProcess) { theStatus = fOpProcess->GetStatus(); break;}
    }
   }

  std::string SiPMC = "SiPMC";
  std::string SiPMS = "SiPMS";
  std::string SiPMdetection;

  //If the particle is an optical photon...
  if(step->GetTrack()->GetDefinition()->GetParticleName() == "opticalphoton"){

    switch (theStatus){

      case TotalInternalReflection: 
        Fiber = PreStepVolume->GetName();

        if(strstr(Fiber.c_str(), C_fiber.c_str())){          //in a Cherenkov fibre
	  int copynumbertower = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(2); 
	  int copynumberslice = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(3); 
	  G4int Cfibercopynumber = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1);

	  std::string LengthFibr =  step->GetPreStepPoint()->GetTouchableHandle()->GetVolume(1)->GetName(); 
				
	  int C_fiber_ID = 0;
	  int c_signal = cher_distribution(generator);
          if (copynumbertower>0){          //in barrel right or endcap right
	    fEventAction->AddVectorCherPER(c_signal, copynumbertower, copynumberslice);
	    //unique Fiber ID needed: 168750000 is the max of Cfibercopynumber
	    C_fiber_ID = Cfibercopynumber+(168750000*copynumberslice);
	  }
	  if (copynumbertower<0){          //in barrel left or endcap left 
            fEventAction->AddVectorCherPEL(c_signal, copynumbertower, copynumberslice);
	    //I want unique Fiber ID: 168750000 is the max of Cfibercopynumber
            C_fiber_ID = Cfibercopynumber-(168750000*copynumberslice);
	  }
	  //fEventAction->AddCherenkov(); // add one photoelectron from Cherenkov process in Cherenkov fibers                  
				
	  // Fibers routine: fill the C fibres info 
	  G4TouchableHandle theTouchable = step->GetPreStepPoint()->GetTouchableHandle();
  	  G4ThreeVector origin(0.,0.,0.);
	  G4ThreeVector zdir(0.,0.,1.);
  	  G4ThreeVector vectPos = theTouchable->GetHistory()->
    	  GetTopTransform().Inverse().TransformPoint(origin);
	  G4ThreeVector direction = theTouchable->GetHistory()->
    	  GetTopTransform().Inverse().TransformAxis(zdir);
	  G4double lengthfiber = atof(LengthFibr.c_str());
	  G4ThreeVector Halffibervect = direction*lengthfiber/2;
	  // Fibre tip position
	  G4ThreeVector vectPostip = vectPos-Halffibervect;
	  // SiPM position
	  G4ThreeVector SiPMvecPos = vectPos+Halffibervect;
          if (c_signal>0){
            //fEventAction->WriteFiber_Info(C_fiber_ID,c_signal,0,vectPostip,copynumberslice,copynumbertower);// 1 == S 0 == C
	    //calculate distance from SiPM (mm)
            double distance = sqrt((SiPMvecPos[0]-step->GetTrack()->GetPosition().getX())*(SiPMvecPos[0]-step->GetTrack()->GetPosition().getX())+(SiPMvecPos[1]-step->GetTrack()->GetPosition().getY())*(SiPMvecPos[1]-step->GetTrack()->GetPosition().getY())+(SiPMvecPos[2]-step->GetTrack()->GetPosition().getZ())*(SiPMvecPos[2]-step->GetTrack()->GetPosition().getZ()));
	    //calculate time of arrival 
	    const float speed_s_fiber = 299.792458/1.59;     //mm/ns
	    double time = distance/speed_s_fiber;
	    vector<double> phtimes (c_signal, time);         //create vector of n identical times of arrival
	    fEventAction->appendfiber(C_fiber_ID, 0, copynumberslice, copynumbertower, c_signal, vectPostip, phtimes);
	    //std::ofstream TimeFile;
            //TimeFile.open("Time.txt", std::ios_base::app);
	    //TimeFile<<"Cher "<<std::fixed << std::setprecision(3) <<C_fiber_ID<<" "<<vectPostip.getX()<<" "<<vectPostip.getY()<<" "<<vectPostip.getZ()<<" "<<c_signal<<" "<<sqrt((SiPMvecPos[0]-step->GetTrack()->GetPosition().getX())*(SiPMvecPos[0]-step->GetTrack()->GetPosition().getX())+(SiPMvecPos[1]-step->GetTrack()->GetPosition().getY())*(SiPMvecPos[1]-step->GetTrack()->GetPosition().getY())+(SiPMvecPos[2]-step->GetTrack()->GetPosition().getZ())*(SiPMvecPos[2]-step->GetTrack()->GetPosition().getZ()))<<" "<<step->GetTrack()->GetGlobalTime()<<G4endl;
	    //TimeFile.close();
          }
        step->GetTrack()->SetTrackStatus(fStopAndKill); //kill photon			
        }
      break;

      //case Detection: //in case optical surface is used to detect photons
    
      default: 
        step->GetTrack()->SetTrackStatus(fStopAndKill);
      break;
    }
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
