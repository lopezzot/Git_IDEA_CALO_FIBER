///
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
// $Id: B4aEventAction.cc 75604 2013-11-04 13:17:26Z gcosmo $
// 
/// \file B4aEventAction.cc
/// \brief Implementation of the B4aEventAction class

#include "B4aEventAction.hh"
#include "B4RunAction.hh"
#include "B4Analysis.hh"
#include "B4PodioManager.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"
#include "G4AutoLock.hh"

#include "Randomize.hh"
#include <iomanip>
#include <vector>

// podio includes

#include "podio/ROOTWriter.h"

namespace { G4Mutex B4aEventActionMutex = G4MUTEX_INITIALIZER; }

B4aEventAction::B4aEventAction() //definition of B4aEvenctAction contructor
 : G4UserEventAction(),
   Energyem(0.),                 //data members set to 0.
   neutrinoleakage(0.),
   leakage(0.),
   EnergyScin(0.),
   EnergyCher(0.),
   EnergyTot(0.),
   PrimaryParticleEnergy(0.),
   VectorSignalsR(0.),          //vectors set to 0-size vector
   VectorSignalsL(0.),
   VectorSignalsCherR(0.),
   VectorSignalsCherL(0.),
   VectorR(0.),
   VectorL(0.),
   VectorR_loop(0.),
   VectorL_loop(0.),
   Fiber_Hits{0.},
   Fibers(0.),                 //vector of Fibers
   FiberIDs(0),                //vector of Fibers IDs
   s_caloHits(NULL),
   c_caloHits(NULL),
   aux_infoHits(NULL),
   s_caloHitContributions(NULL),
   c_caloHitContributions(NULL)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4aEventAction::~B4aEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4aEventAction::BeginOfEventAction(const G4Event* /*event*/)
{  
  //Time_distribution event
  //std::ofstream TimeFile;
  //TimeFile.open("Time.txt", std::ios_base::app);
  //TimeFile<<"Event "<<G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID()<<" % % %"<<G4endl;
  //TimeFile.close();
	
  // initialisation per event
  Energyem = 0.;
  EnergyScin = 0.;
  EnergyCher = 0.;
  EnergyTot = 0.;
  neutrinoleakage = 0.;
  leakage = 0.;
  PrimaryParticleEnergy = 0;  
  
  int fNbOfBarrel = 40;
  int fNbOfEndcap = 35;
  int fNbOfZRot = 36;

  Fibers.clear();

  FiberIDs.clear(); 
	
  for (int i=0;i<1000000;i++) {Fiber_Hits[i]={0};}
	
  for (int i=0;i<VectorR.size();i++){VectorR.at(i)=0.;}

  for (int i=0;i<VectorL.size();i++){VectorL.at(i)=0.;}
	
  for (int i=0;i<VectorR_loop.size();i++){VectorR_loop.at(i)=0.;}

  for (int i=0;i<VectorL_loop.size();i++){VectorL_loop.at(i)=0.;}
       
  for (int i=0;i<VectorSignalsR.size();i++){VectorSignalsR.at(i)=0.;}
    
  for (int i=0;i<VectorSignalsL.size();i++){VectorSignalsL.at(i)=0.;}

  for (int i=0;i<VectorSignalsCherR.size();i++){VectorSignalsCherR.at(i)=0.;}

  for (int i=0;i<VectorSignalsCherL.size();i++){VectorSignalsCherL.at(i)=0.;}
    
  for(int i=0;i<=fNbOfZRot*(fNbOfBarrel+fNbOfEndcap);i++){
    if(VectorR.size() <= fNbOfZRot*(fNbOfBarrel+fNbOfEndcap)){
      VectorR.push_back(0.);
    }
  }

  for(int i=0;i<=fNbOfZRot*(fNbOfBarrel+fNbOfEndcap);i++){
    if(VectorL.size() <= fNbOfZRot*(fNbOfBarrel+fNbOfEndcap)){
      VectorL.push_back(0.);
    }
  }
    
  for(int i=0;i<=fNbOfZRot*(fNbOfBarrel+fNbOfEndcap);i++){
    if(VectorR_loop.size() <= fNbOfZRot*(fNbOfBarrel+fNbOfEndcap)){
      VectorR_loop.push_back(0.);
    }
  }
    
  for(int i=0;i<=fNbOfZRot*(fNbOfBarrel+fNbOfEndcap);i++){
    if(VectorL_loop.size() <= fNbOfZRot*(fNbOfBarrel+fNbOfEndcap)){
      VectorL_loop.push_back(0.);
    }
  }

  for(int i=0;i<=fNbOfZRot*(fNbOfBarrel+fNbOfEndcap);i++){
    if(VectorSignalsR.size() <= fNbOfZRot*(fNbOfBarrel+fNbOfEndcap)){
      VectorSignalsR.push_back(0.);
    }
  }
    
  for(int i=0;i<=fNbOfZRot*(fNbOfBarrel+fNbOfEndcap);i++){
    if(VectorSignalsL.size() <= fNbOfZRot*(fNbOfBarrel+fNbOfEndcap)){
    VectorSignalsL.push_back(0.);
    }
  }

  for(int k=0;k<=fNbOfZRot*(fNbOfBarrel+fNbOfEndcap);k++){
    if(VectorSignalsCherR.size() <= fNbOfZRot*(fNbOfBarrel+fNbOfEndcap)){
      VectorSignalsCherR.push_back(0.);
    }
  }
  
  for(int k=0;k<=fNbOfZRot*(fNbOfBarrel+fNbOfEndcap);k++){
    if(VectorSignalsCherL.size() <= fNbOfZRot*(fNbOfBarrel+fNbOfEndcap)){
      VectorSignalsCherL.push_back(0.);
    }
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4aEventAction::EndOfEventAction(const G4Event* event)
{
  //accumulate statistics over events

  //get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  B4PodioManager * l_podioManager = B4PodioManager::Instance();
  podio::EventStore * l_store = l_podioManager->GetEvtStore();
  podio::ROOTWriter * l_writer = l_podioManager->GetWriter();
    
  //write fiber-by-fiber into file and store in podio Hits
  /*std::ofstream eventFile;
  eventFile.open("Event.txt", std::ios_base::app);
  int v=0;
  G4double E=0.;
  if(G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID()==0) eventFile<<"EvtID\tFiberID\tEt\tXt\tYt\tZt\tFlagt\tslicet\ttowert"<<std::endl;
  while(Fiber_Hits[v].F_ID!=0){
    E = E+Fiber_Hits[v].F_E;

    eventFile<<G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID()<<"\t"<<std::fixed << std::setprecision(3) <<Fiber_Hits[v].F_ID<<"\t"<<Fiber_Hits[v].F_E<<"\t"<<Fiber_Hits[v].F_X<<"\t"<<Fiber_Hits[v].F_Y<<"\t"<<Fiber_Hits[v].F_Z<<"\t"<<Fiber_Hits[v].F_Type<<"\t"<<Fiber_Hits[v].F_slice<<"\t"<<Fiber_Hits[v].F_tower<<std::endl;

    if (Fiber_Hits[v].F_Type == 1){
      auto l_hit = s_caloHits->create();
      l_hit.setCellID((Fiber_Hits[v].F_ID));
      l_hit.setEnergy(Fiber_Hits[v].F_E);
      l_hit.setPosition({Fiber_Hits[v].F_X,Fiber_Hits[v].F_Y,Fiber_Hits[v].F_Z});
    } 
    else if (Fiber_Hits[v].F_Type == 0){
      auto l_hit = c_caloHits->create();              
      l_hit.setCellID((Fiber_Hits[v].F_ID));
      l_hit.setEnergy(Fiber_Hits[v].F_E);
      l_hit.setPosition({Fiber_Hits[v].F_X,Fiber_Hits[v].F_Y,Fiber_Hits[v].F_Z});
    }
    v++;
  }
  eventFile.close();*/

  //order time of arrival vector in Fibers before storing info
  for (unsigned int i = 0; i<Fibers.size(); ++i){Fibers[i].orderphtimes();}  

  //new way to store fiber-by-fiber info in edm4hep/podio
  for (Fiber fiber : Fibers) {
    if (fiber.Type == 1){                         //scintillating
      auto l_hit = s_caloHits->create();
      l_hit.setCellID(fiber.ID);
      l_hit.setEnergy(fiber.E);
      l_hit.setPosition({fiber.Pos.X,fiber.Pos.Y,fiber.Pos.Z});
      //std::cout<<"X: "<<fiber.Pos.X<<" Y: "<<fiber.Pos.Y<< " Z: "<<fiber.Pos.Z<<std::endl;
      // Create the CaloHitContributions
      for (auto ph_time : fiber.phtimes){
	auto l_hitContrib = s_caloHitContributions->create();
	l_hitContrib.setTime(ph_time);
	l_hit.addToContributions(l_hitContrib);
      }
    } 
    else if (fiber.Type == 0){                    //cherenkov
      auto l_hit = c_caloHits->create();              
      l_hit.setCellID(fiber.ID);
      l_hit.setEnergy(fiber.E);
      l_hit.setPosition({fiber.Pos.X,fiber.Pos.Y,fiber.Pos.Z});
      for (auto ph_time : fiber.phtimes){
	auto l_hitContrib = c_caloHitContributions->create();
	l_hitContrib.setTime(ph_time);
	l_hit.addToContributions(l_hitContrib);
      }
      //store me -> fiber.phtimes
    }
  }

  // fill ntuple event by event
  analysisManager->FillNtupleDColumn(0, Energyem);
  analysisManager->FillNtupleDColumn(1, EnergyScin);
  analysisManager->FillNtupleDColumn(2, EnergyCher);
  analysisManager->FillNtupleDColumn(3, EnergyTot);
  analysisManager->FillNtupleDColumn(4, PrimaryParticleEnergy);
  analysisManager->FillNtupleSColumn(5, PrimaryParticleName);
  analysisManager->FillNtupleDColumn(6, neutrinoleakage);
  analysisManager->FillNtupleDColumn(7, leakage);
  analysisManager->AddNtupleRow();                               //columns with vector are automatically filled with this function

  auto l_hit = aux_infoHits->create();
  l_hit.setCellID(0);
  l_hit.setEnergy(Energyem);
  l_hit = aux_infoHits->create();
  l_hit.setCellID(1);
  l_hit.setEnergy(EnergyScin);
  l_hit = aux_infoHits->create();
  l_hit.setCellID(2);
  l_hit.setEnergy(EnergyCher);
  l_hit = aux_infoHits->create();
  l_hit.setCellID(3);
  l_hit.setEnergy(0); //previsouly allocated for NumberofCherenkovDetected (now storing 0, to be removed)
  l_hit = aux_infoHits->create();
  l_hit.setCellID(4);
  l_hit.setEnergy(EnergyTot);
  l_hit = aux_infoHits->create();
  l_hit.setCellID(5);
  l_hit.setEnergy(PrimaryParticleEnergy);
  l_hit = aux_infoHits->create();
  l_hit.setCellID(6);
  l_hit.setEnergy(0); // Not sure what PrimaryParticleEnergy actually is (if needed PDGID can be stored here)
  l_hit = aux_infoHits->create();
  l_hit.setCellID(7);
  l_hit.setEnergy(neutrinoleakage);
  l_hit = aux_infoHits->create();
  l_hit.setCellID(8);
  l_hit.setEnergy(leakage);
  G4AutoLock lock(&B4aEventActionMutex);
  //print here if you need event by event some information of the screen
  std::cout<<"--->Event number of activated fibers: "<<Fibers.size()<<"<---"<<std::endl;

    if (l_writer != NULL) l_writer->writeEvent();
    if (l_store != NULL) l_store->clearCollections();
  
}

void B4aEventAction::PrepareForRun()
{
  B4PodioManager* podioManager = B4PodioManager::Instance();
  
  G4AutoLock lock(&B4aEventActionMutex);
  
  podio::EventStore * l_evtstore = podioManager->GetEvtStore();
  if (l_evtstore == NULL) return;
  podio::ROOTWriter * l_writer = podioManager->GetWriter();

  s_caloHits = new edm4hep::SimCalorimeterHitCollection();
  l_evtstore->registerCollection("S_caloHits",s_caloHits);
  l_writer->registerForWrite("S_caloHits");

  c_caloHits = new edm4hep::SimCalorimeterHitCollection();
  l_evtstore->registerCollection("C_caloHits",c_caloHits);
  l_writer->registerForWrite("C_caloHits");

  aux_infoHits = new edm4hep::SimCalorimeterHitCollection();
  l_evtstore->registerCollection("Auxiliary_infoHits",aux_infoHits);
  l_writer->registerForWrite("Auxiliary_infoHits");

  s_caloHitContributions = new edm4hep::CaloHitContributionCollection();
  l_evtstore->registerCollection("S_caloHitContrib",s_caloHitContributions);
  l_writer->registerForWrite("S_caloHitContrib");

  c_caloHitContributions = new edm4hep::CaloHitContributionCollection();
  l_evtstore->registerCollection("C_caloHitContrib",c_caloHitContributions);
  l_writer->registerForWrite("C_caloHitContrib");
  
}

