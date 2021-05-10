//
// ********************************************************************
// * license and disclaimer                                           *
// *                                                                  *
// * the  geant4 software  is  copyright of the copyright holders  of *
// * the geant4 collaboration.  it is provided  under  the terms  and *
// * conditions of the geant4 software license,  included in the file *
// * license and available at  http://cern.ch/geant4/license .  these *
// * include a list of copyright holders.                             *
// *                                                                  *
// * neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  please see the license in the file  license  and url above *
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
// $Id: B4aEventAction.hh 75215 2013-10-29 16:07:06Z gcosmo $
// 
/// \file B4aEventAction.hh
/// \brief Definition of the B4aEventAction class

#ifndef B4aEventAction_h
#define B4aEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include <vector>
#include "G4ThreeVector.hh"

/// podio includes

#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/CaloHitContributionCollection.h"


/// Event action class 

class B4aEventAction : public G4UserEventAction
{
  public:
    B4aEventAction();
    virtual ~B4aEventAction();

    virtual void  BeginOfEventAction(const G4Event* event);
    virtual void  EndOfEventAction(const G4Event* event);
    
    //to save double
    void Addneutrinoleakage(G4double de);          //add energy of neutrinos in the leakage counter containing the calorimeter
    void Addleakage(G4double de);                  //add energy of all particles that are not neutrinos (or anti_neutrinos) in the leakage counter
    void Addem(G4double de);                       //add em component
    void AddScin(G4double de);                     //add energy in scintillating fibers
    void AddCher(G4double de);                     //add energy in Cherenkov fibers
    void AddCherenkov();                           //add cherenkov photoelectron
    void Addenergy(G4double de);                   //add energy deposited in the calorimeter
    void SavePrimaryParticle(G4int primaryparticlePDGID);
    void SavePrimaryEnergy(G4double primaryparticleenergy);

    //to save vectors in ntuple
    std::vector<G4double>& GetVectorSignalsR() {return VectorSignalsR;}
    std::vector<G4double>& GetVectorSignalsL() {return VectorSignalsL;} 
    std::vector<G4double>& GetVectorSignalsCherR() {return VectorSignalsCherR;}
    std::vector<G4double>& GetVectorSignalsCherL() {return VectorSignalsCherL;}
    std::vector<G4double>& GetVectorR() {return VectorR;}
    std::vector<G4double>& GetVectorL() {return VectorL;}
    std::vector<G4double>& GetVectorR_loop() {return VectorR_loop;}
    std::vector<G4double>& GetVectorL_loop() {return VectorL_loop;}

    //to fill vectors
    void AddVectorScinEnergyR(G4double de, G4int tower, G4int slice); //fill vector of scintillating fibers with signal (tower by tower)
    void AddVectorScinEnergyL(G4double de, G4int tower, G4int slice); //fill vector left side
    void AddVectorCherPER(G4int c_signal, G4int tower, G4int slice);  //fill vector of cherenkov fibers with chernekov signal (tower by tower) 
    void AddVectorCherPEL(G4int c_signal, G4int tower, G4int slice);  //fill vector left side
    void AddVectorR(G4double de, G4int tower, G4int slice);           //fill vector of energy deposited in tower (tower by tower)
    void AddVectorL(G4double de, G4int tower, G4int slice);           //fill vector left side
    void AddVectorR_loop(G4double de, G4int tower, G4int slice);
    void AddVectorL_loop(G4double de, G4int tower, G4int slice);

    //define fiber struct: contain info of fiber integrated signal and fiber location, include methods to fill it hit-by-hit
    struct Fiber{
      int ID, Type, Slice, Tower;           //fiber ID, type (S or C), slice, tower
      int E;                                //number of p.e. in fiber 
      G4ThreeVector Pos;                    //inner tip position (X,Y,Z) (mm)
      std::vector<double> phtimes;          //vector of float, each float is a single photon time of arrival 
      
      void addfiber(Fiber f1){             //function to add single hit to fiber integrated signal (single hit passed by reference)
	E = E+f1.E;                        //add hit photons
	phtimes.insert(phtimes.end(), f1.phtimes.begin(), f1.phtimes.end()); //append photons time of arrival
      };	
      
      void orderphtimes(){ //function to order photon times from first to last
	sort(phtimes.begin(), phtimes.end());
      };
    };

    void appendfiber(int ID, int Type, int Slice, int Tower, int E, G4ThreeVector Pos, std::vector<double> phtimes){//function to search if fiber already exists
      Fiber f{ID, Type, Slice, Tower, E, Pos, phtimes};
      auto it = find(FiberIDs.begin(), FiberIDs.end(), f.ID); //return iterator to fiber ID if exhists or FiberIDs.end() if fiber ID is not found
      if ( it == FiberIDs.end() ){       //fiber not found
        FiberIDs.push_back(f.ID);       //push_back fiber ID
	Fibers.push_back(f);            //push_back fiber
      }										
      else {                             //fiber found
        Fibers[distance(FiberIDs.begin(), it)].addfiber(f); //add fiber contribution with addfiber()
      }								
    };

    typedef struct FiberInfo {
        G4double F_ID, F_E, F_X, F_Y, F_Z; //fiber saturated energy
        G4int F_Type, F_slice, F_tower; //C==0 S==1;
    } Fiber_Info;
    
    void WriteFiber_Info(G4double FID, G4double FE, G4int FType, G4ThreeVector Fpos, G4int slice, G4int tower);
    
    void PrepareForRun();
    
private:
    std::vector<Fiber> Fibers;           //vector of Fibers
    std::vector<int> FiberIDs;           //vector of Fibers IDs
    Fiber_Info Fiber_Hits[1000000];
    G4double  Energyem;             //energy of em component
    G4double  EnergyScin;           //energy in scintillating fibers
    G4double  EnergyCher;           //energy in Cherenkov fibers
    G4int     NofCherenkovDetected; //number of Cherenkov photons detected (in cherenkov fibers)
    G4double  EnergyTot;            //total energy deposited (does not count invisibile energy)
    G4int PrimaryParticleID;        //PDGID of primary particle
    G4double PrimaryParticleEnergy; //primary particle energy
    G4double neutrinoleakage;       //leakage neutrino
    G4double leakage;               //leakage non neutrino

    std::vector<G4double> VectorR_loop;
    std::vector<G4double> VectorL_loop;
    std::vector<G4double> VectorSignalsR;     //vector filled with scintillating signal tower-by-tower
    std::vector<G4double> VectorSignalsL;     //vector filled for right side
    std::vector<G4double> VectorSignalsCherR; //Vector filled with Cherenkov signal tower-by-tower
    std::vector<G4double> VectorSignalsCherL; //vector filled for left side
    std::vector<G4double> VectorR;            //vector with energy deposited in towers 
    std::vector<G4double> VectorL;            //vector fille for left side
    
    //define edm4hep calo hits + auxiliary info
    edm4hep::SimCalorimeterHitCollection * s_caloHits;
    edm4hep::SimCalorimeterHitCollection * c_caloHits;
    edm4hep::SimCalorimeterHitCollection * aux_infoHits;

  // detailed information for timing of the signal

  edm4hep::CaloHitContributionCollection * s_caloHitContributions;
  edm4hep::CaloHitContributionCollection * c_caloHitContributions;
};

// inline functions
inline void B4aEventAction::Addneutrinoleakage(G4double de){neutrinoleakage += de;}

inline void B4aEventAction::Addleakage(G4double de){leakage += de;}

inline void B4aEventAction::AddVectorR(G4double de, G4int tower, G4int slice){VectorR.at(tower+(slice*75)) += de;}

inline void B4aEventAction::AddVectorL(G4double de, G4int tower, G4int slice){
	tower = -1*tower;
	VectorL.at(tower+(slice*75)) += de;
}

inline void B4aEventAction::AddVectorR_loop(G4double de, G4int tower, G4int slice){
    VectorR_loop.at(tower+(slice*75)) = de; 
}

inline void B4aEventAction::AddVectorL_loop(G4double de, G4int tower, G4int slice){
    tower = -1*tower;
    VectorL_loop.at(tower+(slice*75)) = de;
}

inline void B4aEventAction::WriteFiber_Info(G4double FID, G4double FE, G4int FType, G4ThreeVector Fpos, G4int slice, G4int tower){
    int k=0;
    while (Fiber_Hits[k].F_ID!=0 && Fiber_Hits[k].F_ID!=FID){
        k++;}
    Fiber_Hits[k].F_ID = FID;
    Fiber_Hits[k].F_E += FE;
    Fiber_Hits[k].F_Type = FType;
    Fiber_Hits[k].F_X = Fpos[0];
    Fiber_Hits[k].F_Y = Fpos[1];
    Fiber_Hits[k].F_Z = Fpos[2];
    Fiber_Hits[k].F_slice = slice;
    Fiber_Hits[k].F_tower = tower;
}

inline void B4aEventAction::SavePrimaryParticle(G4int primaryparticlePDGID){PrimaryParticleID = primaryparticlePDGID;}

inline void B4aEventAction::SavePrimaryEnergy(G4double primaryparticleenergy){PrimaryParticleEnergy = primaryparticleenergy;}

inline void B4aEventAction::AddVectorScinEnergyR(G4double de, G4int tower, G4int slice) {
    VectorSignalsR.at(tower+(slice*75)) += de;
}

inline void B4aEventAction::AddVectorScinEnergyL(G4double de, G4int tower, G4int slice) {
    tower = -1*tower;
    VectorSignalsL.at(tower+(slice*75)) += de;
}

inline void B4aEventAction::AddVectorCherPEL(G4int c_signal, G4int tower, G4int slice) {
	tower = -1*tower;
    VectorSignalsCherL.at(tower+(slice*75)) = VectorSignalsCherL.at(tower+(slice*75))+c_signal;
}

inline void B4aEventAction::AddVectorCherPER(G4int c_signal, G4int tower, G4int slice) {
    VectorSignalsCherR.at(tower+(slice*75)) = VectorSignalsCherR.at(tower+(slice*75))+c_signal;
}

inline void B4aEventAction::Addem(G4double de) {Energyem += de;}

inline void B4aEventAction::AddScin(G4double de){EnergyScin += de;}

inline void B4aEventAction::AddCher(G4double de){EnergyCher += de;}

inline void B4aEventAction::Addenergy(G4double de){EnergyTot += de;}
                     
#endif

    
