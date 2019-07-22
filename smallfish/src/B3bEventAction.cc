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
//
/// \file B3bEventAction.cc
/// \brief Implementation of the B3bEventAction class

#include "B3bEventAction.hh"
#include "B3bRunAction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"

#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4THitsMap.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "B3bAnalysis.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3bEventAction::B3bEventAction(B3bRunAction* runAction)
 : G4UserEventAction(), 
   fRunAction(runAction),
   fCollID_cryst(-1),
   fCollID_patient(-1),
   fPrintModulo(10000)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3bEventAction::~B3bEventAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3bEventAction::BeginOfEventAction(const G4Event* /*evt*/)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3bEventAction::EndOfEventAction(const G4Event* evt )
{
  G4int evtNb = evt->GetEventID();
  
  if (evtNb%fPrintModulo == 0) { 
    G4cout << G4endl << "---> end of event: " << evtNb << G4endl;
  }      
   
   auto analysisManager = G4AnalysisManager::Instance();
   //Hits collections
  //  
  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  if(!HCE) return;
               
   // Get hits collections IDs
  if (fCollID_cryst < 0) {
    G4SDManager* SDMan = G4SDManager::GetSDMpointer();  
    fCollID_cryst   = SDMan->GetCollectionID("crystal/edep");
    fCollID_patient = SDMan->GetCollectionID("tank/edep_pat");    
  }
  
  std::map<G4int,G4double*>::iterator itr;
  
  //Energy deposit in patient
  //
     
  G4THitsMap<G4double>* evtMapP = (G4THitsMap<G4double>*)(HCE->GetHC(fCollID_patient));
               
  G4double pat = 0.0;
  
  for (itr = evtMapP->GetMap()->begin(); itr != evtMapP->GetMap()->end(); itr++) {
    G4double edep_pat = *(itr->second);
	if(edep_pat > 0.0){
		pat = edep_pat;
	}
	analysisManager->FillNtupleDColumn(1,0, edep_pat);
    analysisManager->AddNtupleRow(1);
}  
  

  //Energy in crystals : identify 'good events'
  //
  const G4double eGood = 511*keV;
  G4int nbOfFired = 0;
  G4int nbGood = 0;
  
   
  G4THitsMap<G4double>* evtMap = 
                     (G4THitsMap<G4double>*)(HCE->GetHC(fCollID_cryst));
               
  for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++) {
    G4double edep = *(itr->second);
	analysisManager->FillNtupleDColumn(0,0, edep);
	analysisManager->AddNtupleRow(0);
    if (edep >= 400*keV) {
		nbOfFired++;
	}
	if (edep == eGood){
		nbGood++;
	}
  }  
    
  if (nbOfFired == 2) {
	  if(nbGood == 2){
		fRunAction->CountEvent();
	  }
	  if(pat > 0.0){
	    fRunAction->CountScatter();
	  }
  }
  
  
 
  
  
  
  
		 
  }
    
  
  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

