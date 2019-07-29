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
/// \file B3bRunAction.cc
/// \brief Implementation of the B3bRunAction class

#include "B3bRunAction.hh"
#include "B3PrimaryGeneratorAction.hh"
#include "B3DetectorConstruction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "B3bAnalysis.hh"
#include "G4AccumulableManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3bRunAction::B3bRunAction()
 : G4UserRunAction(),
   fGoodEvents(0),
   fScatterEvents(0)

{  
  
  //set up analysis manager
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetNtupleMerging(true);

  //tuple

  analysisManager->CreateNtuple("E-deposited Crystal", "Edep");
  analysisManager->CreateNtupleDColumn("Eabs_cryst");
  analysisManager->FinishNtuple();
  analysisManager->CreateNtuple("E-deposited Apparatus", "Edep_pat");
  analysisManager->CreateNtupleDColumn("Eabs_pat");
  analysisManager->FinishNtuple();
  
  //make histogram to store patient radius with scatter fracton
  //format: CreateH1(name, title, nbins, from, to)
  
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->RegisterAccumulable(fGoodEvents);
  accumulableManager->RegisterAccumulable(fScatterEvents);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3bRunAction::~B3bRunAction()
{ 
 delete G4AnalysisManager::Instance(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3bRunAction::BeginOfRunAction(const G4Run* run)
{ 

  G4cout << "### Run " << run->GetRunID() << " start." << G4endl;
  
  
  // reset accumulables to their initial values
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Reset();
  
  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Open an output file
  //
  G4String fileName = "B3b";
  analysisManager->OpenFile(fileName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3bRunAction::EndOfRunAction(const G4Run* run)
{
  G4double nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;
  
  
  // Merge accumulables 
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Merge();
  
  //details about the primary generator (source)
  const B3PrimaryGeneratorAction* generatorAction
    = static_cast<const B3PrimaryGeneratorAction*>(
        G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String partName;
  if (generatorAction) 
  {
    G4ParticleDefinition* particle = generatorAction->GetParticleGun()->GetParticleDefinition();
    partName = particle->GetParticleName();
  }  
  
  //calculate the radius of source offset for output
  G4ThreeVector axial_position = generatorAction->GetPosition();
  
  G4double axial_x = axial_position.getX();
  G4double axial_y = axial_position.getY();
  G4double axial_z = axial_position.getZ();
  
  G4double axial_r = sqrt(pow(axial_x, 2) + pow(axial_y, 2) );
  
  
  //results
  //sensitivity will be the true counts over the initiated events - ratio of true 
  //counts to decays essentiall
  G4double sensitivity = fGoodEvents.GetValue()/nofEvents;
  //print number of photons absorbed in the patient and crystals

  //
  if (IsMaster())
  {
    G4cout
     << G4endl
     << "--------------------End of Global Run-----------------------"
     << G4endl
     << "  The run was " << nofEvents << " events " << G4endl;
  }
  else
  {
    G4cout
     << G4endl
     << "--------------------End of Local Run------------------------"
     << G4endl
     << "  The run was " << nofEvents << " "<< partName << G4endl;
  }      

  G4cout
     << "; Axial Offset: " << axial_r << G4endl
	 << "; Z Offset: " << axial_z << G4endl
     << "; Good Counts: " << fGoodEvents.GetValue()  << G4endl
	 << "; Sensitivity: " << sensitivity << G4endl
     << G4endl;
	 
  //writing specific results to output file
  std::ofstream outfile;
  outfile.open("output.log", std::ios::app);
  outfile << sensitivity << " " << axial_z << std::endl;
  outfile.close();
  auto analysisManager = G4AnalysisManager::Instance();
  // save histograms & ntuple
  //
  analysisManager->Write();
  analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
