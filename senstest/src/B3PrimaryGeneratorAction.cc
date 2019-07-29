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
/// \file B3PrimaryGeneratorAction.cc
/// \brief Implementation of the B3PrimaryGeneratorAction class

#include "B3PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ChargedGeantino.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include "G4GenericMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3PrimaryGeneratorAction::B3PrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(),
   fParticleGun1(0),
   fParticleGun2(0),
   fMessenger(0)
{
  DefineCommands();
  G4int n_particle = 1;
  fParticleGun1  = new G4ParticleGun(n_particle);
  fParticleGun2 = new G4ParticleGun(n_particle);
  
  //Set a default axial offset and z offset.
  fAx_position = G4ThreeVector();
  fZ_position = 0*mm;
  
  
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;

  //----------------------------------
  // Photon 1
  fParticleGun1->SetParticleDefinition(particleTable->FindParticle(particleName="gamma"));
  fParticleGun1->SetParticleEnergy(511*keV);
  fParticleGun1->SetParticlePosition(fAx_position);
  
  //----------------------------------
  // Photon 2
  fParticleGun2->SetParticleDefinition(particleTable->FindParticle(particleName="gamma"));
  fParticleGun2->SetParticleEnergy(511*keV);
  fParticleGun2->SetParticlePosition(fAx_position);
       
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3PrimaryGeneratorAction::~B3PrimaryGeneratorAction()
{
  delete fParticleGun1;
  delete fParticleGun2;
  delete fMessenger;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  
   //set random beam axis
  G4ThreeVector beam_axis;
  G4double theta    = 999999.9;
  G4double cosTheta = 99999.9;
  G4double sinTheta = 9999.9;
  G4double phi      = 999.9;
  G4double thetaMax = 90.; // determines the arc that you want to generate the photons

 
	
  // thetaMax degrees is just outside of 
  // outer crystal outer edge
  while(theta > thetaMax ){
   cosTheta  = 2.*G4UniformRand()-1.;
   G4double sinTheta2 = 1. - cosTheta*cosTheta;
   if( sinTheta2 < 0.)  sinTheta2 = 0.;
    sinTheta  = std::sqrt(sinTheta2); 
    phi       = 2*M_PI*G4UniformRand();
    theta = std::acos(cosTheta)* 180/(M_PI);
   }
      
  // G4cout << " theta = " << theta << G4endl;
      
  // theta wrt x-axis (fixed beam in x)
  beam_axis.set(cosTheta,sinTheta*std::cos(phi), sinTheta*std::sin(phi));
  beam_axis = beam_axis.unit(); 
  
  
  //you want the particles direction and position Sto be randomised each event
  //so the beam_axis calc goes in generateprimaries which generates
  //each event.
  
  fParticleGun1->SetParticleMomentumDirection(beam_axis);
  fParticleGun2->SetParticleMomentumDirection(-beam_axis);
  
  fParticleGun1->GeneratePrimaryVertex(anEvent);
  fParticleGun2->GeneratePrimaryVertex(anEvent);

}

void B3PrimaryGeneratorAction::SetAxialPosition(G4double axial_offset)
{
	//turn axial_offset from double to vector 
	G4double x_value = axial_offset/sqrt(2);
	G4double y_value = axial_offset/sqrt(2);
	G4double z_value = 0;
	G4ThreeVector particle_position = G4ThreeVector(x_value, y_value, z_value);
	fParticleGun1->SetParticlePosition(particle_position);
	fParticleGun2->SetParticlePosition(particle_position);
	fAx_position = particle_position;
}

void B3PrimaryGeneratorAction::SetZPosition(G4double z_offset)
{
	//turn z_offset from double to vector 
	G4double x_value = 0;
	G4double y_value = 0;
	G4double z_value = z_offset;
	G4ThreeVector particle_position = G4ThreeVector(x_value, y_value, z_value);
	fParticleGun1->SetParticlePosition(particle_position);
	fParticleGun2->SetParticlePosition(particle_position);
	fAx_position = particle_position;
	fZ_position = z_offset;
}

void B3PrimaryGeneratorAction::DefineCommands()
{
  // Define command directory using generic messenger class
  fMessenger = new G4GenericMessenger(this, 
                                      "/geomchange/update/", 
                                      "Detector control");

  // cylinder radius command
  auto& axlRadiusCmd
    = fMessenger->DeclareMethodWithUnit("axlR","mm",
                                &B3PrimaryGeneratorAction::SetAxialPosition, 
                                "Set axial offset of source");
  axlRadiusCmd.SetParameterName("length", true);
  //axlRadiusCmd.SetRange("length>=0.");
  axlRadiusCmd.SetDefaultValue("0.");
  
  auto& zCmd
    = fMessenger->DeclareMethodWithUnit("zOff","mm",
                                &B3PrimaryGeneratorAction::SetZPosition, 
                                "Set axial offset of source");
  zCmd.SetParameterName("length", true);
  //zCmd.SetRange("length>=0.");
  zCmd.SetDefaultValue("0.");
  
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

