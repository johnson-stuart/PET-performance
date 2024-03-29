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
   fP1_position(0)
{
  G4int n_particle = 1;
  fParticleGun1  = new G4ParticleGun(n_particle);
  fParticleGun2 = new G4ParticleGun(n_particle);
  
  //use of the RandFlat function to set up a random distribution of decay in a target area.
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;

  //----------------------------------
  // Photon 1
  fParticleGun1->SetParticleDefinition(particleTable->FindParticle(particleName="gamma"));
  fParticleGun1->SetParticleEnergy(511*keV);

  //----------------------------------
  // Photon 2
  fParticleGun2->SetParticleDefinition(particleTable->FindParticle(particleName="gamma"));
  fParticleGun2->SetParticleEnergy(511*keV);
       
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3PrimaryGeneratorAction::~B3PrimaryGeneratorAction()
{
  delete fParticleGun1;
  delete fParticleGun2;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
   
  //set random particle position based on the geometry! :)
  G4double p1_x0 = fP1_position.getX();
  G4double p1_y0 = fP1_position.getY();
  G4double p1_z0 = fP1_position.getZ();
   
  G4double x0 = G4RandFlat::shoot((-p1_x0-1)*mm, (-p1_x0+1)*mm); 
  G4double y0 = G4RandFlat::shoot((p1_y0-1)*mm, (p1_y0+1)*mm);
  G4double z0 = G4RandFlat::shoot((p1_z0-3)*mm, (p1_z0+3)*mm); 
  
  G4ThreeVector new_p1_position = G4ThreeVector(x0,y0,z0);
  
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
  
  fParticleGun1->SetParticlePosition(new_p1_position);
  fParticleGun2->SetParticlePosition(new_p1_position);
  
  fParticleGun1->GeneratePrimaryVertex(anEvent);
  fParticleGun2->GeneratePrimaryVertex(anEvent);

}

void B3PrimaryGeneratorAction::SetPosition(G4ThreeVector particle_position)
{
	fParticleGun1->SetParticlePosition(particle_position);
	fParticleGun2->SetParticlePosition(particle_position);
	fP1_position = particle_position;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

