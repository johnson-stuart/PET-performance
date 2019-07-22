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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3PrimaryGeneratorAction::B3PrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(),
   fParticleGun(0)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3PrimaryGeneratorAction::~B3PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4double x0 = G4RandFlat::shoot(-2.5*mm, 2.5*mm); 
  //G4double x0 = 0*mm;
  G4double z0 = G4RandFlat::shoot(-3*mm, 3*mm); 
  G4double y0 = -17*mm + G4RandFlat::shoot(-2.5*mm, 2.5*mm);

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


  //----------------------------------
  // Photon 1
  fParticleGun->SetParticleDefinition(particleTable->FindParticle(particleName="gamma"));
  fParticleGun->SetParticleEnergy(511*keV);
  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  fParticleGun->SetParticleMomentumDirection(beam_axis);
    
  G4ThreeVector random = G4RandomDirection();
    
  // Restrict to YZ plane
  G4ThreeVector randomYZ = random.cross(beam_axis).unit();
    
  G4ThreeVector y_axis   = G4ThreeVector(0,1,0);
  //fParticleGun->SetParticlePolarization(randomYZ);
  fParticleGun->GeneratePrimaryVertex(anEvent);

  //----------------------------------
  // Photon 2
  fParticleGun->SetParticleEnergy(511*keV);
  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  fParticleGun->SetParticleMomentumDirection(-beam_axis);

  // randomise second photon polarisation wrt first
  G4ThreeVector randomYZ_2 = G4RandomDirection().cross(beam_axis).unit();
  // 
  G4ThreeVector z_axis = G4ThreeVector(0,0,1);
  // perpendicular to first photon
  //G4ThreeVector perpRandomYZ = beam_axis.cross(randomYZ).unit();
  //fParticleGun->SetParticlePolarization(randomYZ_2);

            
  //create vertex
  //
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

