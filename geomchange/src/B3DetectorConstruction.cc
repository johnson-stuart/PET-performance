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
/// \file B3DetectorConstruction.cc
/// \brief Implementation of the B3DetectorConstruction class

#include "B3DetectorConstruction.hh"
#include "B3PrimaryGeneratorAction.hh"
#include "G4Region.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSDoseDeposit.hh"
#include "G4VisAttributes.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4GeometryManager.hh"
#include "Randomize.hh"
#include "G4RegionStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"
#include "G4ios.hh"
#include "G4GenericMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3DetectorConstruction::B3DetectorConstruction(B3PrimaryGeneratorAction* p)
: G4VUserDetectorConstruction(),
  fMessenger(0),
  fPhysWorld(0), 
  fLogicWorld(0), 
  fLogicCrys(0),
  fLogicRing(0), 
  fLogicDet(0), 
  fLogicCyl(0),
  fCryst(0),
  fTank(0),
  fPrimaryGenerator(p) //order of these must follow the order listed in B3DetectorConstruction.hh
{
  DefineMaterials();
  fCylRadius = 10*mm;
  DefineCommands();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3DetectorConstruction::~B3DetectorConstruction()
{ 
	delete fMessenger;
	
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3DetectorConstruction::DefineMaterials()
{
  
  //Build LYSO
  G4double a; // mass of a mole
  G4double z; // mean number of protons
  G4String name, symbol;

  a=16.00*g/mole;
  G4Element*  elO = new G4Element(name="Oxygen",
				  symbol="O",
				  z=8., a);
  a=28.09*g/mole;
  G4Element*  elSi = new G4Element(name="Silicon",
				   symbol="Si",
				   z=14., a);
  a=174.97*g/mole;
  G4Element* elLu = new G4Element(name="Lutetium",
				  symbol="Lu",
				  z=71., a); 
  a=89.91*g/mole;
  G4Element* elY = new G4Element(name="Yttrium",
				 symbol="Y",
				 z=39., a);

  
  
G4Material* LYSO;
G4double density = 7.4*g/cm3;

 LYSO = new G4Material("Lu2Y2SiO5",
			density,
            4);

  LYSO->AddElement(elSi, 7*perCent);
  LYSO->AddElement(elLu, 71*perCent);
  LYSO->AddElement(elY, 4*perCent);
  LYSO->AddElement(elO, 18*perCent);


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B3DetectorConstruction::Construct()
{  

   //Cleanup old geometry
   
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  
  //Since I am changing the phantom dimensions, I want the source to remain at the same 
  //relative position within my phantom. I am required to make use of a function within 
  // my PrimarGeneratorAction class.
  G4double source_radius = 0.4*fCylRadius;
  G4double x0 = G4RandFlat::shoot((-source_radius-1)*mm, (-source_radius+1)*mm); 
  G4double y0 = G4RandFlat::shoot((source_radius-1)*mm, (source_radius+1)*mm);
  G4double z0 = G4RandFlat::shoot(-3*mm, 3*mm); 
  G4ThreeVector source_position = G4ThreeVector(x0,y0,z0);
  fPrimaryGenerator->SetPosition(source_position);


  //Defining the Detector Ring parameters
  //------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------
  
  G4double cryst_dX =40*mm, cryst_dY = 40*mm, cryst_dZ = 10*mm;
  G4int nb_cryst = 12;
  G4int nb_rings = 2;
  //
  G4double dPhi = twopi/nb_cryst, half_dPhi = 0.5*dPhi;
  G4double cosdPhi = std::cos(half_dPhi);
  G4double tandPhi = std::tan(half_dPhi);
  // 
  G4double ring_R1 = 0.5*cryst_dY/tandPhi;
  G4double ring_R2 = (ring_R1+cryst_dZ)/cosdPhi;
  G4double detector_dZ = nb_rings*cryst_dX; 
  //
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* default_mat = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* cryst_mat   = nist->FindOrBuildMaterial("Lu2Y2SiO5");
   
  //------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------
     
  // World
  //------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------
  
  G4double world_sizeXY = 0.5*ring_R1*cm;
  G4double world_sizeZ  = 0.5*ring_R1*cm;
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ); //its size
      
  fLogicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        default_mat,         //its material
                        "World");            //its name
                                
   fPhysWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      fLogicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      fCheckOverlaps);       // checking overlaps 
                 
				 
  //------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------

  // Define the Ring to house the crystals
  //------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------
  
  G4Tubs* solidRing =
    new G4Tubs("Ring", ring_R1, ring_R2, 0.5*cryst_dX, 0., twopi);
      
  fLogicRing =                         
    new G4LogicalVolume(solidRing,           //its solid
                        default_mat,         //its material
                        "Ring");             //its name
						
						
  //------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------   
 
 
  // Define crystal
  //------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------
  
  G4double gap = 0.5*mm;        //a gap for wrapping
  G4double dX = cryst_dX - gap, dY = cryst_dY - gap;
  G4Box* solidCryst = new G4Box("crystal", dX/2, dY/2, cryst_dZ/2);
                     
  fLogicCrys = 
    new G4LogicalVolume(solidCryst,          //its solid
                        cryst_mat,           //its material
                        "CrystalLV");        //its name
						
  
  for (G4int icrys = 0; icrys < nb_cryst ; icrys++) {
    G4double phi = icrys*dPhi;
    G4RotationMatrix rotm  = G4RotationMatrix();
    rotm.rotateY(90*deg); 
    rotm.rotateZ(phi);
    G4ThreeVector uz = G4ThreeVector(std::cos(phi),  std::sin(phi),0.);     
    G4ThreeVector position = (ring_R1+0.5*cryst_dZ)*uz;
    G4Transform3D transform = G4Transform3D(rotm,position);
                                    
    new G4PVPlacement(transform,             //rotation,position
                      fLogicCrys,            //its logical volume
                      "crystal",             //its name
                      fLogicRing,             //its mother  volume
                      false,                 //no boolean operation
                      icrys,                 //copy number
                      false);       // checking overlaps 
  }
  //------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------   
  
  
  // Full detector
  //------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------
  
  G4Tubs* solidDetector =
    new G4Tubs("Detector", ring_R1, ring_R2, 0.5*detector_dZ, 0., twopi);
      
  fLogicDet =                         
    new G4LogicalVolume(solidDetector,       //its solid
                        default_mat,         //its material
                        "Detector");         //its name
                                 
  G4double OG = -0.5*(detector_dZ + cryst_dX);
  for (G4int iring = 0; iring < nb_rings ; iring++) {
    OG += cryst_dX;
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(0,0,OG), //position
                      fLogicRing,             //its logical volume
                      "ring",                //its name
                      fLogicDet,         //its mother  volume
                      false,                 //no boolean operation
                      iring,                 //copy number
                      fCheckOverlaps);       // checking overlaps 
  }
                       
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    fLogicDet,           //its logical volume
                    "Detector",              //its name
                    fLogicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    fCheckOverlaps);         // checking overlaps 
                 
  //
  // Define phantom and containers
  //------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------
  
  // Phantom Container
  //------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------
  //Need to reconstruct the logical volume each time it runs with fLogicCyl.
  
  G4double tank_radius_inner = 0*mm;
  G4double tank_radius_outer = fCylRadius;
  G4double tank_length = 55*mm;
  G4RotationMatrix* tank_rotate = new G4RotationMatrix();
  
  G4Material* tank_mat = nist->FindOrBuildMaterial("G4_WATER");
  G4Tubs* solidTank = new G4Tubs("Tank", tank_radius_inner, tank_radius_outer, 0.5*tank_length, 0., 2*M_PI);
  fLogicCyl = new G4LogicalVolume(solidTank, tank_mat, "tankLV");
  new G4PVPlacement(tank_rotate, G4ThreeVector(), fLogicCyl, "Tank", fLogicWorld, false, 0, fCheckOverlaps);
  
  //Fish/Patient
  
  //will define the fish/patient by radiation distribution since it's homogenous density relative to water.
  
  // Visualization attributes
  //------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------
  fLogicRing->SetVisAttributes (G4VisAttributes::GetInvisible());
  fLogicDet->SetVisAttributes (G4VisAttributes::GetInvisible());    
  auto WaterVis= new G4VisAttributes(G4Colour(0.,0.,1.0,0.5)); //make the cylinder blue
  WaterVis->SetVisibility(true);
  WaterVis->SetForceSolid(true);
  fLogicCyl->SetVisAttributes(WaterVis);

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl; 

  //PHYSICAL WORLD RETURN
  //------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------
  return fPhysWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3DetectorConstruction::ConstructSDandField()
{
  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);
 
  // declare crystal as a MultiFunctionalDetector scorer
  //  
  if(!fCryst){
	fCryst = new G4MultiFunctionalDetector("crystal");
	G4SDManager::GetSDMpointer()->AddNewDetector(fCryst);
	G4VPrimitiveScorer* primitiv1 = new G4PSEnergyDeposit("edep");
	fCryst->RegisterPrimitive(primitiv1);
	SetSensitiveDetector("CrystalLV",fCryst);
  }
  // declare cylinder as the same
  if(!fTank){
	fTank = new G4MultiFunctionalDetector("tank");
	G4SDManager::GetSDMpointer()->AddNewDetector(fTank);
	G4VPrimitiveScorer* primitiv2 = new G4PSEnergyDeposit("edep_pat");
	fTank->RegisterPrimitive(primitiv2);
	SetSensitiveDetector("tankLV", fTank);
  }
}

void B3DetectorConstruction::SetCylinderRadius(G4double val)
{
  fCylRadius = val;
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}  

void B3DetectorConstruction::DefineCommands()
{
  // Define command directory using generic messenger class
  fMessenger = new G4GenericMessenger(this, 
                                      "/geomchange/update/", 
                                      "Detector control");

  // cylinder radius command
  auto& cylRadiusCmd
    = fMessenger->DeclareMethodWithUnit("cylR","mm",
                                &B3DetectorConstruction::SetCylinderRadius, 
                                "Set cylinder radius");
  cylRadiusCmd.SetParameterName("length", true);
  cylRadiusCmd.SetRange("length>0.");
  cylRadiusCmd.SetDefaultValue("10.");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
