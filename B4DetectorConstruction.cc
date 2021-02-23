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
/// \file B4DetectorConstruction.cc
/// \brief Implementation of the B4DetectorConstruction class

//************* This code creates a hexagonal cone formed of three layers (two scintillator and a middle carbon layer) that will be used as
//************* polarimeter for neutrons (and protons).

#include "B4DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "B4Analysis.hh"

#include "G4Box.hh"
#include "G4Polyhedra.hh"
#include "G4Polycone.hh"
#include "G4Trap.hh"
#include "G4Trd.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UserLimits.hh"
#include <sstream>
#include <iostream>


//.....ooo0000oooo VARIABLES TO BE CHANGED EXTERNALLY....ooo0000ooo....

extern G4double gplasticThickness;
//extern G4double gcarbonThickness;

// Other variables
G4double density, massfraction;
G4int natoms, nel;

//....oooOO0OOooo........oooOO0OOoooint........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal
G4GlobalMagFieldMessenger* B4DetectorConstruction::fMagFieldMessenger = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4DetectorConstruction::B4DetectorConstruction()
 : G4VUserDetectorConstruction(),
   fAbsorberPV(),
   //fGapPV(nullptr),
   fGapPV(),
   //fGapPV_2(nullptr),
   //fGapPV_3(nullptr),
   //fGapPV_4(nullptr),
   //fGapPV_5(nullptr),
   //fGapPV_6(nullptr),
   fAbsorber1PV(),
   fCheckOverlaps(true)
   //fStepLimit(NULL)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4DetectorConstruction::~B4DetectorConstruction()
{
  //delete fStepLimit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B4DetectorConstruction::Construct()
{
  // Define materials
  DefineMaterials();

  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo....DEFINING MATERIALS....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4DetectorConstruction::DefineMaterials()
{
  // Material defined using NIST Manager
  auto nistManager = G4NistManager::Instance();
  //Define a carbon material using nistManager

  nistManager->FindOrBuildMaterial("G4_C");
  //nistManager->FindOrBuildMaterial("G4_Pb");

 //Define a scintillator (plastic) material
  G4Material* Sci = new G4Material("Scintillator", density=1.032*g/cm3, nel=2);

  //Define carbon and hyrogen material for building scintillator material
  G4Element* elC = nistManager->FindOrBuildElement("C");
  G4Element* elH = nistManager->FindOrBuildElement("H");
  G4Element* elN  = nistManager->FindOrBuildElement("N");
  G4Element* elO  = nistManager->FindOrBuildElement("O");
  G4Element* elAr  = nistManager->FindOrBuildElement("Ar");

  Sci->AddElement(elC, natoms=9);
  Sci->AddElement(elH, natoms=10);

  // Define a liquid argon material (not in use)
  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;
  G4double density;
  new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
         // The argon by NIST Manager is a gas with a different density

  // Vacuum
  new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);

 // air
  //density = 1.2929e-03 *g/cm3;  // at 20 degree
  const G4double expTemp = STP_Temperature+20.*kelvin; //20 degrees, e.g. approx room temperature
  G4Material* Air = new G4Material("Air", density = 1.2929e-03 *g/cm3, nel=3,  kStateGas, expTemp);
  G4double ttt = 75.47+23.20+1.28;
  Air-> AddElement(elN,  massfraction= 75.47/ttt);
  Air-> AddElement(elO,  massfraction= 23.20/ttt);
  Air-> AddElement(elAr, massfraction=  1.28/ttt);


  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo...DEFINING AND PLACING LAYERS.....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B4DetectorConstruction::DefineVolumes()
{

  //G4double maxStep = .001*mm;
  //fStepLimit = new G4UserLimits(maxStep);  //can define a step size, but not in use currently

  // Geometry parameters
  auto worldSizeXY = 1.*m; //define worldsize
  auto worldSizeZ  = 1.*m;

  // Get materials
  auto defaultMaterial = G4Material::GetMaterial("Galactic");
  auto plasticMaterial = G4Material::GetMaterial("Scintillator");
  auto carbonMaterial = G4Material::GetMaterial("G4_C");

  if ( ! defaultMaterial || ! plasticMaterial || ! carbonMaterial ) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined.";
    G4Exception("B4DetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  }

  //
  // World
  //
  auto worldS
    = new G4Box("World",           // its name
                 worldSizeXY/2, worldSizeXY/2, worldSizeZ/2); // its size

  auto worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 defaultMaterial,  // its material
                 "World");         // its name

  auto worldPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 worldLV,          // its logical volume
                 "World",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps

///Defining Cone shape
            //Defining the variables for creating the Hexagonal Cone
                 G4double thick 			= 2*cm; //thickness of the carbon layer
                 G4double openingangle =45*deg;
         		     G4double lengthwithnose 		= 15*cm;  //distance from front to back horizontally of the cone
                 G4double gap 			= 2*mm; //thickness of the carbon layer
                 G4double length 		= lengthwithnose-thick; //length of the entire cone minus the length of the nose to give the length of the "main" part of the cone
                 G4double aperture     = 2*cm;   //diameter of the aperture of the carbon aka the width of the "cone side" at the small/back aperture
                 G4double bottomRad     = (sqrt(3)*thick/3)+aperture;
                 G4double topRad        = lengthwithnose*tan(openingangle);     //the width of the "cone side" at the large/front aperture

                 G4double x= ((0.5*sqrt(3)/3)*(topRad+(sqrt(3)*bottomRad))); //radius of cut out prism for carbon
                 G4double yangle = 180*deg;  //determines direction of cone (e.g. along z axis) -needs to be 180*deg so cone is facing correct direction (positive z)
                 G4double zangle = 0*deg;  //optional variable that can be defined differently  (zero - not in use)
                 G4double x_coords[6] = {0,0,0,0,0,0};  //optional variable that can be defined differently  (zero - not in use)
                 G4double y_coords[6] = {0,0,0,0,0,0}; //optional variable that can be defined differently  (zero - not in use)
                 G4double y_angles[6] = {0,0,0,0,0,0};  //optional variable that can be defined differently  (zero - not in use)
                 G4double z_angles[6] = {0,0,0,0,0,0};  //optional variable that can be defined differently  (zero - not in use)
                 G4RotationMatrix* rot;
                 G4RotationMatrix* rot2;
                 G4RotationMatrix* rot3;
                 rot = new G4RotationMatrix(0*deg,openingangle,0*deg);    //angle for creating the nose (e.g. diagonal half of a box) needs to be the same as tan(angle)=y/15cm where angle is the opening angle and y is the height of the cone (e.g. topRad)
                 rot2 = new G4RotationMatrix(30*deg,0*deg,0*deg);   //angle for cutting out the side of the cone
                 rot3 = new G4RotationMatrix(180*deg,270*deg,0*deg);   //angle to rotate nose by when adding it to the cone

        //Defining the simple shape for each side of the cone
            //Cut out shape dimensions
                G4double zplanesg[] = {-200*cm, 210*cm};     //define two planes for the subtraction solid that are an arbitary long distance
                G4double rinnerg[] = {0.*mm, 0*mm};          //define the inner radius as zero (aka solid shape) for the subtraction solid (a triangular prism)
                G4double routerg[] = {x, x};          //define outer radius  (distance to the middle of the sides) for triangular prism
            //Cone dimensions
                 G4double zPositions[] 	= {-1/2.*length, 1/2.*length};  //the front and back positions of the cone
         		     G4double rInner[] 		= {bottomRad-thick/2., topRad-thick/2.};   //the inner radius of the cone
         		     G4double rOuter[]		= {bottomRad + thick/2., topRad+thick/2.};  //the outer radius of the cone


         		     G4Polyhedra * pconebig = new G4Polyhedra("pcone", 0, 2*pi, 6, 2,
         									zPositions, rInner, rOuter);   //Define a cone that will then be split into pieces
                 G4Polyhedra* pconesubtract = new G4Polyhedra("CutOutShape",
                                        0.,
                                        2*pi,
                                        3,
                                        2,
                                        zplanesg,
                                        rinnerg,
                                        routerg)	;              //Define a triangular prism that can be used to "cut out" each side of the cone

//Creating the Carbon layer
      //Creating the "nose" of the cone (for the inner aperture) for the carbon layer- a triangular prism with a triangle with 45 degrees
G4double y= ((thick/1)/tan(openingangle))/2;  //
G4double z= (thick/1)/2;  //
G4double h= z/sin(openingangle);  //
G4double e = sin((90*deg)-openingangle)*h;
G4double d = (h*cos((90*deg)-openingangle));


                auto box1
                     = new G4Box("box1",     // its name
                                            bottomRad, y, z );        //define a box that is twice the size of the triangular prism needed for the cone
                                            //the nose needs to be as wide as the side on the inner apertures and as thick as the side (to get 45 degrees y=z)
                auto box2
                     = new G4Box("box2",     // its name
                                           bottomRad+100, h, h);     //define another box twice as wide and arbitarily longer to subtract from box1 to make the nose
//bottomRad+100, ((thick/cos(openingangle))/tan(openingangle))/(2*sin(45)), (thick/cos(openingangle))/(2*sin(45)));

                G4VSolid* pconenose =
                          new G4SubtractionSolid("solid_case", // name
                          box1,  // solid A
                          box2,  // solid B
                          rot,            // rotation
                          G4ThreeVector(0*cm,d,-e)); // subtract box2 from box1 to make the nose //G4ThreeVector(0*cm,tan(openingangle)*(thick),-tan((90*deg)-openingangle)*(thick/2))
//(sin(openingangle)*h*cos(openingangle)),-(cos(openingangle)*h*cos(openingangle))))
             G4VSolid* pconemain =
                            new G4UnionSolid("solid_case",
                            pconebig,
                            pconenose,
                            rot3,
                            G4ThreeVector(0*cm,(bottomRad), -((length + 2*y)/2)));   //add the nose onto the cone



              G4VSolid* pcone = new G4IntersectionSolid("solid_case", // name
                                       pconemain,  // solid A
                                       pconesubtract,  // solid B
                                       rot2,            // rotation
                                       G4ThreeVector(0,2*x, 0)); //use the large triangular prism to cut the side with the nose out of the cone
//2*(topRad+(1.16*thick))
//((2*sqrt(3)/3)*topRad)+bottomRad+thick






              auto pconeLV
                   = new G4LogicalVolume(
                                pcone,     // its solid
                                carbonMaterial,  // its material
                                "pcone");   // its name             //Create the logical volume of the carbon layer side

                                std::stringstream converter;
                                std::stringstream converter1;

                for (int i = 0; i < 6; i++) {
                                //G4double y_angle = y_angles[i];
                                //G4double z_angle = z_angles[i];
                                converter.str(" ");
                                converter << "" << i << "";
                                converter1.str(" ");
                                converter1 << "carbon_" << i << "";
                                rot = new G4RotationMatrix((30+(i*60))*deg,yangle,zangle);
                                G4double x_coord = x_coords[i];
                                G4double y_coord = y_coords[i];
                                fGapPV[i] = new G4PVPlacement(                                           //auto pconeside_physical[i] =
                                rot,       // rotation
                                G4ThreeVector(x_coord,y_coord,0*m),  // at Ex. G4ThreeVector((i-3)*m,(i-3)*m,2.5*m)
                                pconeLV,          // its logical volume
                                converter1.str(),    // its name
                                worldLV,          // its mother  volume
                                false,            // no boolean operation
                                0,                // copy number
                                fCheckOverlaps);  // checking overlaps
                              }                                       //Place the carbon sides in the world

                              // G4double length 		= lengthwithnose-thick; //length of the entire cone minus the length of the nose to give the length of the "main" part of the cone
                              // G4double aperture     = 2*cm;   //diameter of the aperture of the carbon aka the width of the "cone side" at the small/back aperture
                              // G4double bottomRad     = (sqrt(3)*thick/3)+aperture;
                              // G4double topRad        = lengthwithnose*tan(openingangle);     //the width of the "cone side" at the large/front aperture

//Scintillator layers
//Defining the values for the outer scintillator layer
G4double innerscintlength = lengthwithnose - 2*z;
            G4double scintthick = 0.5*cm;   //scintillator layer thickness, 0.5cm, 0.7cm, 1cm (will need to adjust other values to remove gaps for anything other than 0.5cm)
            //G4double diamouterscint =
G4double bottomscintRad1     = (sqrt(3)*thick/3)+aperture;     //the width of the "cone side" at the small/back aperture for inner layer
            G4double topscintRad1        = (innerscintlength*tan(openingangle));     //the width of the "cone side" at the large/front aperture for inner layer
            G4double outerscintlength = length+(2*z);
            G4double bottomscintRad     = (sqrt(3)*thick/3)+aperture;    //the width of the "cone side" at the small/back aperture for outer layer
            G4double topscintRad        = (outerscintlength*tan(openingangle));     //the width of the "cone side" at the large/front aperture for outer layer
            G4double x1        =((0.5*sqrt(3)/3)*(topscintRad+(sqrt(3)*bottomscintRad)));
            G4double zPositions1[] 	= {-1/2.*(length+(2*y)), 1/2.*(length+(2*y))};  //the front and back positions of the top layer of scint cone
            G4double rInnerscint[] 		= {bottomRad -thick/2 +gap , topRad + thick/2 +gap };   //the inner radius of the cone for outer layer
            G4double rOuterscint[]		= {bottomRad-thick/2 + scintthick +gap , topRad + thick/2+ scintthick +gap};  //outer radius of the cone for outer layer
            G4double rscintinnerg[] = {0.*mm, 0*mm};            //define as zero to get a solid shape for cutting out
            G4double rscintouterg[] = {x1, x1};       //define radius for cutting out shape
            G4double a = 0.*mm ;
            G4double b = 0.*mm  ;
            G4double a_coords[6] = {0,-a,0,a,a,a};
            G4double b_coords[6] = {b,0,-b,-b,0,b};
//Defining the values for the inner scintillator layer
// G4double rInnerscint[] 		= {bottomRad -thick/2 , topscintRad1 + (bottomRad - thick/2) };   //the inner radius of the cone for outer layer
// G4double rOuterscint[]		= {bottomRad-thick/2 + scintthick, topscintRad1 + (bottomRad - thick/2 + scintthick) };  //outer radius of the cone for outer layer

            G4double x2       =((0.5*sqrt(3)/3)*(topscintRad1+(sqrt(3)*bottomscintRad1)));
            //G4double bottomscintlength = length - (((scintthick+thick)/2)/cos(openingangle));
            G4double zPositions2[] 	= {-1/2.*(length- (2*scintthick/tan(openingangle))), 1/2.*(length- (2*scintthick/tan(openingangle)))};  //the front and back positions of the top layer of scint cone
            G4double rInnerscint1[] 		= {bottomRad-thick/2 -scintthick -gap + (2*scintthick), topRad - thick/2 -scintthick -gap};   //the inner radius of the cone for inner layer
            G4double rOuterscint1[]		= {bottomRad-thick/2 -gap + (2*scintthick), topRad - thick/2 -gap};    //outer radius of the cone for inner layer
            // G4double rInnerscint1[] 		= {bottomRad-thick/2 -scintthick, topscintRad1+ (bottomRad-thick/2 - scintthick)  };   //the inner radius of the cone for inner layer
            // G4double rOuterscint1[]		= {bottomRad-thick/2, topscintRad1+(bottomRad-thick/2) };    //outer radius of the cone for inner layer
            G4double rscintinnerg1[] = {0.*mm, 0*mm};           //define as zero to get a solid shape for cutting out
            G4double rscintouterg1[] = {x2, x2};       //define radius for cutting out shape (smaller for inner layer)
            // G4double rInner[] 		= {bottomRad-thick/2., topRad-thick/2.};   //the inner radius of the cone
            // G4double rOuter[]		= {bottomRad + thick/2., topRad+thick/2.};  //the outer radius of the cone
//Outer layer
            G4Polyhedra * scintconebig = new G4Polyhedra("scintconebig", 0, 2*pi, 6, 2,
                        									zPositions1, rInnerscint, rOuterscint);    //define the cone for the outer scintillator layer
            G4Polyhedra* scintconesubtract = new G4Polyhedra("scintconesub",
                                                       0.,
                                                       2*pi,
                                                       3,
                                                       2,
                                                       zplanesg,
                                                       rscintinnerg,
                                                       rscintouterg)	;     //define the cutting out shape for outer layer

              G4VSolid* scintcone = new G4IntersectionSolid("scint_case", // name
                                                        scintconebig,  // solid A
                                                        scintconesubtract,  // solid B
                                                        rot2,            // rotation
                                                        G4ThreeVector(0,2*x1, 0)); // translation  //cut out the side of the cone for the outer layer

              auto scintconeLV
                                                             = new G4LogicalVolume(
                                                                          scintcone,     // its solid
                                                                          plasticMaterial,  // its material
                                                                          "scintconeLV");   // its name        //define logical volume for the outer layer

              for (int i = 0; i < 6; i++) {
                                                        //G4double y_angle = y_angles[i];
                                                        //G4double z_angle = z_angles[i];
                                                        rot = new G4RotationMatrix((30+(i*60))*deg,yangle,zangle);
                                                        G4double a_coord = a_coords[i];
                                                        G4double b_coord = b_coords[i];
                                                        fAbsorber1PV[i] = new G4PVPlacement(                                           //auto pconeside_physical[i] =
                                                        rot,       // rotation
                                                        G4ThreeVector(a_coord,b_coord,y),  // at Ex. G4ThreeVector((i-3)*m,(i-3)*m,2.5*m)
                                                        scintconeLV,          // its logical volume
                                                        "scintcone",    // its name
                                                        worldLV,          // its mother  volume
                                                        false,            // no boolean operation
                                                        0,                // copy number
                                                        fCheckOverlaps);  // checking overlaps
                                                      }            //Place the different sides of the outer layer


//Inner layer
            G4Polyhedra * scintconebig1 = new G4Polyhedra("scintconebig", 0, 2*pi, 6, 2,
                                          zPositions2, rInnerscint1, rOuterscint1);    //define the inner layer cone
            G4Polyhedra* scintconesubtract1 = new G4Polyhedra("scintconesub",
                                                       0.,
                                                       2*pi,
                                                       3,
                                                       2,
                                                       zplanesg,
                                                       rscintinnerg1,
                                                       rscintouterg1)	;         //define the cutting out shape

              G4VSolid* scintcone1 = new G4IntersectionSolid("scint_case", // name
                                                        scintconebig1,  // solid A
                                                        scintconesubtract1,  // solid B
                                                        rot2,            // rotation
                                                        G4ThreeVector(0,2*x2,0)); // translation    //cut out the side of the cone for the inner layer

              auto scintcone1LV
                                                             = new G4LogicalVolume(
                                                                          scintcone1,     // its solid
                                                                          plasticMaterial,  // its material
                                                                          "scintcone1LV");   // its name           //define a logical volume


              for (int i = 0; i < 6; i++) {
                                                        //G4double y_angle = y_angles[i];
                                                        //G4double z_angle = z_angles[i];
                                                        rot = new G4RotationMatrix((30+(i*60))*deg,yangle,zangle);
                                                        G4double a_coord = a_coords[i];
                                                        G4double b_coord = b_coords[i];
                                                        fAbsorberPV[i] = new G4PVPlacement(                                           //auto pconeside_physical[i] =
                                                        rot,       // rotation
                                                        G4ThreeVector(a_coord,b_coord,-(scintthick/tan(openingangle))),  // at Ex. G4ThreeVector((i-3)*m,(i-3)*m,2.5*m)
                                                        scintcone1LV,          // its logical volume
                                                        "scintcone1",    // its name
                                                        worldLV,          // its mother  volume
                                                        false,            // no boolean operation
                                                        0,                // copy number
                                                        fCheckOverlaps);  // checking overlaps
                                                      }       //place the inner sides


// Defining a rectangular polarimeter  (not in use currently)
  //Defining Geometry parameters
  G4int nofLayers = 2;
  G4double plasticThickness = 0.01*m; //7*mm-1*cm
  G4double carbonThickness =  0.02*m; //2*cm
  G4double polarSizeXY  = 0.5*m;


  auto layerThickness = plasticThickness + carbonThickness;
  auto polarThickness = (2 * plasticThickness) + (carbonThickness);

  auto polarimeterS
    = new G4Box("polarimeter",     // its name
                 polarSizeXY/2, polarSizeXY/2, polarThickness/2); // its size

  auto polarLV
    = new G4LogicalVolume(
                 polarimeterS,     // its solid
                 defaultMaterial,  // its material
                 "polarimeter");   // its name

  // new G4PVPlacement(
  //                0,                // no rotation
  //                G4ThreeVector(),  // at (0,0,0)
  //                polarLV,          // its logical volume
  //                "polarimeter",    // its name
  //                worldLV,          // its mother  volume
  //                false,            // no boolean operation
  //                0,                // copy number
  //                fCheckOverlaps);  // checking overlaps

  //
  // Creating a Generic Cuboid Layer
  //
  // auto layerS
  //   = new G4Box("Layer",           // its name
  //                polarSizeXY/2, polarSizeXY/2, layerThickness/2); // its size
  //
  // auto layerLV
  //   = new G4LogicalVolume(
  //                layerS,           // its solid
  //                defaultMaterial,  // its material
  //                "Layer");         // its name
  //
  // new G4PVReplica(
  //                "Layer",          // its name
  //                layerLV,          // its logical volume
  //                polarLV,          // its mother
  //                kZAxis,           // axis of replication
  //                nofLayers,        // number of replica
  //                layerThickness);  // witdth of replica

  //
  // plastic Layer 1
  //
  auto plasticS
    = new G4Box("plastic",            // its name
                 polarSizeXY/2, polarSizeXY/2, plasticThickness/2); // its size

  auto plasticLV
    = new G4LogicalVolume(
                 plasticS,        // its solid
                 plasticMaterial, // its material
                 "plastic");          // its name

  //plasticLV->SetUserLimits(fStepLimit);
  //
  // fAbsorberPV
  //   = new G4PVPlacement(
  //                0,                // no rotation
  //                G4ThreeVector(0., 0., (plasticThickness/2)-(polarThickness/2)), // its position
  //                plasticLV,       // its logical volume
  //                "plastic",           // its name
  //                polarLV,          // its mother  volume
  //                false,            // no boolean operation
  //                0,                // copy number
  //                fCheckOverlaps);  // checking overlaps


  //
  // plastic Layer 2

  auto plastic1S
    = new G4Box("plastic 1",            // its name
                 polarSizeXY/2, polarSizeXY/2, plasticThickness/2); // its size

  auto plastic1LV
    = new G4LogicalVolume(
                 plastic1S,        // its solid
                 plasticMaterial, // its material
                 "plastic 1");          // its name

  //plastic1LV->SetUserLimits(fStepLimit);
  //
  // fAbsorber1PV
  //   = new G4PVPlacement(
  //                0,                // no rotation
  //                G4ThreeVector(0., 0., -(plasticThickness/2)+(polarThickness/2)), // its position
  //                plastic1LV,       // its logical volume
  //                "plastic1",           // its name
  //                polarLV,          // its mother  volume
  //                false,            // no boolean operation
  //                0,                // copy number
  //                fCheckOverlaps);  // checking overlaps

  //
  // carbon
  //
  auto carbonS
    = new G4Box("carbon",             // its name
                 polarSizeXY/2, polarSizeXY/2, carbonThickness/2); // its size

  auto carbonLV
    = new G4LogicalVolume(
                 carbonS,             // its solid
                 carbonMaterial,      // its material
                 "carbon");           // its name

  //carbonLV->SetUserLimits(fStepLimit);
  //
  // fGapPV
  //   = new G4PVPlacement(
  //                0,                // no rotation
  //                G4ThreeVector(0., 0., 0.), // its position
  //                carbonLV,            // its logical volume
  //                "carbon",            // its name
  //                polarLV,          // its mother  volume
  //                false,            // no boolean operation
  //                0,                // copy number
  //                fCheckOverlaps);  // checking overlaps
  //
  // //
  // print parameters
  //
  G4cout
    << G4endl
    << "------------------------------------------------------------" << G4endl
    << "---> The polarimeter is " << nofLayers << " layers of: [ "
    << plasticThickness/mm << "mm of " << plasticMaterial->GetName()
    << " ] +  1 layer of ["
    << carbonThickness/mm << "mm of " << carbonMaterial->GetName() << " ] " << G4endl
    << "------------------------------------------------------------" << G4endl;

  //
  // Visualization attributes
  //
  worldLV->SetVisAttributes (G4Colour(0,0,0.5));

  auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(0,0,1.0));   //blue
  simpleBoxVisAtt->SetVisibility(true);
  simpleBoxVisAtt->SetForceSolid(true);
  // carbonLV->SetVisAttributes(simpleBoxVisAtt);  //for rectangle polarimeter
  scintconeLV->SetVisAttributes(simpleBoxVisAtt);
  scintcone1LV->SetVisAttributes(simpleBoxVisAtt);


  auto redsimpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,0,0));    //red
  redsimpleBoxVisAtt->SetVisibility(true);
  redsimpleBoxVisAtt->SetForceSolid(true);
  // plasticLV->SetVisAttributes(redsimpleBoxVisAtt);    //for rectangle polarimeter
  // plastic1LV->SetVisAttributes(redsimpleBoxVisAtt);    //for rectangle polarimeter

  pconeLV->SetVisAttributes(redsimpleBoxVisAtt);


  //
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4DetectorConstruction::ConstructSDandField()
{
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue;
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);

  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
