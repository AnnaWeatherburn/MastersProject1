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
/// \file B4DetectorConstruction.hh
/// \brief Definition of the B4DetectorConstruction class

#ifndef B4DetectorConstruction_h
#define B4DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4GlobalMagFieldMessenger;
class G4UserLimits;

/// Detector construction class to define materials and geometry.
/// The calorimeter is a box made of a given number of layers. A layer consists
/// of an absorber plate and of a detection gap. The layer is replicated.
///
/// Four parameters define the geometry of the calorimeter :
///
/// - the thickness of an absorber plate,
/// - the thickness of a gap,
/// - the number of layers,
/// - the transverse size of the calorimeter (the input face is a square).
///
/// In addition a transverse uniform magnetic field is defined
/// via G4GlobalMagFieldMessenger class.

class B4DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    B4DetectorConstruction();
    virtual ~B4DetectorConstruction();

  public:
    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

    // get methods
    //
    const G4VPhysicalVolume* GetAbsorberPV() const;
    const G4VPhysicalVolume* GetGapPV() const;
    const G4VPhysicalVolume* GetAbsorber1PV() const;
    const G4VPhysicalVolume* GetAbsorberPV_1() const;
    const G4VPhysicalVolume* GetAbsorberPV_2() const;
    const G4VPhysicalVolume* GetAbsorberPV_3() const;
    const G4VPhysicalVolume* GetAbsorberPV_4() const;
    const G4VPhysicalVolume* GetAbsorberPV_5() const;
    const G4VPhysicalVolume* GetAbsorberPV_6() const;
    const G4VPhysicalVolume* GetGapPV_1() const;
    const G4VPhysicalVolume* GetGapPV_2() const;
    const G4VPhysicalVolume* GetGapPV_3() const;
    const G4VPhysicalVolume* GetGapPV_4() const;
    const G4VPhysicalVolume* GetGapPV_5() const;
    const G4VPhysicalVolume* GetGapPV_6() const;
    const G4VPhysicalVolume* GetAbsorber1PV_1() const;
    const G4VPhysicalVolume* GetAbsorber1PV_2() const;
    const G4VPhysicalVolume* GetAbsorber1PV_3() const;
    const G4VPhysicalVolume* GetAbsorber1PV_4() const;
    const G4VPhysicalVolume* GetAbsorber1PV_5() const;
    const G4VPhysicalVolume* GetAbsorber1PV_6() const;


  private:
    // methods
    //
    void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes();

    // data members
    //
    static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger;
                                      // magnetic field messenger
    G4VPhysicalVolume*   fAbsorberPV[6];
    //G4VPhysicalVolume*   fGapPV;
    G4VPhysicalVolume*   fGapPV[6];
    G4VPhysicalVolume*   fAbsorber1PV[6];
    // G4VPhysicalVolume*   fAbsorberPV_1; // the absorber physical volume
    // G4VPhysicalVolume*   fAbsorberPV_2; // the absorber physical volume
    // G4VPhysicalVolume*   fAbsorberPV_3; // the absorber physical volume
    // G4VPhysicalVolume*   fAbsorberPV_4; // the absorber physical volume
    // G4VPhysicalVolume*   fAbsorberPV_5; // the absorber physical volume
    // G4VPhysicalVolume*   fAbsorberPV_6; // the absorber physical volume
    //G4VPhysicalVolume*   fGapPV[0];      // the gap physical volume
    //G4VPhysicalVolume*   fGapPV[1];      // the gap physical volume
    //G4VPhysicalVolume*   fGapPV[2];      // the gap physical volume
    //G4VPhysicalVolume*   fGapPV[3];      // the gap physical volume
    //G4VPhysicalVolume*   fGapPV[4];      // the gap physical volume
    //G4VPhysicalVolume*   fGapPV[5];      // the gap physical volume
    // G4VPhysicalVolume*   fAbsorber1PV_1; // the absorber second layer physical volume
    // G4VPhysicalVolume*   fAbsorber1PV_2; // the absorber second layer physical volume
    // G4VPhysicalVolume*   fAbsorber1PV_3; // the absorber second layer physical volume
    // G4VPhysicalVolume*   fAbsorber1PV_4; // the absorber second layer physical volume
    // G4VPhysicalVolume*   fAbsorber1PV_5; // the absorber second layer physical volume
    // G4VPhysicalVolume*   fAbsorber1PV_6; // the absorber second layer physical volume

    G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps
    //G4UserLimits* fStepLimit;
};

// inline functions
// inline const G4VPhysicalVolume* B4DetectorConstruction::GetAbsorberPV() const {
//   return fAbsorberPV;
// }

inline const G4VPhysicalVolume* B4DetectorConstruction::GetAbsorberPV_1() const {
  return fAbsorberPV[0];
}
inline const G4VPhysicalVolume* B4DetectorConstruction::GetAbsorberPV_2() const {
  return fAbsorberPV[1];
}
inline const G4VPhysicalVolume* B4DetectorConstruction::GetAbsorberPV_3() const {
  return fAbsorberPV[2];
}
inline const G4VPhysicalVolume* B4DetectorConstruction::GetAbsorberPV_4() const {
  return fAbsorberPV[3];
}
inline const G4VPhysicalVolume* B4DetectorConstruction::GetAbsorberPV_5() const {
  return fAbsorberPV[4];
}
inline const G4VPhysicalVolume* B4DetectorConstruction::GetAbsorberPV_6() const {
  return fAbsorberPV[5];
}

// inline const G4VPhysicalVolume* B4DetectorConstruction::GetGapPV() const  {
//   return fGapPV;
// }

inline const G4VPhysicalVolume* B4DetectorConstruction::GetGapPV_1() const  {
  return fGapPV[0];
}
inline const G4VPhysicalVolume* B4DetectorConstruction::GetGapPV_2() const  {
  return fGapPV[1];
}
inline const G4VPhysicalVolume* B4DetectorConstruction::GetGapPV_3() const  {
  return fGapPV[2];
}
inline const G4VPhysicalVolume* B4DetectorConstruction::GetGapPV_4() const  {
  return fGapPV[3];
}
inline const G4VPhysicalVolume* B4DetectorConstruction::GetGapPV_5() const  {
  return fGapPV[4];
}
inline const G4VPhysicalVolume* B4DetectorConstruction::GetGapPV_6() const  {
  return fGapPV[5];
}

// inline const G4VPhysicalVolume* B4DetectorConstruction::GetAbsorber1PV() const {
//
//   return fAbsorber1PV;
// }

inline const G4VPhysicalVolume* B4DetectorConstruction::GetAbsorber1PV_1() const {

  return fAbsorber1PV[0];
}

inline const G4VPhysicalVolume* B4DetectorConstruction::GetAbsorber1PV_2() const {

  return fAbsorber1PV[1];
}

inline const G4VPhysicalVolume* B4DetectorConstruction::GetAbsorber1PV_3() const {

  return fAbsorber1PV[2];
}

inline const G4VPhysicalVolume* B4DetectorConstruction::GetAbsorber1PV_4() const {

  return fAbsorber1PV[3];
}

inline const G4VPhysicalVolume* B4DetectorConstruction::GetAbsorber1PV_5() const {

  return fAbsorber1PV[4];
}

inline const G4VPhysicalVolume* B4DetectorConstruction::GetAbsorber1PV_6() const {

  return fAbsorber1PV[5];
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
