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
/// \file B4aSteppingAction.cc
/// \brief Implementation of the B4aSteppingAction class

#include "B4aSteppingAction.hh"
#include "B4aEventAction.hh"
#include "B4DetectorConstruction.hh"


#include "G4Step.hh"
#include "G4RunManager.hh"
#include "iostream"
#include "fstream"
using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4aSteppingAction::B4aSteppingAction(
                      const B4DetectorConstruction* detectorConstruction,
                      B4aEventAction* eventAction)
  : G4UserSteppingAction(),
    fDetConstruction(detectorConstruction),
    fEventAction(eventAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4aSteppingAction::~B4aSteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4aSteppingAction::UserSteppingAction(const G4Step* step)
{
// Collect energy and track length step by step
//G4cout <<" Stepping ACtion" << G4endl;
  // get volume of the current step
  auto volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

  // energy deposit
  auto edep = step->GetTotalEnergyDeposit();

  // step length
  G4double stepLength = 0.;

  G4Track* aTrack = step->GetTrack();

  if (aTrack->GetDefinition()->GetParticleName() == "neutron" || aTrack->GetDefinition()->GetParticleName() == "proton" ) { //&&
          stepLength = step->GetStepLength();
          G4cout << "Step Length=" << stepLength << G4endl;
         } //changed this to be dependent on particle type as needed
//stepLength = step->GetStepLength();
//G4cout << "Step Length=" << stepLength << G4endl;

  // if ( step->GetTrack()->GetDefinition()->GetPDGCharge() != 1. ) {
  //   stepLength = step->GetStepLength();
  // }

  //stepLength = step->GetStepLength();

  // if ( volume == fDetConstruction->GetAbsorberPV() ) {
  //   fEventAction->AddAbs(edep,stepLength);
  //   G4cout <<"First Scintillator 0" << G4endl;
  // }

  if ( volume == fDetConstruction->GetAbsorberPV_1() ) {
    fEventAction->AddAbs_1(edep,stepLength);
    G4cout <<"First (Inner) Scintillator 1" << G4endl;
  }
  if ( volume == fDetConstruction->GetAbsorberPV_2() ) {
    fEventAction->AddAbs_2(edep,stepLength);
    G4cout <<"First Scintillator 2" << G4endl;
  }
  if ( volume == fDetConstruction->GetAbsorberPV_3() ) {
    fEventAction->AddAbs_3(edep,stepLength);
    G4cout <<"First Scintillator 3" << G4endl;
  }
  if ( volume == fDetConstruction->GetAbsorberPV_4() ) {
    fEventAction->AddAbs_4(edep,stepLength);
    G4cout <<"First Scintillator 4" << G4endl;
  }
  if ( volume == fDetConstruction->GetAbsorberPV_5() ) {
    fEventAction->AddAbs_5(edep,stepLength);
    G4cout <<"First Scintillator 5" << G4endl;
  }
  if ( volume == fDetConstruction->GetAbsorberPV_6() ) {
    fEventAction->AddAbs_6(edep,stepLength);
    G4cout <<"First Scintillator 6" << G4endl;
  }

  // if ( volume == fDetConstruction->GetAbsorber1PV() ) {
  //   fEventAction->AddAbs1(edep,stepLength);
  //   G4cout <<"Second Scintillator 0" << G4endl;
  // }

  if ( volume == fDetConstruction->GetAbsorber1PV_1() ) {
    fEventAction->AddAbs1_1(edep,stepLength);
    G4cout <<"Second Scintillator 1" << G4endl;
  }
  if ( volume == fDetConstruction->GetAbsorber1PV_2() ) {
    fEventAction->AddAbs1_2(edep,stepLength);
    G4cout <<"Second Scintillator 2" << G4endl;
  }
  if ( volume == fDetConstruction->GetAbsorber1PV_3() ) {
    fEventAction->AddAbs1_3(edep,stepLength);
    G4cout <<"Second Scintillator 3" << G4endl;
  }
  if ( volume == fDetConstruction->GetAbsorber1PV_4() ) {
    fEventAction->AddAbs1_4(edep,stepLength);
    G4cout <<"Second Scintillator 4" << G4endl;
  }
  if ( volume == fDetConstruction->GetAbsorber1PV_5() ) {
    fEventAction->AddAbs1_5(edep,stepLength);
    G4cout <<"Second Scintillator 5" << G4endl;
  }
  if ( volume == fDetConstruction->GetAbsorber1PV_6() ) {
    fEventAction->AddAbs1_6(edep,stepLength);
    G4cout <<"Second Scintillator 6" << G4endl;
  }

  // if ( volume == fDetConstruction->GetGapPV() ) {
  //   fEventAction->AddGap(edep,stepLength);
  //   G4cout <<"Carbon 0" << G4endl;
  // }

  if ( volume == fDetConstruction->GetGapPV_1() ) {
    fEventAction->AddGap_1(edep,stepLength);
    G4cout <<"Carbon 1" << G4endl;
  }
  if ( volume == fDetConstruction->GetGapPV_2() ) {
    fEventAction->AddGap_2(edep,stepLength);
    G4cout <<"Carbon 2" << G4endl;
  }
  if ( volume == fDetConstruction->GetGapPV_3() ) {
    fEventAction->AddGap_3(edep,stepLength);
    G4cout <<"Carbon 3" << G4endl;
  }
  if ( volume == fDetConstruction->GetGapPV_4() ) {
    fEventAction->AddGap_4(edep,stepLength);
    G4cout <<"Carbon 4" << G4endl;
  }
  if ( volume == fDetConstruction->GetGapPV_5() ) {
    fEventAction->AddGap_5(edep,stepLength);
    G4cout <<"Carbon 5" << G4endl;
  }
  if ( volume == fDetConstruction->GetGapPV_6() ) {
    fEventAction->AddGap_6(edep,stepLength);
    G4cout <<"Carbon 6" << G4endl;
  }

  //	Check if particle is about to leave the polarimeter and save the momentum



  // Create and open a text file
  //ofstream MyFile("filename.txt");
//  G4String fileName1 = "POutput"; **********************
//  analysisManager->OpenFile(fileName1);**********************

	if (aTrack->GetDefinition()->GetParticleName() == "neutron") { //*******************************************

    		if (aTrack->GetTrackStatus() != fStopAndKill) { //*******************************************
    			if (step->GetPreStepPoint()->GetMaterial()->GetName() == "Scintillator") { //*******************************************
            if (step->GetPostStepPoint()->GetMaterial()->GetName() == "Galactic") { //*******************************************
              G4ThreeVector MomentumVector = step->GetPreStepPoint()->GetMomentum(); //*******************************************
              G4double momentumx = MomentumVector(0);
              G4double momentumy = MomentumVector(1);
              G4double momentumz = MomentumVector(2);
              G4ThreeVector particleposVector = step->GetPreStepPoint()->GetPosition (); //****************
              G4double particleposx = particleposVector(0);
              G4double particleposy = particleposVector(1);
              G4double particleposz = particleposVector(2);
              G4double energyparticle = step->GetPreStepPoint()->GetKineticEnergy (); //****************
              fEventAction->AddLorentz(momentumx, momentumy, momentumz, energyparticle, particleposx, particleposy, particleposz); //*******************************************
              G4cout <<"P= "<< MomentumVector << " ///  px= "<< momentumx << " ///  E= " << energyparticle << G4endl << G4endl; //*******************************************

    			}

          // else  {
          //   G4cout << "did something else" << G4endl << G4endl;
          // }
        }
        }

  }
  // Close the file
  //MyFile.close();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
