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
/// \file B4aEventAction.cc
/// \brief Implementation of the B4aEventAction class

#include "B4aEventAction.hh"
#include "B4RunAction.hh"
#include "B4Analysis.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4aEventAction::B4aEventAction()
 : G4UserEventAction(),
   fEnergyAbs(0.),
   fEnergyGap(0.),
   fTrackLAbs(0.),
   fTrackLGap(0.),
   fEnergyAbs1(0.),
   fTrackLAbs1(0.),
   fEnergyAbs_1(0.),
   fEnergyGap_1(0.),
   fTrackLAbs_1(0.),
   fTrackLGap_1(0.),
   fEnergyAbs1_1(0.),
   fTrackLAbs1_1(0.),
   fEnergyAbs_2(0.),
   fEnergyGap_2(0.),
   fTrackLAbs_2(0.),
   fTrackLGap_2(0.),
   fEnergyAbs1_2(0.),
   fTrackLAbs1_2(0.),
   fEnergyAbs_3(0.),
   fEnergyGap_3(0.),
   fTrackLAbs_3(0.),
   fTrackLGap_3(0.),
   fEnergyAbs1_3(0.),
   fTrackLAbs1_3(0.),
   fEnergyAbs_4(0.),
   fEnergyGap_4(0.),
   fTrackLAbs_4(0.),
   fTrackLGap_4(0.),
   fEnergyAbs1_4(0.),
   fTrackLAbs1_4(0.),
   fEnergyAbs_5(0.),
   fEnergyGap_5(0.),
   fTrackLAbs_5(0.),
   fTrackLGap_5(0.),
   fEnergyAbs1_5(0.),
   fTrackLAbs1_5(0.),
   fEnergyAbs_6(0.),
   fEnergyGap_6(0.),
   fTrackLAbs_6(0.),
   fTrackLGap_6(0.),
   fEnergyAbs1_6(0.),
   fTrackLAbs1_6(0.),
   fMomentumX(0.), //*******************************************
   fMomentumY(0.), //*******************************************
   fMomentumZ(0.), //*******************************************
   fEnergyParticle(0.), //*******************************************
   fParticlePosX(0.), //*******************************************
   fParticlePosY(0.), //*******************************************
   fParticlePosZ(0.) //***************************************
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4aEventAction::~B4aEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4aEventAction::BeginOfEventAction(const G4Event* /*event*/)
{
  // initialisation per event
  fEnergyAbs = 0.;
  fEnergyGap = 0.;
  fTrackLAbs = 0.;
  fTrackLGap = 0.;
  fEnergyAbs1 = 0.;
  fTrackLAbs1 = 0.;
  fEnergyAbs_1 = 0.;
  fEnergyGap_1 = 0.;
  fTrackLAbs_1 = 0.;
  fTrackLGap_1 = 0.;
  fEnergyAbs1_1 = 0.;
  fTrackLAbs1_1 = 0.;
  fEnergyAbs_2 = 0.;
  fEnergyGap_2 = 0.;
  fTrackLAbs_2 = 0.;
  fTrackLGap_2 = 0.;
  fEnergyAbs1_2 = 0.;
  fTrackLAbs1_2 = 0.;
  fEnergyAbs_3 = 0.;
  fEnergyGap_3 = 0.;
  fTrackLAbs_3 = 0.;
  fTrackLGap_3 = 0.;
  fEnergyAbs1_3 = 0.;
  fTrackLAbs1_3 = 0.;
  fEnergyAbs_4 = 0.;
  fEnergyGap_4 = 0.;
  fTrackLAbs_4 = 0.;
  fTrackLGap_4 = 0.;
  fEnergyAbs1_4 = 0.;
  fTrackLAbs1_4 = 0.;
  fEnergyAbs_5 = 0.;
  fEnergyGap_5 = 0.;
  fTrackLAbs_5 = 0.;
  fTrackLGap_5 = 0.;
  fEnergyAbs1_5 = 0.;
  fTrackLAbs1_5 = 0.;
  fEnergyAbs_6 = 0.;
  fEnergyGap_6 = 0.;
  fTrackLAbs_6 = 0.;
  fTrackLGap_6 = 0.;
  fEnergyAbs1_6 = 0.;
  fTrackLAbs1_6 = 0.;
  fMomentumX = 0.; //*******************************************
  fMomentumY = 0.; //*******************************************
  fMomentumZ = 0.; //*******************************************
  fEnergyParticle = 0.; //*******************************************
  fParticlePosX = 0.; //*******************************************
  fParticlePosY = 0.; //*******************************************
  fParticlePosZ = 0.; //*********************************************
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4aEventAction::EndOfEventAction(const G4Event* event)
{
  // Accumulate statistics
  //

  // get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  G4double totalEnergyfirstScint = fEnergyAbs_1 + fEnergyAbs_2 + fEnergyAbs_3 + fEnergyAbs_4 + fEnergyAbs_5 + fEnergyAbs_6 ;
  G4double totalEnergyCarbon = fEnergyGap_1 + fEnergyGap_2 + fEnergyGap_3 + fEnergyGap_4 + fEnergyGap_5 + fEnergyGap_6 ;
  G4double totalEnergysecondScint = fEnergyAbs1_1 + fEnergyAbs1_2 + fEnergyAbs1_3 + fEnergyAbs1_4 + fEnergyAbs1_5 + fEnergyAbs1_6 ;


  // fill histograms
  analysisManager->FillH1(0, fEnergyAbs);
  analysisManager->FillH1(1, fEnergyGap);
  analysisManager->FillH1(2, fTrackLAbs);
  analysisManager->FillH1(3, fTrackLGap);
  analysisManager->FillH1(4, fEnergyAbs1);
  analysisManager->FillH1(5, fTrackLAbs1);
  analysisManager->FillH1(6, fMomentumX); //*******************************************
  analysisManager->FillH1(7, fMomentumY); //*******************************************
  analysisManager->FillH1(8, fMomentumZ); //*******************************************
  analysisManager->FillH1(9, fEnergyParticle); //*******************************************
  analysisManager->FillH1(10, fParticlePosX); //*******************************************
  analysisManager->FillH1(11, fParticlePosY); //*******************************************
  analysisManager->FillH1(12, fParticlePosZ); //*******************************************
  analysisManager->FillH1(13, fEnergyAbs_1);
  analysisManager->FillH1(14, fEnergyGap_1);
  analysisManager->FillH1(15, fTrackLAbs_1);
  analysisManager->FillH1(16, fTrackLGap_1);
  analysisManager->FillH1(17, fEnergyAbs1_1);
  analysisManager->FillH1(18, fTrackLAbs1_1);
  analysisManager->FillH1(19, fEnergyAbs_2);
  analysisManager->FillH1(20, fEnergyGap_2);
  analysisManager->FillH1(21, fTrackLAbs_2);
  analysisManager->FillH1(22, fTrackLGap_2);
  analysisManager->FillH1(23, fEnergyAbs1_2);
  analysisManager->FillH1(24, fTrackLAbs1_2);
  analysisManager->FillH1(25, fEnergyAbs_3);
  analysisManager->FillH1(26, fEnergyGap_3);
  analysisManager->FillH1(27, fTrackLAbs_3);
  analysisManager->FillH1(28, fTrackLGap_3);
  analysisManager->FillH1(29, fEnergyAbs1_3);
  analysisManager->FillH1(30, fTrackLAbs1_3);
  analysisManager->FillH1(31, fEnergyAbs_4);
  analysisManager->FillH1(32, fEnergyGap_4);
  analysisManager->FillH1(33, fTrackLAbs_4);
  analysisManager->FillH1(34, fTrackLGap_4);
  analysisManager->FillH1(35, fEnergyAbs1_4);
  analysisManager->FillH1(36, fTrackLAbs1_4);
  analysisManager->FillH1(37, fEnergyAbs_5);
  analysisManager->FillH1(38, fEnergyGap_5);
  analysisManager->FillH1(39, fTrackLAbs_5);
  analysisManager->FillH1(40, fTrackLGap_5);
  analysisManager->FillH1(41, fEnergyAbs1_5);
  analysisManager->FillH1(42, fTrackLAbs1_5);
  analysisManager->FillH1(43, fEnergyAbs_6);
  analysisManager->FillH1(44, fEnergyGap_6);
  analysisManager->FillH1(45, fTrackLAbs_6);
  analysisManager->FillH1(46, fTrackLGap_6);
  analysisManager->FillH1(47, fEnergyAbs1_6);
  analysisManager->FillH1(48, fTrackLAbs1_6);

  if (( fEnergyAbs <= 0.5*MeV ) && ( fEnergyAbs1 >= 0.5*MeV) ) {
    G4cout << "Fits the criteria of < 0.5MeV deposited in the first layer and more than 0.5MeV deposited in the second layer" << G4endl;
    analysisManager->FillH1(49, fEnergyAbs);
    analysisManager->FillH1(50, fEnergyAbs1);
        }

  else {
      }

  analysisManager->FillH1(51, totalEnergyfirstScint);
  analysisManager->FillH1(52, totalEnergyCarbon);
  analysisManager->FillH1(53, totalEnergysecondScint);
  // fill ntuple
  analysisManager->FillNtupleDColumn(0, fEnergyAbs);
  analysisManager->FillNtupleDColumn(1, fEnergyGap);
  analysisManager->FillNtupleDColumn(2, fTrackLAbs);
  analysisManager->FillNtupleDColumn(3, fTrackLGap);
  analysisManager->FillNtupleDColumn(4, fEnergyAbs1);
  analysisManager->FillNtupleDColumn(5, fTrackLAbs1);
  analysisManager->FillNtupleDColumn(6, fMomentumX); //*******************************************
  analysisManager->FillNtupleDColumn(7, fMomentumY); //*******************************************
  analysisManager->FillNtupleDColumn(8, fMomentumZ); //*******************************************
  analysisManager->FillNtupleDColumn(9, fEnergyParticle); //*******************************************
  analysisManager->FillNtupleDColumn(10, fParticlePosX); //*******************************************
  analysisManager->FillNtupleDColumn(11, fParticlePosY); //*******************************************
  analysisManager->FillNtupleDColumn(12, fParticlePosZ); //*******************************************
  analysisManager->FillNtupleDColumn(13, fEnergyAbs_1);
  analysisManager->FillNtupleDColumn(14, fEnergyGap_1);
  analysisManager->FillNtupleDColumn(15, fTrackLAbs_1);
  analysisManager->FillNtupleDColumn(16, fTrackLGap_1);
  analysisManager->FillNtupleDColumn(17, fEnergyAbs1_1);
  analysisManager->FillNtupleDColumn(18, fTrackLAbs1_1);
  analysisManager->FillNtupleDColumn(19, fEnergyAbs_2);
  analysisManager->FillNtupleDColumn(20, fEnergyGap_2);
  analysisManager->FillNtupleDColumn(21, fTrackLAbs_2);
  analysisManager->FillNtupleDColumn(22, fTrackLGap_2);
  analysisManager->FillNtupleDColumn(23, fEnergyAbs1_2);
  analysisManager->FillNtupleDColumn(24, fTrackLAbs1_2);
  analysisManager->FillNtupleDColumn(25, fEnergyAbs_3);
  analysisManager->FillNtupleDColumn(26, fEnergyGap_3);
  analysisManager->FillNtupleDColumn(27, fTrackLAbs_3);
  analysisManager->FillNtupleDColumn(28, fTrackLGap_3);
  analysisManager->FillNtupleDColumn(29, fEnergyAbs1_3);
  analysisManager->FillNtupleDColumn(30, fTrackLAbs1_3);
  analysisManager->FillNtupleDColumn(31, fEnergyAbs_4);
  analysisManager->FillNtupleDColumn(32, fEnergyGap_4);
  analysisManager->FillNtupleDColumn(33, fTrackLAbs_4);
  analysisManager->FillNtupleDColumn(34, fTrackLGap_4);
  analysisManager->FillNtupleDColumn(35, fEnergyAbs1_4);
  analysisManager->FillNtupleDColumn(36, fTrackLAbs1_4);
  analysisManager->FillNtupleDColumn(37, fEnergyAbs_5);
  analysisManager->FillNtupleDColumn(38, fEnergyGap_5);
  analysisManager->FillNtupleDColumn(39, fTrackLAbs_5);
  analysisManager->FillNtupleDColumn(40, fTrackLGap_5);
  analysisManager->FillNtupleDColumn(41, fEnergyAbs1_5);
  analysisManager->FillNtupleDColumn(42, fTrackLAbs1_5);
  analysisManager->FillNtupleDColumn(43, fEnergyAbs_6);
  analysisManager->FillNtupleDColumn(44, fEnergyGap_6);
  analysisManager->FillNtupleDColumn(45, fTrackLAbs_6);
  analysisManager->FillNtupleDColumn(46, fTrackLGap_6);
  analysisManager->FillNtupleDColumn(47, fEnergyAbs1_6);
  analysisManager->FillNtupleDColumn(48, fTrackLAbs1_6);
  analysisManager->FillNtupleDColumn(49, fEnergyAbs); //*******************************************
  analysisManager->FillNtupleDColumn(50, fEnergyAbs1); //*******************************************
  analysisManager->FillNtupleDColumn(51, totalEnergyfirstScint);
  analysisManager->FillNtupleDColumn(52, totalEnergyCarbon);
  analysisManager->FillNtupleDColumn(53, totalEnergysecondScint);
  //analysisManager->AddNtupleRow();******************************************
  analysisManager->AddNtupleRow();

  // Print per event (modulo n)
  //
  auto eventID = event->GetEventID();
  auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {
    //G4cout << "---> End of event: " << eventID << G4endl;
    //
    // G4cout
    //    << "   Plastic 1: total energy: " << std::setw(7)
    //                                     << G4BestUnit(fEnergyAbs,"Energy")
    //    << "       total track length: " << std::setw(7)
    //                                     << G4BestUnit(fTrackLAbs,"Length")
    //    << G4endl
    //    << "   Plastic 2: total energy: " << std::setw(7)
    //                                     << G4BestUnit(fEnergyAbs1,"Energy")
    //    << "       total track length: " << std::setw(7)
    //                                     << G4BestUnit(fTrackLAbs1,"Length")
    //    << G4endl
    //    << "        Carbon Layer: total energy: " << std::setw(7)
    //                                     << G4BestUnit(fEnergyGap,"Energy")
    //    << "       total track length: " << std::setw(7)
    //                                     << G4BestUnit(fTrackLGap,"Length")
    //    << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
