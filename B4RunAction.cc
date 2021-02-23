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
/// \file B4RunAction.cc
/// \brief Implementation of the B4RunAction class

#include "B4RunAction.hh"
#include "B4Analysis.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4RunAction::B4RunAction()
 : G4UserRunAction()
{
  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);

  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in B4Analysis.hh
  auto analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Create directories
  //analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetNtupleMerging(true);
    // Note: merging ntuples is available only with Root output

  // Book histograms, ntuple
  //

  // Creating histograms
  analysisManager->CreateH1("Eabs","Edep in plastic layer 1", 2000, 0.*MeV, 10*MeV);
  analysisManager->CreateH1("Egap","Edep in carbon layer", 2000, 0.*MeV, 100*MeV);
  analysisManager->CreateH1("Labs","trackL in plastic layer 1", 5000, 0.*m, 0.1*m);
  analysisManager->CreateH1("Lgap","trackL in carbon layer", 5000, 0.*m, 0.2*m);
  analysisManager->CreateH1("Eabs1","Edep in plastic layer 2", 2000, -0.05*MeV, 100*MeV);
  analysisManager->CreateH1("Labs1","trackL in plastic layer 2", 5000, -0.05*m, 0.1*m);
  analysisManager->CreateH1("MomentumX", "X Momentum leaving polarimeter", 500, 0., 50.); //*******************************************
  analysisManager->CreateH1("MomentumY", "X Momentum leaving polarimeter", 500, 0., 50.); //*******************************************
  analysisManager->CreateH1("MomentumZ", "X Momentum leaving polarimeter", 500, 0., 50.);  //*******************************************
  analysisManager->CreateH1("ParticleEnergy", "Kinetic Energy leaving polarimeter", 500, 0., 100*MeV);  //*******************************************
  analysisManager->CreateH1("ParticleExitPositionX", "X position leaving polarimeter", 500, 0., 50.);
  analysisManager->CreateH1("ParticleExitPositionY", "Y position leaving polarimeter", 500, 0., 50.);
  analysisManager->CreateH1("ParticleExitPositionZ", "Z position leaving polarimeter", 500, 0., 100.);
  analysisManager->CreateH1("Eabs_1","Edep in plastic section 1 layer 1", 200, 0.*MeV, 100*MeV);
  analysisManager->CreateH1("Egap_1","Edep in carbon section 1", 200, 0.*MeV, 100*MeV);
  analysisManager->CreateH1("Labs_1","trackL in plastic section 1 layer 1", 500, 0.*m, 0.1*m);
  analysisManager->CreateH1("Lgap_1","trackL in carbon section 1", 500, 0.*m, 0.2*m);
  analysisManager->CreateH1("Eabs1_1","Edep in plastic section 1 layer 2", 200, -0.05*MeV, 100*MeV);
  analysisManager->CreateH1("Labs1_1","trackL in plastic section 1 layer 2", 500, -0.05*m, 0.1*m);
  analysisManager->CreateH1("Eabs_2","Edep in plastic section 2 layer 1", 200, 0.*MeV, 100*MeV);
  analysisManager->CreateH1("Egap_2","Edep in carbon section 2", 200, 0.*MeV, 100*MeV);
  analysisManager->CreateH1("Labs_2","trackL in plastic section 2 layer 1", 500, 0.*m, 0.1*m);
  analysisManager->CreateH1("Lgap_2","trackL in carbon section 2 ", 500, 0.*m, 0.2*m);
  analysisManager->CreateH1("Eabs1_2","Edep in plastic section 2 layer 2", 200, -0.05*MeV, 100*MeV);
  analysisManager->CreateH1("Labs1_2","trackL in plastic section 2 layer 2", 500, -0.05*m, 0.1*m);
  analysisManager->CreateH1("Eabs_3","Edep in plastic section 3 layer 1", 200, 0.*MeV, 100*MeV);
  analysisManager->CreateH1("Egap_3","Edep in carbon section 3", 200, 0.*MeV, 100*MeV);
  analysisManager->CreateH1("Labs_3","trackL in plastic section 3 layer 1", 500, 0.*m, 0.1*m);
  analysisManager->CreateH1("Lgap_3","trackL in carbon section 3", 500, 0.*m, 0.2*m);
  analysisManager->CreateH1("Eabs1_3","Edep in plastic section 3 layer 2", 200, -0.05*MeV, 100*MeV);
  analysisManager->CreateH1("Labs1_3","trackL in plastic section 3 layer 2", 500, -0.05*m, 0.1*m);
  analysisManager->CreateH1("Eabs_4","Edep in plastic section 4 layer 1", 200, 0.*MeV, 100*MeV);
  analysisManager->CreateH1("Egap_4","Edep in carbon section 4", 200, 0.*MeV, 100*MeV);
  analysisManager->CreateH1("Labs_4","trackL in plastic section 4 layer 1", 500, 0.*m, 0.1*m);
  analysisManager->CreateH1("Lgap_4","trackL in carbon section 4", 500, 0.*m, 0.2*m);
  analysisManager->CreateH1("Eabs1_4","Edep in plastic section 4 layer 2", 200, -0.05*MeV, 100*MeV);
  analysisManager->CreateH1("Labs1_4","trackL in plastic section 4 layer 2", 500, -0.05*m, 0.1*m);
  analysisManager->CreateH1("Eabs_5","Edep in plastic section 5 layer 1", 200, 0.*MeV, 100*MeV);
  analysisManager->CreateH1("Egap_5","Edep in carbon section 5", 200, 0.*MeV, 100*MeV);
  analysisManager->CreateH1("Labs_5","trackL in plastic section 5 layer 1", 500, 0.*m, 0.1*m);
  analysisManager->CreateH1("Lgap_5","trackL in carbon section 5", 500, 0.*m, 0.2*m);
  analysisManager->CreateH1("Eabs1_5","Edep in plastic section 5 layer 2", 200, -0.05*MeV, 100*MeV);
  analysisManager->CreateH1("Labs1_5","trackL in plastic section 5 layer 2", 500, -0.05*m, 0.1*m);
  analysisManager->CreateH1("Eabs_6","Edep in plastic section 6 layer 1", 200, 0.*MeV, 100*MeV);
  analysisManager->CreateH1("Egap_6","Edep in carbon section 6", 200, 0.*MeV, 100*MeV);
  analysisManager->CreateH1("Labs_6","trackL in plastic section 6 layer 1", 500, 0.*m, 0.1*m);
  analysisManager->CreateH1("Lgap_6","trackL in carbon section 6", 500, 0.*m, 0.2*m);
  analysisManager->CreateH1("Eabs1_6","Edep in plastic section 6 layer 2", 200, -0.05*MeV, 100*MeV);
  analysisManager->CreateH1("Labs1_6","trackL in plastic section 6 layer 2", 500, -0.05*m, 0.1*m);
  analysisManager->CreateH1("EScintillatorCriteria", "E dep in first layer (E<0.5MeV)", 500, 0., 50.);
  analysisManager->CreateH1("EScintillator1Criteria", "E dep in second layer (E>0.5MeV)", 500, 0., 1000.);
  analysisManager->CreateH1("totalEnergyfirstScint", "Total Energy First Scintillator", 500, 0., 1000.);
  analysisManager->CreateH1("totalEnergyCarbon", "Total Energy Carbon", 500, 0., 1000.);
  analysisManager->CreateH1("totalEnergysecondScint", "Total Energy Second Scintillator", 500, 0., 1000.);

  // Creating ntuple
  //
  analysisManager->CreateNtuple("B4", "Edep and TrackL");
  analysisManager->CreateNtupleDColumn("Eabs");
  analysisManager->CreateNtupleDColumn("Egap");
  analysisManager->CreateNtupleDColumn("Labs");
  analysisManager->CreateNtupleDColumn("Lgap");
  analysisManager->CreateNtupleDColumn("Eabs1");
  analysisManager->CreateNtupleDColumn("Labs1");
  analysisManager->CreateNtupleDColumn("Momentum"); //*
  analysisManager->CreateNtupleDColumn("MomentumX"); //*******************************************
  analysisManager->CreateNtupleDColumn("MomentumY"); //*******************************************
  analysisManager->CreateNtupleDColumn("MomentumZ"); //*******************************************
  analysisManager->CreateNtupleDColumn("ParticleEnergy"); //*******************************************
  analysisManager->CreateNtupleDColumn("ParticleExitPositionX"); //*******************************************
  analysisManager->CreateNtupleDColumn("ParticleExitPositionY"); //*******************************************
  analysisManager->CreateNtupleDColumn("ParticleExitPositionZ"); //*******************************************
  analysisManager->CreateNtupleDColumn("Eabs_1");
  analysisManager->CreateNtupleDColumn("Egap_1");
  analysisManager->CreateNtupleDColumn("Labs_1");
  analysisManager->CreateNtupleDColumn("Lgap_1");
  analysisManager->CreateNtupleDColumn("Eabs1_1");
  analysisManager->CreateNtupleDColumn("Labs1_1");
  analysisManager->CreateNtupleDColumn("Eabs_2");
  analysisManager->CreateNtupleDColumn("Egap_2");
  analysisManager->CreateNtupleDColumn("Labs_2");
  analysisManager->CreateNtupleDColumn("Lgap_2");
  analysisManager->CreateNtupleDColumn("Eabs1_2");
  analysisManager->CreateNtupleDColumn("Labs1_2");
  analysisManager->CreateNtupleDColumn("Eabs_3");
  analysisManager->CreateNtupleDColumn("Egap_3");
  analysisManager->CreateNtupleDColumn("Labs_3");
  analysisManager->CreateNtupleDColumn("Lgap_3");
  analysisManager->CreateNtupleDColumn("Eabs1_3");
  analysisManager->CreateNtupleDColumn("Labs1_3");
  analysisManager->CreateNtupleDColumn("Eabs_4");
  analysisManager->CreateNtupleDColumn("Egap_4");
  analysisManager->CreateNtupleDColumn("Labs_4");
  analysisManager->CreateNtupleDColumn("Lgap_4");
  analysisManager->CreateNtupleDColumn("Eabs1_4");
  analysisManager->CreateNtupleDColumn("Labs1_4");
  analysisManager->CreateNtupleDColumn("Eabs_5");
  analysisManager->CreateNtupleDColumn("Egap_5");
  analysisManager->CreateNtupleDColumn("Labs_5");
  analysisManager->CreateNtupleDColumn("Lgap_5");
  analysisManager->CreateNtupleDColumn("Eabs1_5");
  analysisManager->CreateNtupleDColumn("Labs1_5");
  analysisManager->CreateNtupleDColumn("Eabs_6");
  analysisManager->CreateNtupleDColumn("Egap_6");
  analysisManager->CreateNtupleDColumn("Labs_6");
  analysisManager->CreateNtupleDColumn("Lgap_6");
  analysisManager->CreateNtupleDColumn("Eabs1_6");
  analysisManager->CreateNtupleDColumn("Labs1_6");
  analysisManager->CreateNtupleDColumn("EScintillatorCriteria"); //*******************************************
  analysisManager->CreateNtupleDColumn("EScintillator1Criteria"); //********************************************
  analysisManager->CreateNtupleDColumn("totalEnergyfirstScint"); //********************************************
  analysisManager->CreateNtupleDColumn("totalEnergyCarbon"); //********************************************
  analysisManager->CreateNtupleDColumn("totalEnergysecondScint"); //********************************************
  analysisManager->FinishNtuple();

  //analysisManager->CreateNtuple("Polarimeter", "Momentum Vector"); //*******************************************
  //analysisManager->FinishNtuple(); //*******************************************
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4RunAction::~B4RunAction()
{
  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4RunAction::BeginOfRunAction(const G4Run* /*run*/)
{
  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);

  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Open an output file
  //
  G4String fileName = "B4";
  analysisManager->OpenFile(fileName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4RunAction::EndOfRunAction(const G4Run* /*run*/)
{
  // print histogram statistics
  //
  auto analysisManager = G4AnalysisManager::Instance();
  if ( analysisManager->GetH1(1) ) {
    G4cout << G4endl << " ----> print histograms statistic ";
    if(isMaster) {
    //  G4cout << "for the entire run " << G4endl << G4endl;
    }
    else {
    //  G4cout << "for the local thread " << G4endl << G4endl;
    }

    // G4cout << " EAbs : mean = "
    //    << G4BestUnit(analysisManager->GetH1(0)->mean(), "Energy")
    //    << " rms = "
    //    << G4BestUnit(analysisManager->GetH1(0)->rms(),  "Energy") << G4endl;
    //
    // G4cout << " EGap : mean = "
    //    << G4BestUnit(analysisManager->GetH1(1)->mean(), "Energy")
    //    << " rms = "
    //    << G4BestUnit(analysisManager->GetH1(1)->rms(),  "Energy") << G4endl;
    //
    // G4cout << " LAbs : mean = "
    //   << G4BestUnit(analysisManager->GetH1(2)->mean(), "Length")
    //   << " rms = "
    //   << G4BestUnit(analysisManager->GetH1(2)->rms(),  "Length") << G4endl;
    //
    // G4cout << " LGap : mean = "
    //   << G4BestUnit(analysisManager->GetH1(3)->mean(), "Length")
    //   << " rms = "
    //   << G4BestUnit(analysisManager->GetH1(3)->rms(),  "Length") << G4endl;
    //
    //   G4cout << " EAbs1 : mean = "
    //      << G4BestUnit(analysisManager->GetH1(4)->mean(), "Energy")
    //      << " rms = "
    //      << G4BestUnit(analysisManager->GetH1(4)->rms(),  "Energy") << G4endl;
    //
    //   G4cout << " LAbs1 : mean = "
    //     << G4BestUnit(analysisManager->GetH1(5)->mean(), "Length")
    //     << " rms = "
    //     << G4BestUnit(analysisManager->GetH1(5)->rms(),  "Length") << G4endl;

  }

  // save histograms & ntuple
  //
  analysisManager->Write();
  analysisManager->CloseFile();

//  analysisManager->OpenFile("Polarimeter"); //*******************************************
//  analysisManager->Write("Polarimeter"); //*******************************************
//  analysisManager->CloseFile(); //*******************************************
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
