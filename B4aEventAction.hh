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
/// \file B4aEventAction.hh
/// \brief Definition of the B4aEventAction class

#ifndef B4aEventAction_h
#define B4aEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "B4PrimaryGeneratorAction.hh" //->  v included myself probs not all necessary

#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh" //-> ^ included myself probs not all necessary

/// Event action class
///
/// It defines data members to hold the energy deposit and track lengths
/// of charged particles in Absober and Gap layers:
/// - fEnergyAbs, fEnergyGap, fTrackLAbs, fTrackLGap, fEnergyAbs1, fTrackLAbs1
/// which are collected step by step via the functions
/// - AddAbs(), AddGap(), AddAbs1()
class B4aEventAction : public G4UserEventAction
{
  public:
    B4aEventAction();
    virtual ~B4aEventAction();

    virtual void  BeginOfEventAction(const G4Event* event);
    virtual void    EndOfEventAction(const G4Event* event);

    void AddAbs(G4double de, G4double dl);
    void AddAbs_1(G4double de, G4double dl);
    void AddAbs_2(G4double de, G4double dl);
    void AddAbs_3(G4double de, G4double dl);
    void AddAbs_4(G4double de, G4double dl);
    void AddAbs_5(G4double de, G4double dl);
    void AddAbs_6(G4double de, G4double dl);
    void AddAbs1(G4double de, G4double dl);
    void AddAbs1_1(G4double de, G4double dl);
    void AddAbs1_2(G4double de, G4double dl);
    void AddAbs1_3(G4double de, G4double dl);
    void AddAbs1_4(G4double de, G4double dl);
    void AddAbs1_5(G4double de, G4double dl);
    void AddAbs1_6(G4double de, G4double dl);
    void AddGap(G4double de, G4double dl);
    void AddGap_1(G4double de, G4double dl);
    void AddGap_2(G4double de, G4double dl);
    void AddGap_3(G4double de, G4double dl);
    void AddGap_4(G4double de, G4double dl);
    void AddGap_5(G4double de, G4double dl);
    void AddGap_6(G4double de, G4double dl);
    void AddLorentz(G4double dpx, G4double dpy, G4double dpz, G4double dE, G4double dx, G4double dy, G4double dz); //*******************************************

  private:
    G4double  fEnergyAbs;
    G4double  fEnergyGap;
    G4double  fTrackLAbs;
    G4double  fTrackLGap;
    G4double  fEnergyAbs1;
    G4double fTrackLAbs1;
    G4double  fEnergyAbs_1;
    G4double  fEnergyGap_1;
    G4double  fTrackLAbs_1;
    G4double  fTrackLGap_1;
    G4double  fEnergyAbs1_1;
    G4double fTrackLAbs1_1;
    G4double  fEnergyAbs_2;
    G4double  fEnergyGap_2;
    G4double  fTrackLAbs_2;
    G4double  fTrackLGap_2;
    G4double  fEnergyAbs1_2;
    G4double fTrackLAbs1_2;
    G4double  fEnergyAbs_3;
    G4double  fEnergyGap_3;
    G4double  fTrackLAbs_3;
    G4double  fTrackLGap_3;
    G4double  fEnergyAbs1_3;
    G4double fTrackLAbs1_3;
    G4double  fEnergyAbs_4;
    G4double  fEnergyGap_4;
    G4double  fTrackLAbs_4;
    G4double  fTrackLGap_4;
    G4double  fEnergyAbs1_4;
    G4double fTrackLAbs1_4;
    G4double  fEnergyAbs_5;
    G4double  fEnergyGap_5;
    G4double  fTrackLAbs_5;
    G4double  fTrackLGap_5;
    G4double  fEnergyAbs1_5;
    G4double fTrackLAbs1_5;
    G4double  fEnergyAbs_6;
    G4double  fEnergyGap_6;
    G4double  fTrackLAbs_6;
    G4double  fTrackLGap_6;
    G4double  fEnergyAbs1_6;
    G4double fTrackLAbs1_6;
    G4double fMomentumX; //*******************************************
    G4double fMomentumY; //*******************************************
    G4double fMomentumZ; //*******************************************
    G4double fEnergyParticle; //*************************************
    G4double fParticlePosX; //*************************************
    G4double fParticlePosY; //*************************************
    G4double fParticlePosZ; //*******************************************

};

// inline functions

inline void B4aEventAction::AddAbs(G4double de, G4double dl) {
  fEnergyAbs += de;
  fTrackLAbs += dl;
}

inline void B4aEventAction::AddGap(G4double de, G4double dl) {
  fEnergyGap += de;
  fTrackLGap += dl;
}

inline void B4aEventAction::AddAbs1(G4double de, G4double dl) {
  fEnergyAbs1 += de;
  fTrackLAbs1 += dl;
}

inline void B4aEventAction::AddAbs_1(G4double de, G4double dl) {
  fEnergyAbs_1 += de;
  fTrackLAbs_1 += dl;
}

inline void B4aEventAction::AddGap_1(G4double de, G4double dl) {
  fEnergyGap_1 += de;
  fTrackLGap_1 += dl;
}

inline void B4aEventAction::AddAbs1_1(G4double de, G4double dl) {
  fEnergyAbs1_1 += de;
  fTrackLAbs1_1 += dl;
}

inline void B4aEventAction::AddAbs_2(G4double de, G4double dl) {
  fEnergyAbs_2 += de;
  fTrackLAbs_2 += dl;
}

inline void B4aEventAction::AddGap_2(G4double de, G4double dl) {
  fEnergyGap_2 += de;
  fTrackLGap_2 += dl;
}

inline void B4aEventAction::AddAbs1_2(G4double de, G4double dl) {
  fEnergyAbs1_2 += de;
  fTrackLAbs1_2 += dl;
}

inline void B4aEventAction::AddAbs_3(G4double de, G4double dl) {
  fEnergyAbs_3 += de;
  fTrackLAbs_3 += dl;
}

inline void B4aEventAction::AddGap_3(G4double de, G4double dl) {
  fEnergyGap_3 += de;
  fTrackLGap_3 += dl;
}

inline void B4aEventAction::AddAbs1_3(G4double de, G4double dl) {
  fEnergyAbs1_3 += de;
  fTrackLAbs1_3 += dl;
}

inline void B4aEventAction::AddAbs_4(G4double de, G4double dl) {
  fEnergyAbs_4 += de;
  fTrackLAbs_4 += dl;
}

inline void B4aEventAction::AddGap_4(G4double de, G4double dl) {
  fEnergyGap_5 += de;
  fTrackLGap_5 += dl;
}

inline void B4aEventAction::AddAbs1_4(G4double de, G4double dl) {
  fEnergyAbs1_5 += de;
  fTrackLAbs1_5 += dl;
}

inline void B4aEventAction::AddAbs_5(G4double de, G4double dl) {
  fEnergyAbs_6 += de;
  fTrackLAbs_6 += dl;
}

inline void B4aEventAction::AddGap_5(G4double de, G4double dl) {
  fEnergyGap_6 += de;
  fTrackLGap_6 += dl;
}

inline void B4aEventAction::AddAbs1_5(G4double de, G4double dl) {
  fEnergyAbs1_6 += de;
  fTrackLAbs1_6 += dl;
}

inline void B4aEventAction::AddAbs_6(G4double de, G4double dl) {
  fEnergyAbs_6 += de;
  fTrackLAbs_6 += dl;
}

inline void B4aEventAction::AddGap_6(G4double de, G4double dl) {
  fEnergyGap_6 += de;
  fTrackLGap_6 += dl;
}

inline void B4aEventAction::AddAbs1_6(G4double de, G4double dl) {
  fEnergyAbs1_6 += de;
  fTrackLAbs1_6 += dl;
}

inline void B4aEventAction::AddLorentz(G4double dpx, G4double dpy, G4double dpz, G4double dE, G4double dx, G4double dy, G4double dz) { //*******************************************
  fMomentumX += dpx; //*******************************************
  fMomentumY += dpy; //*******************************************
  fMomentumZ += dpz; //*******************************************
  fEnergyParticle += dE;
  fParticlePosX += dx;
  fParticlePosY += dy;
  fParticlePosZ += dz;
} //*******************************************

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
