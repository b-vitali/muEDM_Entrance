/// \file  EventAction.cc
/// \brief Implementation of the EventAction class

#include "TTree.h"

#include "RunAction.hh"
#include "EventAction.hh"
#include "ScintHit.hh"
#include "SiPMHit.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"

#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"
#include "G4THitsMap.hh"
#include "G4SystemOfUnits.hh"

/**
 	TODO: 	Ntuple
**/
G4bool scint = true;
G4bool sipm = true;

EventAction::EventAction(RunAction* runAction) : 
	G4UserEventAction(), fRunAction(runAction), fCollIDScint_gate(-1), fCollIDScint_telescope(-1), fCollIDSiPM_gate(-1), fCollIDSiPM_telescope(-1)
	, fEvID(-1){
		if(scint && sipm){
			tmp_scint_gate = new ScintSD("Scint_gate",2);
			tmp_sipm_gate = new SiPMSD("SiPM_gate",3);
		}
		else if(!scint && sipm){
			tmp_sipm_gate = new SiPMSD("SiPM_gate",2);
		}
		else G4cout<<"No sensitive detector aside VD"<<G4endl;
	}

EventAction::~EventAction(){}

void EventAction::BeginOfEventAction(const G4Event*){}

void EventAction::EndOfEventAction(const G4Event* event){

	G4AnalysisManager *man = G4AnalysisManager::Instance();

	// Hits collections
	G4HCofThisEvent*HCE = event->GetHCofThisEvent();
	if(!HCE) return;
	
	G4int N;
	G4bool gate = true;
	G4bool telescope = false;

	//###################################
	ScintHit* scintHit;
	if(scint && gate){
		if(fCollIDScint_gate < 0){
			G4SDManager* SDMan = G4SDManager::GetSDMpointer();
			fCollIDScint_gate = SDMan->GetCollectionID("Scint_gate/scintCollection");
		}

		if(true){
		ScintHitsCollection* ScintHitCollection_gate = (ScintHitsCollection*) (HCE->GetHC(fCollIDScint_gate));

		N = ScintHitCollection_gate->entries();
		for(int i = 0; i < N; i++){
			scintHit = (*ScintHitCollection_gate)[i];

			fEvID = event->GetEventID();

			tmp_scint_gate->FillNtupla(man, scintHit, 2);

			scintHit->Clear();
		}
		}
	}

	//###################################
	if(sipm && gate){
		if(fCollIDSiPM_gate < 0){
			G4SDManager* SDMan = G4SDManager::GetSDMpointer();
			fCollIDSiPM_gate = SDMan->GetCollectionID("SiPM_gate/sipmCollection");
		}
		SiPMHitsCollection* SiPMHitCollection_gate = (SiPMHitsCollection*) (HCE->GetHC(fCollIDSiPM_gate));


		SiPMHit* sipmHit;
		N = SiPMHitCollection_gate->entries();
		for(int i = 0; i < N; i++){
			sipmHit = (*SiPMHitCollection_gate)[i];
			G4cout<<"subhits : "<<sipmHit->GetEvent().size()<<G4endl;
			fEvID = event->GetEventID();

			if(scint)tmp_sipm_gate->FillNtupla(man, sipmHit, 3);
			else tmp_sipm_gate->FillNtupla(man, sipmHit, 2);
			sipmHit->Clear();
		}
	}

	if(fEvID % 100 == 0 || (fEvID & (fEvID - 1)) == 0 ) 
	{std::cout << "Filled second ntupla" << std::endl;
	std::cout << "Event n. " << fEvID << std::endl;}
}
