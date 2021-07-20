//ACP 2021
//implementation of Heed processes in Geant4 to code electron drift in TPC
//based off method described in https://arxiv.org/pdf/1806.05880.pdf
//and implemented in https://github.com/lennertdekeukeleere/Geant4GarfieldDegradInterface/tree/master/ALICE

#include "A2HeedModel.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4SystemOfUnits.hh"
#include "A2DetectorConstruction.hh"
#include "G4RunManager.hh"
//#include <stdio.h>
#include "A2Target.hh"
#include "A2SD.hh"
//#include "DriftLineTrajectory.hh"
#include "G4TrackingManager.hh"
#include "G4EventManager.hh"
#include "G4TransportationManager.hh"
#include "G4VVisManager.hh"


A2HeedModel::A2HeedModel(G4String modelName, G4Region* actVol, A2Target* target, A2SD* anode)
: G4VFastSimulationModel(modelName, actVol), fA2Target(target), fA2SD(anode) { 
	//empty for now
	fFakeStep = new G4Step();
	fFakePreStepPoint  = fFakeStep->GetPreStepPoint();
  	fFakePostStepPoint = fFakeStep->GetPostStepPoint();
  	fTouchableHandle   = new G4TouchableHistory();
  	fpNavigator        = new G4Navigator();
	fNaviSetup = false;
}

A2HeedModel::~A2HeedModel(){
	delete fFakeStep;
  	delete fpNavigator;
}

G4bool A2HeedModel::IsApplicable(const G4ParticleDefinition& particleType){
	G4String particleName = particleType.GetParticleName();
	if(particleName=="e-")return true; //for now set only applicable to electrons
	//G4cout<<"Heed model not applicable"<<G4endl;
	return false;
}

G4bool A2HeedModel::ModelTrigger(const G4FastTrack& fastTrack){
	G4double ekin = fastTrack.GetPrimaryTrack()->GetKineticEnergy()/keV;
	//G4String particleName = fastTrack.GetPrimaryTrack()->GetParticleDefinition()->GetParticleName();
	if (ekin <=1) { //for low energy electrons
		//G4cout<<"Heed model triggered!"<<G4endl; //debugging message
		return true;
	}
	//G4cout<<"Heed Model not triggered"<<G4endl; //debugging message
	return false;
}

//now things get complicated... need to do the physics and plot stuff
void A2HeedModel::DoIt(const G4FastTrack& fastTrack, G4FastStep& fastStep){
	//get all relevant step data from the track
	G4ThreeVector direction = fastTrack.GetPrimaryTrack()->GetMomentumDirection();
	G4ThreeVector worldPosition = fastTrack.GetPrimaryTrack()->GetPosition()/cm;
	G4double ekin = fastTrack.GetPrimaryTrack()->GetKineticEnergy()/keV;
	G4double time = fastTrack.GetPrimaryTrack()->GetGlobalTime();
	G4String particleName = fastTrack.GetPrimaryTrack()->GetParticleDefinition()->GetParticleName();
	//pass step data to transportation 
	//G4cout<<"Heed model processing..."<<G4endl;
	Transport(fastStep, fastTrack, particleName, ekin, time, worldPosition.x(), worldPosition.y(), worldPosition.z(), direction.x(), direction.y(), direction.z());
}

//take just the delta electron run instructions
void A2HeedModel::Transport(G4FastStep& fastStep,const G4FastTrack& fastTrack, G4String particleName, double ekin_keV, double t, double x_mm, double y_mm, double z_mm, double dx, double dy, double dz){
	G4cout<<"Transporting delta electron of energy "<< ekin_keV <<" keV"<<G4endl; //debugging message
	
	//calculate track parameters
	//everything here is incorrect - I'm just trying to have the simulation do something, then I'll get it to do the right thing
	//idea: first compute drift velocity. then drift the particle until a) it hits anode, or b) it hits the walls
	//find time taken to drift
	//x,y positions will be initial velocity times travel time, although some energy will be lost to the gas
	//z position will be drift velocity times travel time
	//might need to do a loop and keep processing it?
	G4double drift_vel = 200000; //cm/s... this is a lot... citing this paper here: https://journals.aps.org/pr/pdf/10.1103/PhysRev.117.470
	//probably need to implement a for loop here, or some alternate kind of stepping action
	//consider sideways motion to be negligible as a first approximation and just drift electrons to the anode	
	G4double x_pos = x_mm*mm;
	G4double y_pos = y_mm*mm;

	G4double pathLength = -115.5 - z_mm; //propogate to anode Z position
	G4double radius = sqrt(x_mm*x_mm + y_mm+y_mm); //calculate radius - will it hit the anode?
	G4double time = pathLength/drift_vel; //time taken to drift to stopping place


	
	G4double z_pos = -115.5*mm;
	//G4double z_pos = -113*mm;
	G4ThreeVector position = G4ThreeVector(x_pos,y_pos,z_pos);

	//create any necessary secondary particles
	//fastStep.CreateSecondaryTrack();

	//set final track parameters
	fastStep.KillPrimaryTrack(); //kill the step
	fastStep.SetPrimaryTrackFinalProperTime(time*s);
	//fastStep.SetPrimaryTrackFinalKineticEnergy(0); //end with no energy
	fastStep.SetPrimaryTrackPathLength(pathLength*mm); //travel calculated distance
	fastStep.SetPrimaryTrackFinalPosition(position); //set to final calculated position
	fastStep.SetTotalEnergyDeposited(ekin_keV*keV); //deposit all energy
	
	//if(radius <=50){//if it hits the anode radius
	//	G4cout<<"Anode Hit!"<<G4endl;
	//	fSensitive->Hit(fastStep);

	//} else { 
	//	G4cout<<"Electron misses anode"<<G4endl;
	//}
	ProcessHit(position,ekin_keV,time);
	
	//end
	//fastStep.UpdateStepForPostStep();
	//fastStep.DumpInfo(); //print all step information
	//fastStep.SetKineticEnergy(0);
}

void A2HeedModel::ProcessSecondaries(){
	//do nothing for now
	//later deal with secondary electrons created
}

void A2HeedModel::ProcessHit(G4ThreeVector position, G4double ekin_keV,G4double time){ //based on examples/extended/Par01
	//fill fake step to prepare for hit
	if (!fNaviSetup)
    {
      fpNavigator->SetWorldVolume(G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume());
      //G4cout<<G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume()->GetName()<<" "<<position<<G4endl;
      fpNavigator->LocateGlobalPointAndUpdateTouchableHandle(position,G4ThreeVector(0.,0.,0.),fTouchableHandle,true);
      fNaviSetup = true;
      //G4cout<<"Creating navigation setup"<<G4endl;
    } else {
      fpNavigator->
        LocateGlobalPointAndUpdateTouchableHandle(position,G4ThreeVector(0.,0.,0.),fTouchableHandle);
      //G4cout<<G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume()->GetName()<< " "<<position<<G4endl;
      //G4cout<<"Updating navigator"<<G4endl;
     }
  //fill track with attributes needed by A2Hit
  G4DynamicParticle* electron = new G4DynamicParticle(G4Electron::ElectronDefinition(),G4ThreeVector(0,0,-1),ekin_keV*1000); 
  fFakeTrack= new G4Track(electron,time,position);
  //fFakeTrack->SetStep(fFakeStep);
  //attach step to track
  //fFakeStep->InitializeStep(fFakeTrack);
  //--------------------------------------
  // Fills attribute of the G4Step needed
  // by our sensitive detector:
  //-------------------------------------
  // set touchable volume at PreStepPoint:
  fFakePreStepPoint->SetTouchableHandle(fTouchableHandle);
  // set total energy deposit:
  //fFakeStep->SetTotalEnergyDeposit(fFakeStep->GetTrack()->GetKineticEnergy());
  //G4cout<<fFakeStep->GetTotalEnergyDeposit()<<" "<<fFakeTrack->GetKineticEnergy()*keV<<G4endl;
  //fFakeStep->SetStepLength(2.*mm); //for now: checking what this does
  //fFakeStep->SetTotalEnergyDeposit(ekin_keV);
  //fFakeTrack->SetStepLength(2.*mm); //for now
  G4cout<<fFakeStep->GetStepLength()<<" "<<fFakeStep->GetTotalEnergyDeposit()<<G4endl;
 //fFakeStep->UpdateTrack(); 
	
	G4cout<<"Processing Hit..."<<G4endl;
	G4VPhysicalVolume* fCurrentVolume = fFakeStep->GetPreStepPoint()->GetPhysicalVolume();
	//G4VPhysicalVolume* fCurrentVolume = fFakePreStepPoint->GetPhysicalVolume();
	G4VSensitiveDetector* fSensitive;
	if( fCurrentVolume != 0 ) {
	fSensitive = fCurrentVolume->GetLogicalVolume()->GetSensitiveDetector();
	if( fSensitive != 0 ) {
		//troubleshooting
		if(!fSensitive->GetROgeometry()){} else{
		G4cout<<fSensitive->GetROgeometry()->GetName()<<G4endl;
		}
		if(!fSensitive->GetFilter()){} else{
		G4cout<<fSensitive->GetFilter()->GetName()<<G4endl;
		}
		//problem found: no RO geometry, therefore no hit
		//G4VReadOutGeometry *anodeRO = new G4VReadOutGeometry(); //create one???
		//anodeRO->BuildROGeometry();
		//fSensitive->SetROgeometry(anodeRO); //add one???
		//I think this should be done in the detector construction???
		fSensitive->Hit(fFakeStep);
		//G4cout<<"Hitting sensitive detector"<<G4endl;
		//fSensitive->ProcessHits();
	} else {
		//G4cout<<"Sensitive detector not found"<<G4endl;
	}
	} else {
		//G4cout<<"Current volume not found"<<G4endl;
	}
}

void A2HeedModel::GenerateDetectorResponse(){
	//do nothing for now
	//later make sure electrons are picked up by anode
}

//these are instantiated in DeltaElectronHeedModel
void A2HeedModel::ProcessEvent(){
	//do nothing
	//later make sure that things are passed back to Geant4 and relevant data is recorded
}

void A2HeedModel::Reset(){
	//do nothing
	//later reset class variables
}

