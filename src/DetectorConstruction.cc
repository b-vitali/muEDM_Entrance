/// \file  DetectorConstruction.hh
/// \brief Implementation of the class to define the experimental setup

#include "DetectorConstruction.hh"
#include <math.h>  
//#include "G4MagneticField.hh"

#include "G4AutoDelete.hh"

G4ThreadLocal 
G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = 0; 

DetectorConstruction::DetectorConstruction()
{
	fVDOn = true;

	// Scintillator dimensions
    fScintSizeX_telescope = 3*cm;
    fScintSizeY_telescope = 1*cm;
    fScintSizeZ_telescope = 20*cm;

    fScintSizeX_gate = 2*cm;
    fScintSizeY_gate = 2*cm;
    fScintSizeZ_gate = 0.1*mm;

	// VirtualDetector dimensions
    fVDSizeX = 5*cm;
    fVDSizeY = 5*cm;
    fVDSizeZ = 1*mm;

	// World dimentions
    fWorldSizeX = 3*std::max({fScintSizeX_telescope,fVDSizeX});
    fWorldSizeY = 3*std::max({fScintSizeY_telescope,fVDSizeY});
    fWorldSizeZ = 2*std::max({fScintSizeZ_telescope,fVDSizeZ});

	// At creation it calls for the function creating the materials 
	fDetectorMessenger = new DetectorMessenger(this);
	DefineMaterials();
	DefineOpticalProperties();
}

DetectorConstruction :: ~DetectorConstruction()
{
	delete fDetectorMessenger;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
	// At construction the DefineVolumes describes the geometry
	return Entrance();
}

void DetectorConstruction::DefineMaterials()
{
	G4double a; // atomic mass
	G4double z; // atomic number
	G4double density;

	//? You can take materials already defined 
	//G4NistManager* nist = G4NistManager::Instance();
	
	//? Or you can define your elements and materials
	//! Elements
	fH  = new G4Element( "H",  "H", z =  1., a =   1.01*g/mole);
	fC  = new G4Element( "C",  "C", z =  6., a =  12.01*g/mole);
	fN  = new G4Element( "N",  "N", z =  7., a =  14.01*g/mole);
	fO  = new G4Element( "O",  "O", z =  8., a =  16.00*g/mole);
	fSie= new G4Element("Si", "Si", z = 14., a = 28.0855*g/mole);
	fY  = new G4Element( "Y",  "Y", z = 39., a = 88.90585*g/mole);
	fCe = new G4Element("Ce", "Ce", z = 58., a = 140.116*g/mole);
	fLu = new G4Element("Lu", "Lu", z = 71., a = 174.967*g/mole);

	//! Materials
	// BC400
	fBC400 = new G4Material("BC400", density = 1.023*g/cm3, 2);
	fBC400->AddElement(fC, 1000);
	fBC400->AddElement(fH, 1103);

	fBC400_noscint = new G4Material("BC400", density = 1.023*g/cm3, 2);
	fBC400_noscint->AddElement(fC, 1000);
	fBC400_noscint->AddElement(fH, 1103);

	// LYSO
	fLYSO = new G4Material("LYSO", density = 7.1*g/cm3, 5);
	fLYSO->AddElement(fLu,  9);
	fLYSO->AddElement( fY, 10);
	fLYSO->AddElement(fSie, 5);
	fLYSO->AddElement( fO, 25);
	fLYSO->AddElement(fCe,  5);

	// Vacuuum
	fVacuum = new G4Material("Vacuum",z=1.,a=1.01*g/mole, 
		density = universe_mean_density, 
		kStateGas, 0.1 * kelvin, 1.e-19 * pascal);
	fVacuum_nogamma = new G4Material("Vacuum_nogamma",z=1.,a=1.01*g/mole, 
		density = universe_mean_density, 
		kStateGas, 0.1 * kelvin, 1.e-19 * pascal);

	// Air
	fAir = new G4Material("Air", density = 0.0010*g/cm3, 2);
	fAir->AddElement(fN, 70 * perCent);
	fAir->AddElement(fO, 30 * perCent);

	// Optical grease
	fOG = new G4Material("OpticalGrease",z=1.,a=1.01*g/mole, 
		     		 density = universe_mean_density, kStateGas,
				 0.1 * kelvin, 1.e-19 * pascal);

	// Silicium
	G4NistManager* NISTman = G4NistManager::Instance();
	fSi = NISTman->FindOrBuildMaterial("G4_Si");

	// Silicon resin
	fSiResin = new G4Material("SiResin",z=1.,a=1.01*g/mole, 
		     		 density = universe_mean_density, kStateGas,
				 0.1 * kelvin, 1.e-19 * pascal);

	// Assign default materials
	fScintMaterial = fBC400;
	fSiPMMaterial  = fSiResin;
}

void DetectorConstruction::DefineOpticalProperties()
{
	//? Material properties tables
	// ----------------------------------------------------------	
	//  BC400 optics
	// ----------------------------------------------------------

	std::vector<G4double> energy, scint;
	G4double tempe, tempscint;
	std::ifstream myfile;
	myfile.open("../tables/BC400_light_out.txt");
	while(true){
		myfile >> tempe >> tempscint;
		energy.push_back(1239.84197/tempe);
		scint.push_back(tempscint);
		
		if(myfile.eof()) break;
	}
	myfile.close();

	std::reverse(energy.begin(), energy.end());
  	std::reverse(scint.begin(), scint.end());

	assert(energy.size() == scint.size());
	const G4int bc400 = int(energy.size());

	G4double* BC400_Energy = new G4double[bc400];
	G4double* BC400_SCINT = new G4double[bc400];

	G4double* BC400_RIND = new G4double[bc400];
	G4double* BC400_ABSL = new G4double[bc400];
	
	for(int i = 0; i < bc400; i++){
		BC400_Energy[i] = energy.at(i)*eV;
		BC400_SCINT[i] = scint.at(i);
		BC400_RIND[i] = 1.58;
		BC400_ABSL[i] = 160*cm;
	}
	
	energy.clear();
	scint.clear();
	
	
	assert(sizeof(BC400_SCINT) == sizeof(BC400_Energy));
	
	assert(sizeof(BC400_RIND) == sizeof(BC400_Energy));
	
	assert(sizeof(BC400_ABSL) == sizeof(BC400_Energy));

	fBC400_mt = new G4MaterialPropertiesTable();
	fBC400_mt->AddProperty(       "RINDEX", BC400_Energy,  BC400_RIND, bc400);
	fBC400_mt->AddProperty(    "ABSLENGTH", BC400_Energy,  BC400_ABSL, bc400);
	fBC400_mt->AddProperty("SCINTILLATIONCOMPONENT1", BC400_Energy, BC400_SCINT, bc400);
	
	fBC400_mt->AddConstProperty("SCINTILLATIONYIELD",        11050./MeV);
	fBC400_mt->AddConstProperty(   "RESOLUTIONSCALE",                 1);
	fBC400_mt->AddConstProperty(  "SCINTILLATIONTIMECONSTANT1",            2.4*ns);
	fBC400_mt->AddConstProperty(  "SCINTILLATIONRISETIME1",   0.9*ns);
	fBC400_mt->AddConstProperty(  "SCINTILLATIONTIMECONSTANT2",            2.4*ns);
	fBC400_mt->AddConstProperty(  "SCINTILLATIONRISETIME2",   0.9*ns);
	
	fBC400->SetMaterialPropertiesTable(fBC400_mt);

	//  Set Birks Constant
	fBC400->GetIonisation()->SetBirksConstant(0.117645*mm/MeV);

	// ----------------------------------------------------------
	//  LYSO optics
	// ----------------------------------------------------------

	myfile.open("../tables/LYSO_light_out.txt");
	while(true){
		myfile >> tempe >> tempscint;
		energy.push_back(1239.84197/tempe);
		scint.push_back(tempscint);
		
		if(myfile.eof()) break;
	}
	myfile.close();

	assert(energy.size() == scint.size());
	const G4int lyso = int(energy.size());

	std::reverse(energy.begin(), energy.end());
  	std::reverse(scint.begin(), scint.end());

	G4double* LYSO_Energy = new G4double[lyso];
	G4double* LYSO_SCINT  = new G4double[lyso];

	G4double* LYSO_RIND = new G4double[lyso];
	G4double* LYSO_ABSL = new G4double[lyso];
	
	for(int i = 0; i < lyso; i++){
		LYSO_Energy[i] = energy.at(i)*eV;
		LYSO_SCINT[i] = scint.at(i);
		LYSO_RIND[i] = 1.81;
		LYSO_ABSL[i] = 20*cm;
	}
	
	energy.clear();
	scint.clear();
	
	
	assert(sizeof(LYSO_SCINT) == sizeof(LYSO_Energy));
	
	assert(sizeof(LYSO_RIND) == sizeof(LYSO_Energy));
	
	assert(sizeof(LYSO_ABSL) == sizeof(LYSO_Energy));

	fLYSO_mt = new G4MaterialPropertiesTable();
	fLYSO_mt->AddProperty(       "RINDEX", LYSO_Energy,  LYSO_RIND, lyso);
	fLYSO_mt->AddProperty(    "ABSLENGTH", LYSO_Energy,  LYSO_ABSL, lyso);
	fLYSO_mt->AddProperty("SCINTILLATIONCOMPONENT1", LYSO_Energy, LYSO_SCINT, lyso);
	
	fLYSO_mt->AddConstProperty("SCINTILLATIONYIELD",        33200./MeV);
	fLYSO_mt->AddConstProperty(   "RESOLUTIONSCALE",                 1);
	fLYSO_mt->AddConstProperty(  "SCINTILLATIONTIMECONSTANT1",            36*ns);
	fLYSO_mt->AddConstProperty(  "SCINTILLATIONTIMECONSTANT2",            36*ns);
	
	fLYSO->SetMaterialPropertiesTable(fLYSO_mt);

	//  Set Birks Constant
	//! fLYSO->GetIonisation()->SetBirksConstant(0.117645*mm/MeV);

	// ----------------------------------------------------------	
	// Vacuum & Air
	// ----------------------------------------------------------	
	fPhotonWorldPropagation = true;
	// if(fPhotonWorldPropagation){
		G4double vacuum_Energy[] = {1.5*eV, 4.*eV};
		const G4int vacnum = sizeof(vacuum_Energy) / sizeof(G4double);

		G4double vRIND = 1.;
		G4double vacuum_RIND[] = {vRIND, vRIND};
		assert(sizeof(vacuum_RIND) == sizeof(vacuum_Energy));

		G4MaterialPropertiesTable* vacuum_mt = new G4MaterialPropertiesTable();
		vacuum_mt->AddProperty("RINDEX", vacuum_Energy, vacuum_RIND, vacnum);
		fVacuum->SetMaterialPropertiesTable(vacuum_mt);
		fAir   ->SetMaterialPropertiesTable(vacuum_mt);
	// }

	// ----------------------------------------------------------	
	// Silicium
	// ----------------------------------------------------------	
	G4double Si_Energy[] = {.5*eV, 9.*eV};
	const G4int Sinum = sizeof(vacuum_Energy) / sizeof(G4double);

	G4double Si_RIND[] = {3.4, 3.4};
	assert(sizeof(Si_RIND) == sizeof(Si_Energy));

	G4MaterialPropertiesTable* Si_mt = new G4MaterialPropertiesTable();
	Si_mt->AddProperty("RINDEX", Si_Energy, Si_RIND, Sinum);
	fSi->SetMaterialPropertiesTable(Si_mt);

	// ----------------------------------------------------------	
	// Silicon resin
	// ----------------------------------------------------------	
	G4double SiRes_RIND[] = {1.41, 1.41};
	assert(sizeof(SiRes_RIND) == sizeof(Si_Energy));
	
	G4MaterialPropertiesTable* SiRes_mt = new G4MaterialPropertiesTable();
	SiRes_mt->AddProperty("RINDEX", Si_Energy, SiRes_RIND, Sinum);
	fSiResin->SetMaterialPropertiesTable(SiRes_mt);

	// ----------------------------------------------------------	
	// Optical grease
	// ----------------------------------------------------------	
	//? better if it was higher?
	G4double OG_RIND[] = {1.465, 1.465};
	assert(sizeof(OG_RIND) == sizeof(Si_Energy));
		
	G4MaterialPropertiesTable* OG_mt = new G4MaterialPropertiesTable();
	OG_mt->AddProperty("RINDEX", Si_Energy, OG_RIND, Sinum);
	fOG->SetMaterialPropertiesTable(OG_mt);
}

G4VPhysicalVolume* DetectorConstruction::Entrance()
{
	fCheckOverlaps = true;

	// World dimentions
    fWorldSizeX = 3*5*cm;
    fWorldSizeY = 3*5*cm;
    fWorldSizeZ = 3*5*cm;

	// World Solid (size) -> Logical (material) -> PVPLacement (posiz, rotaz, to interact)
	fSolidWorld	= new G4Box("World", 0.5*fWorldSizeX, 0.5*fWorldSizeY, 0.5*fWorldSizeZ);
    fLogicWorld = new G4LogicalVolume(fSolidWorld, fVacuum_nogamma, "World");
    fLogicWorld	->SetVisAttributes(G4Colour(1, 1, 1, 0.1));
	fPhysWorld	= new G4PVPlacement(0, G4ThreeVector(0,0,0), fLogicWorld, "World", 0, false, 0, fCheckOverlaps);

	// Scint sizes
    fScintSizeX_gate = 2*cm;
    fScintSizeY_gate = 2*cm;
    fScintSizeZ_gate = 0.1*mm;
	
	//? Read Solid and Phys.
	fReadSizeX = 3*mm;
	fReadSizeY = 3*mm;
	fReadSizeZ = 0.4*mm;

	double fSiPMSizeX = 2.5*mm;
	double fSiPMSizeY = 2.5*mm;
	double fSiPMSizeZ = fReadSizeZ;
	
	fSolidGrease = new G4Box("Read", 0.5*fReadSizeX,0.5*fReadSizeY, 0.5*0.5*fReadSizeZ);
    fLogicGrease = new G4LogicalVolume(fSolidGrease, fOG, "Grease");
	fLogicGrease->SetVisAttributes(G4Colour(1,0,0, 0.5));

	//? Put Grease and SiPM in Read
    G4ThreeVector Grease_pos = G4ThreeVector(0, 0, -(0.5*fReadSizeZ*0.5)); 
    G4ThreeVector SiPM_pos = G4ThreeVector(0, 0, 0.5*fReadSizeZ*0.5);
		
	// Position the Element with Scint and Read in
	// Up
	G4Rotate3D 	  rotation  = G4Rotate3D(0*90*deg, G4ThreeVector(0, 0, 1));
	G4Translate3D translate = G4Translate3D(G4ThreeVector(0, 0, 0));
	G4Transform3D transform = G4Translate3D(0,0,0)*rotation*translate;

	G4Rotate3D 	  flip_sipm  = G4Rotate3D(90*deg, G4ThreeVector(1,0,0));
	G4Rotate3D 	  piflip_sipm  = G4Rotate3D(180*deg, G4ThreeVector(1,0,0));
    
	/*
		Gate Scintillator, Element and Read
	*/

	int howmanySiPM_gate = 4;

	// Scint gate
	fSolidScint_gate 	= new G4Box("fSolidScint_gate", 0.5*fScintSizeX_gate, 0.5*fScintSizeY_gate, 0.5*fScintSizeZ_gate);
    fLogicScint_gate 	= new G4LogicalVolume(fSolidScint_gate, fScintMaterial, "fLogicScint_gate");
	fLogicScint_gate->SetVisAttributes(G4Colour(0.34, 0.57, 0.8, 0.5));

	// Element gate
	fElementSizeX_gate = fScintSizeX_gate + 2*fReadSizeZ;
	fElementSizeY_gate = fScintSizeY_gate + 2*fReadSizeZ;
	fElementSizeZ_gate = std::max(fScintSizeZ_gate, fReadSizeX);

	fSolidElement_gate = new G4Box("Element_gate", 0.5*(fElementSizeX_gate),0.5*(fElementSizeY_gate), 0.5*(fElementSizeZ_gate));
    fLogicElement_gate = new G4LogicalVolume(fSolidElement_gate, fVacuum, "fSolidElement_gate");
	fLogicElement_gate->SetVisAttributes(G4Colour(0, 1, 0, 0.2));

	// Read gate
	fSolidRead_gate	= new G4Box("fSolidRead_gate", 0.5*fReadSizeX,0.5*fReadSizeY, 0.5*fReadSizeZ);
    fLogicRead_gate = new G4LogicalVolume(fSolidRead_gate, fSi, "fSolidRead_gate");
	fLogicRead_gate->SetVisAttributes(G4Colour(1.0, 0.0, 1.0, 0.3));

	// SiPM
	fSolidSiPM_gate	= new G4Box("SiPM_gate", 0.5*fSiPMSizeX, 0.5*fSiPMSizeY, 0.5*0.5*fSiPMSizeZ);
    fLogicSiPM_gate = new G4LogicalVolume(fSolidSiPM_gate, fSiPMMaterial, "SiPM_gate");
    fLogicSiPM_gate	->SetVisAttributes(G4Colour(0,0,1, 0.5));

	fPhysGrease	= new G4PVPlacement(0, Grease_pos, fLogicGrease, "Grease", fLogicRead_gate, false, fCheckOverlaps);
	fPhysSiPM_gate		= new G4PVPlacement(0, SiPM_pos, fLogicSiPM_gate, "SiPM_gate", fLogicRead_gate, false, fCheckOverlaps);

	// Position the Element and the Scint and Read inside
	fPhysElement_gate  = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), fLogicElement_gate, "Element_gate", fLogicWorld, false, fCheckOverlaps);
	fPhysScint_gate	   = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), fLogicScint_gate, "Scint_gate", fLogicElement_gate, false, fCheckOverlaps);
	
	flip_sipm  = G4Rotate3D(90*deg, G4ThreeVector(0,1,0));

	double step = fReadSizeX + (fScintSizeX_gate - howmanySiPM_gate*fReadSizeX) / howmanySiPM_gate;

	for(int i = 0; i<4; i++){
	for(int j=0; j<howmanySiPM_gate; j += 1){
		rotation  = G4Rotate3D(i*90*deg, G4ThreeVector(0, 0, 1));
		if(howmanySiPM_gate%2==1){
			translate = G4Translate3D(G4ThreeVector(0.5*fScintSizeX_gate+0.5*fReadSizeZ, (j-howmanySiPM_gate/2)*step, 0));
		}
		else{
			translate = G4Translate3D(G4ThreeVector(0.5*fScintSizeX_gate+0.5*fReadSizeZ, (j-howmanySiPM_gate/2+0.5)*step, 0));
		}
		transform = rotation*translate*flip_sipm;
		fPhysRead_gate	   = new G4PVPlacement(transform, fLogicRead_gate, "Read_gate", fLogicElement_gate, true, i*10+j, fCheckOverlaps);
	}
	}
	if(fVDOn){
		// VirtualDetector Solid (size) -> Logical (material) -> PVPLacement (posiz, rotaz, to interact)
		fVDSizeZ = fElementSizeZ_gate*0.5-0.5*mm; //Updated to fill 
		fSolidVD	= new G4Box("VD", 0.5*fVDSizeX, 0.5*fVDSizeY, 0.5*fVDSizeZ);
    	fLogicVD 	= new G4LogicalVolume(fSolidVD, fVacuum_nogamma, "VD");
    	fPhysVD		= new G4PVPlacement(0, G4ThreeVector(0., 0., fElementSizeZ_gate*0.5 + fVDSizeZ*0.5), fLogicVD, "VD", fLogicWorld, true, 0, fCheckOverlaps);
		fPhysVD		= new G4PVPlacement(0, G4ThreeVector(0., 0., -fElementSizeZ_gate*0.5 - fVDSizeZ*0.5), fLogicVD, "VD", fLogicWorld, true, 2, fCheckOverlaps);
    	fLogicVD	->SetVisAttributes(G4Colour(0.8, 0.34, 0.68, 0.1));
	}
    return fPhysWorld;
}


void DetectorConstruction::ConstructSDandField()
{
	auto sdManager = G4SDManager::GetSDMpointer();
   	G4String SDname;

	if(fLogicScint_gate){
		ScintSD* scint_SD_gate = new ScintSD(SDname="Scint_gate");
  		sdManager->AddNewDetector(scint_SD_gate);
		fLogicScint_gate->SetSensitiveDetector(scint_SD_gate);
	}	
	

	if(fLogicSiPM_gate){
		SiPMSD * SiPM_SD_gate = new SiPMSD("SiPM_gate");
  		sdManager->AddNewDetector(SiPM_SD_gate);
		fLogicSiPM_gate->SetSensitiveDetector(SiPM_SD_gate);
	}

	if(fLogicVD){
		// Create the Sensitive Detector defined in VirtualDetectorSD 
		VirtualDetectorSD * VD_SD = new VirtualDetectorSD("VirtualDetector");
		// Assign the SD to the logial volume
		fLogicVD->SetSensitiveDetector(VD_SD);
	}
}

/*
	From here functions used through the Messenger to modify the detector
*/
void DetectorConstruction :: SetScintSize(G4double size){
	fScintSizeZ_gate = size;
	fSolidScint_gate->SetZHalfLength(size*0.5);
	G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void DetectorConstruction :: SetScintSize(G4ThreeVector size){
	fScintSizeX_gate = size.getX();
	fScintSizeY_gate = size.getY();
	fScintSizeZ_gate = size.getZ();

	fSolidScint_gate->SetXHalfLength(size.getX()*0.5);
	fSolidScint_gate->SetYHalfLength(size.getY()*0.5);
	fSolidScint_gate->SetZHalfLength(size.getZ()*0.5);

	// Element gate
	fSolidElement_gate->SetXHalfLength((size.getX()+1*mm)*0.5);
	fSolidElement_gate->SetYHalfLength((size.getY()+1*mm)*0.5);
	fSolidElement_gate->SetZHalfLength((size.getZ()+1*mm)*0.5);

	// Read gate
	fSolidRead_gate->SetXHalfLength((0.5*mm)*0.5);
	fSolidRead_gate->SetYHalfLength((size.getY())*0.5);
	fSolidRead_gate->SetZHalfLength((size.getZ())*0.5);

	G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void DetectorConstruction::SetScintMaterial(G4String name){
	if(name == "BC400") fScintMaterial = fBC400;
	else if(name == "LYSO") fScintMaterial = fLYSO;
	G4RunManager::GetRunManager()->ReinitializeGeometry();
}