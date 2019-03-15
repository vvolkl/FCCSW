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
// $Id: DreamDetectorConstruction.cc 87359 2014-12-01 16:04:27Z gcosmo $
//
/// \file DreamDetectorConstruction.cc
/// \brief Implementation of the DreamDetectorConstruction class

#include "SimG4Dream/DreamDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
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

#include "G4OpBoundaryProcess.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"

#include "G4Trap.hh"
#include "G4Trd.hh"

#include "SimG4Dream/dimensionB.hh"
#include "SimG4Dream/dimensionE.hh"

#include <string>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal
G4GlobalMagFieldMessenger* DreamDetectorConstruction::fMagFieldMessenger = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DreamDetectorConstruction::DreamDetectorConstruction()
: G4VUserDetectorConstruction(),
modulePV(0),
fCheckOverlaps(true)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DreamDetectorConstruction::~DreamDetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DreamDetectorConstruction::Construct()
{
    // Define materials
	DefineMaterials();

    // Define volumes
	return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void DreamDetectorConstruction::DefineMaterials()
{
    // Copper material defined using NIST Manager
    // I use Cu as default absorber material but you can switch to lead
	G4NistManager* CunistManager = G4NistManager::Instance();
	CunistManager->FindOrBuildMaterial("G4_Cu");

    // Lead material defined using NIST Manager
    // G4NistManager* PbnistManager = G4NistManager::Instance();
    // PbnistManager->FindOrBuildMaterial("G4_Pb");

    // Polystyrene material defined using NIST Manager
    // I use this material for the core of plastic scintillating fibers
    // cannot find any G4_Polystyrene, I build it later
    //G4NistManager* PynistManager = G4NistManager::Instance();
    //PynistManager->FindOrBuildMaterial("G4_Polystyrene");

    // PMMA material, there's no default G4_PMMA, I build it (C502H8)
    G4String name, symbol;    // a=mass of a mole;
    G4double a, z;            // z=mean number of protons;
    
    // create elements
    a = 1.01*g/mole;
    G4Element* elH  = new G4Element(name="Hydrogen",symbol="H" , z= 1., a); //Hidrogen
    
    a = 12.01*g/mole;
    G4Element* elC  = new G4Element(name="Carbon"  ,symbol="C" , z= 6., a); //Carbon
    
    a = 16.00*g/mole;
    G4Element* elO  = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a); //Oxygen
    
    a = 28.09*g/mole;
    G4Element* elSi = new G4Element(name="Silicon", symbol="Si", z=14., a); //Silicon
    
    a = 18.9984*g/mole;
    G4Element* elF  = new G4Element("Fluorine",symbol="F" , z= 9., a); //Fluorine
    
    a = 63.546*g/mole;
    G4Element* elCu = new G4Element("Copper", symbol="Cu", z=29., a); //Copper
    
    a = 65.38*g/mole;
    G4Element* elZn = new G4Element("Zinc", symbol="Zn", z=30., a); //Zinc
    
    // create PMMA
    G4Material* PMMA = new G4Material("PMMA", 1.19*g/cm3, 3); //name, density and number of elements
    PMMA -> AddElement(elC, 5);
    PMMA -> AddElement(elO, 2);
    PMMA -> AddElement(elH, 8); //PMMA building complete
    
    // create Polystyrene (C5H5)
    G4Material* Polystyrene = new G4Material("Polystyrene", 1.05*g/cm3, 2);
    Polystyrene -> AddElement(elC, 8);
    Polystyrene -> AddElement(elH, 8); //Polystyrene building complete
    
    // create Fluorinated Polymer (C2F2)
    // I use it for the cladding of the Cherenkov fibers
    G4Material* fluorinatedPolymer =
    new G4Material("Fluorinated_Polymer", 1.43*g/cm3, 2);
    fluorinatedPolymer->AddElement(elC,2);
    fluorinatedPolymer->AddElement(elF,2);
    //fluorinatedPolymer->AddElement(H,2); //Fluorinated Polymer building complete
    
    // create Glass (SiO2)
    G4Material* Glass = new G4Material("Glass", 2.4*g/cm3, 2);
    Glass -> AddElement(elSi, 1);
    Glass -> AddElement(elO, 2); //Glass building complete
    
    // Vacuum material defined using NIST Manager
    G4NistManager* VanistManager = G4NistManager::Instance();
    VanistManager->FindOrBuildMaterial("G4_Galactic");
    
    // Silicon material defined using NIST Manager
    G4NistManager* SinistManager = G4NistManager::Instance();
    SinistManager->FindOrBuildMaterial("G4_Si");
    
    // create Cu260 (Brass)
    // I use it for the absorber of the real small beam tested module
    double density = 8.53*g/cm3;
    int ncomponentsbrass = 2;
    G4Material* Cu260 = new G4Material(name="Brass", density, ncomponentsbrass);
    Cu260->AddElement(elCu, 70*perCent);
    Cu260->AddElement(elZn, 30*perCent);
    
    // Air material defined using NIST Manager
    // You can use Air instead of vacuum
    G4NistManager* AirnistManager = G4NistManager::Instance();
    AirnistManager->FindOrBuildMaterial("G4_AIR");
    
    // Print materials
    // I don't want to print materials all the times,
    // if you want uncomment it
    //G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DreamDetectorConstruction::DefineVolumes()
{

    // Geometry parameters of the world, world is a box
	G4double worldX = 10*m;
	G4double worldY = 10*m;
	G4double worldZ = 10*m;

    // Get materials for vacuum, absorber, scintillating and cherenkov fibers, SiPM
    G4Material* defaultMaterial = G4Material::GetMaterial("G4_Galactic"); // G4_AIR or G4_Galactic
    
    // Building the calorimeter
    
    // Here I build the world
    
    G4VSolid* worldS
    = new G4Box("World",                        // its name
                worldX/2, worldY/2, worldZ/2); // its size
    
    G4LogicalVolume* worldLV
    = new G4LogicalVolume(
                          worldS,           // its solid
                          defaultMaterial,  // its material (Galactic or Air)
                          "World");         // its name
    
    // I set the world as invisible
    worldLV->SetVisAttributes(G4VisAttributes::Invisible);
    
    G4VPhysicalVolume* worldPV
    = new G4PVPlacement(
                        0,                // no rotation
                        G4ThreeVector(),  // at (0,0,0)
                        worldLV,          // its logical volume
                        "World",          // its name
                        0,                // its mother  volume
                        false,            // no boolean operation
                        0,                // copy number
                        fCheckOverlaps);  // checking overlaps
    
    G4NistManager* nistManager = G4NistManager::Instance();
    
    G4String symbol;             //a=mass of a mole;
    G4double a, z, density;      //z=mean number of protons;
    G4int iz, n;                 //iz=number of protons  in an isotope;
    // n=number of nucleons in an isotope;
    
    G4int ncomponents, natoms;
    G4double abundance, fractionmass;
    G4Material* cu  =new G4Material("Copper"   , z=29., a=63.546*g/mole, density=8.96*g/cm3);
    G4Element* H  = nistManager->FindOrBuildElement(1);
    G4Element* C  = nistManager->FindOrBuildElement(6);
    G4Element* N  = nistManager->FindOrBuildElement(7);
    G4Element* O  = nistManager->FindOrBuildElement(8);
    G4Element* F  = nistManager->FindOrBuildElement(9);
    G4Element* Si = nistManager->FindOrBuildElement(14);
    
    //--- for PMT Cathod ---
    G4Material* Al = new G4Material("Aluminium", z=13., a=26.98*g/mole, density=2.700*g/cm3);
    
    //--- for PMT Glass ---
    G4Material* Glass = new G4Material("Glass", density=1.032*g/cm3,2);
    Glass->AddElement(C,91.533*perCent);
    Glass->AddElement(H,8.467*perCent);
    
    ///--- for scintillation fiber core ---
    G4Material* polystyrene =
    new G4Material("Polystyrene",density= 1.05*g/cm3, ncomponents=2);
    polystyrene->AddElement(C, natoms=8);
    polystyrene->AddElement(H, natoms=8);
    
    ///--- for cladding (scintillation fibers) ---
    G4Material* pmma_clad =
    new G4Material("PMMA_Clad",density= 1.19*g/cm3, ncomponents=3);
    pmma_clad->AddElement(C, natoms=5);
    pmma_clad->AddElement(H, natoms=8);
    pmma_clad->AddElement(O, natoms=2);
    
    ///--- for Cerenkov fiber core ---
    G4Material* pmma =
    new G4Material("PMMA",density= 1.19*g/cm3, ncomponents=3);
    pmma->AddElement(C, natoms=5);
    pmma->AddElement(H, natoms=8);
    pmma->AddElement(O, natoms=2);
    
    ///--- for cladding (Cerenkov fibers) ---
    G4Material* fluorinatedPolymer =
    new G4Material("Fluorinated_Polymer", density= 1.43*g/cm3, ncomponents=2);
    fluorinatedPolymer->AddElement(C,2);
    fluorinatedPolymer->AddElement(F,2);
    
    G4Material* Air = nistManager->FindOrBuildMaterial("G4_AIR",false);
    
    ///--- Material property tables for fiber materials ---
    G4MaterialPropertiesTable* mpAir;
    G4MaterialPropertiesTable* mpPS;
    G4MaterialPropertiesTable* mpPMMA;
    G4MaterialPropertiesTable* mpFS;
    G4MaterialPropertiesTable* mpGlass;
    G4MaterialPropertiesTable* mpPMTPC;
    
    //
    //--- Generate and add material properties table ---
    //
    G4double PhotonEnergy[] = {2.00*eV,2.03*eV,2.06*eV,2.09*eV,2.12*eV,
    	2.15*eV,2.18*eV,2.21*eV,2.24*eV,2.27*eV,
    	2.30*eV,2.33*eV,2.36*eV,2.39*eV,2.42*eV,
    	2.45*eV,2.48*eV,2.51*eV,2.54*eV,2.57*eV,
    	2.60*eV,2.63*eV,2.66*eV,2.69*eV,2.72*eV,
    	2.75*eV,2.78*eV,2.81*eV,2.84*eV,2.87*eV,
    	2.90*eV,2.93*eV,2.96*eV,2.99*eV,3.02*eV,
    	3.05*eV,3.08*eV,3.11*eV,3.14*eV,3.17*eV,
    	3.20*eV,3.23*eV,3.26*eV,3.29*eV,3.32*eV,
    	3.35*eV,3.38*eV,3.41*eV,3.44*eV,3.47*eV};

    	const G4int nEntries = sizeof(PhotonEnergy) / sizeof(G4double);
    //--- PMMA ---
    	G4double RefractiveIndex_PMMA[nEntries] =
    	{
    		1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
    		1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
    		1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
    		1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
    		1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49
    	};
    	mpPMMA = new G4MaterialPropertiesTable();
    	mpPMMA->AddProperty("RINDEX",PhotonEnergy,RefractiveIndex_PMMA,nEntries);
    	pmma->SetMaterialPropertiesTable(mpPMMA);

    //--- Fluorinated Polymer (FS) ---
    	G4double RefractiveIndex_FluorinatedPolymer[nEntries] =
    	{
    		1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42,
    		1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42,
    		1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42,
    		1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42,
    		1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42
    	};
    	mpFS = new G4MaterialPropertiesTable();
    	mpFS->AddProperty("RINDEX",PhotonEnergy,RefractiveIndex_FluorinatedPolymer,nEntries);
    	fluorinatedPolymer->SetMaterialPropertiesTable(mpFS);

    //
    //Glass
    //
    	G4double RefractiveIndex_Glass[nEntries] =
    	{  1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
    		1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
    		1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
    		1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
    		1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49
    	};
    	mpGlass = new G4MaterialPropertiesTable();
    	mpGlass->AddProperty("RINDEX",PhotonEnergy,RefractiveIndex_Glass,nEntries);
    //mpGlass->AddProperty("ABSLENGTH",PhotonEnergy,Glass_AbsLength,nEntries);
    //Glass->SetMaterialPropertiesTable(mpGlass);

    	G4double RefractiveIndex_Air[nEntries] =
    	{
    		1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
    		1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
    		1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
    		1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
    		1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00
    	};

    	mpAir = new G4MaterialPropertiesTable();
    	mpAir->AddProperty("RINDEX", PhotonEnergy, RefractiveIndex_Air, nEntries);
    	Air->SetMaterialPropertiesTable(mpAir);
    //default materials of the World

    //---Materials for Cerenkov fiber---
    	clad_C_Material = fluorinatedPolymer;
    	core_C_Material = pmma;

    //---Materials for Scintillation fiber---
    	clad_S_Material = pmma_clad;
    	core_S_Material = polystyrene;

    //--Material for PMT glass---
    	Glass_Material = Glass;

    //--- Material for PMT Photocathod ---
    	PMTPC_Material = Al;

    //--- Photocathod property ---
    	G4double p_mppc[2] = {2.00*eV, 3.47*eV};
    	G4double refl_mppc[2] = {0.0, 0.0};
    G4double effi_mppc[2] = {0.11, 0.11}; // mimic Quantum Efficiency
    G4double photocath_ReR[] = {1.92, 1.92};
    G4double photocath_ImR[] = {1.69, 1.69};
    
    mpPMTPC = new G4MaterialPropertiesTable();
    mpPMTPC->AddProperty("REFLECTIVITY",p_mppc,refl_mppc,2);
    mpPMTPC->AddProperty("EFFICIENCY",p_mppc,effi_mppc,2);
    
    G4OpticalSurface* photocath_opsurf = new G4OpticalSurface("photocath_opsurf",glisur,polished,dielectric_metal);
    photocath_opsurf->SetMaterialPropertiesTable(mpPMTPC);
    
    //Calorimeter parameters
    innerR = 2500; //inner radius /1800
    tower_height = 2000; //tower height 2500
    NbOfBarrel = 40; //(it was 52 before) number of towers in barrel right (left)
    //NbOfEndcap = 47;
    NbOfEndcap = 46; //number of towers in endcap
    NbOfZRot = 16; //283 //number of Z to round around the center
    //PMTT = 1*mm;
    PMTT = 0*mm; //? giulia, prima =1;
    
    
    fulltheta = 0; //?
    
    //2*pi/number of tower to complete a rotation around the center
    phi_unit = 2*M_PI/(G4double)NbOfZRot;
    
    //Parameters for fibers
    clad_C_rMin = 0.49*mm; //cladding cherenkov minimum radius
    clad_C_rMax = 0.50*mm; //cladding cherenkov max radius
    clad_C_Dz   = 2.5*m;   //cladding cherenkov lenght
    clad_C_Sphi = 0.;      //cladding cherenkov min rotation
    clad_C_Dphi = 2.*M_PI; //cladding chrenkov max rotation
    
    core_C_rMin = 0.*mm;
    core_C_rMax = 0.49*mm;
    core_C_Dz   = 2.5*m;
    core_C_Sphi = 0.;
    core_C_Dphi = 2.*M_PI;
    
    clad_S_rMin = 0.485*mm;
    clad_S_rMax = 0.50*mm;
    clad_S_Dz   = 2.5*m;
    clad_S_Sphi = 0.;
    clad_S_Dphi = 2.*M_PI;
    
    core_S_rMin = 0.*mm;
    core_S_rMax = 0.485*mm;
    core_S_Dz   = 2.5*m;
    core_S_Sphi = 0.;
    core_S_Dphi = 2.*M_PI;
    
    theta_unit=0; //?
    deltatheta=0; //?
    thetaofcenter=0; //?
    
    //creating fibers solids
    G4cout << "r_clad= " << clad_C_rMax << " r_coreC=" << core_C_rMax << " r_coreS=" << core_S_rMax << G4endl;
    fiber = new G4Tubs("fiber",0,clad_C_rMax,tower_height/2.,0*deg,360.*deg);// S is the same
    fiberC = new G4Tubs("fiberC",0,core_C_rMax,tower_height/2.,0*deg,360.*deg);
    fiberS = new G4Tubs("fiberS",0,core_S_rMax,tower_height/2.,0*deg,360.*deg);
    
    //vector for logical volumes of fibers
    //G4LogicalVolume* fiberCLog[2500];
    //G4LogicalVolume* fiberSLog[2500];
    
    // Prepare for logical volume of fiber tower_height=2000 (before it was 2000)
    for(int length=1;length<=tower_height;length++){ //from 1 to 20000
        double half=0.5*length; //half to build objects with proper dimensions
        char name[80];
        sprintf(name,"fiber%d",length);
        fiber = new G4Tubs(name,0,clad_C_rMax,half,0*deg,360.*deg); //creating fibers G4Tubs
        sprintf(name,"fiberC%d",length);
        fiberC = new G4Tubs(name,0,core_C_rMax,half,0*deg,360.*deg);
        sprintf(name,"fiberS%d",length);
        fiberS = new G4Tubs(name,0,core_S_rMax,half,0*deg,360.*deg);
        visAttrC = new G4VisAttributes(G4Colour(0.,0.,1.0));
        visAttrC->SetVisibility(true);
        visAttrC->SetDaughtersInvisible(true);
        visAttrC->SetForceWireframe(true);
        visAttrC->SetForceSolid(true);
        visAttrS = new G4VisAttributes(G4Colour(1.,0.,0.));
        visAttrS->SetVisibility(true);
        visAttrS->SetDaughtersInvisible(true);
        visAttrS->SetForceWireframe(true);
        visAttrS->SetForceSolid(true);
        
        fiberCLog[length] = new G4LogicalVolume(fiber,clad_C_Material,"fiberCladC");
        fiberSLog[length] = new G4LogicalVolume(fiber,clad_S_Material,"fiberCladS");
        fiberCLog[length]->SetVisAttributes(visAttrC);
        fiberSLog[length]->SetVisAttributes(visAttrS);
        G4LogicalVolume* fiberCoreCLog = new G4LogicalVolume(fiberC,core_C_Material,"fiberCoreC");
        G4LogicalVolume* fiberCoreSLog = new G4LogicalVolume(fiberS,core_S_Material,"fiberCoreS");
        // new G4PVPlacement(0,G4ThreeVector(0,0,0),fiberCoreCLog,"fiberCoreCherePhys",fiberCLog[length],false,0,fCheckOverlaps);
        //new G4PVPlacement(0,G4ThreeVector(0,0,0),fiberCoreSLog,"fiberCoreScintPhys",fiberSLog[length],false,0,fCheckOverlaps);
        /*if(sd){
         fiberCoreCLog->SetSensitiveDetector(sd);
         fiberCoreSLog->SetSensitiveDetector(sd);
         }*/
    }
    //Final logical volumes of fibers
    fiberCladCLog = fiberCLog[2500];
    fiberCladSLog = fiberSLog[2500];
    
    //Counter on number of volumes
    G4int volnum=0;
    
    //G4double deltatheta_barrel[] = {0.02222,0.02220,0.02217,0.02214,0.02209,0.02203,0.02196,0.02188,0.02179,0.02169,0.02158,0.02146,0.02133,0.02119,0.02105,0.02089,0.02073,0.02056,0.02039,0.02020,0.02002,0.01982,0.01962,0.01941,0.01920,0.01898,0.01876,0.01854,0.01831,0.01808,0.01785,0.01761,0.01738,0.01714,0.01689,0.01665,0.01641,0.01616,0.01592,0.01567,0.01543,0.01518,0.01494,0.01470,0.01445,0.01421,0.01397,0.01373,0.01350,0.01326,0.01303,0.01280};
    G4double deltatheta_barrel[] = {0.02222,0.02220,0.02217,0.02214,0.02209,0.02203,0.02196,0.02188,0.02179,0.02169,0.02158,0.02146,0.02133,0.02119,0.02105,0.02089,0.02073,0.02056,0.02039,0.02020,0.02002,0.01982,0.01962,0.01941,0.01920,0.01898,0.01876,0.01854,0.01831,0.01808,0.01785,0.01761,0.01738,0.01714,0.01689,0.01665,0.01641,0.01616,0.01592,0.01567};
    G4double deltatheta_endcap = deltatheta_barrel[NbOfBarrel-1];
    
    double thetaB = 0;
    for(int i=0;i<NbOfBarrel;i++) thetaB += deltatheta_barrel[i];
    	double thetaE = 0;
    for(int i=0;i<NbOfEndcap;i++) thetaE += deltatheta_endcap;
    	G4cout<<"theta"<< thetaB<< " "<<thetaE<<G4endl;
    double length = tower_height; //(was tower height + 100)length of physical volumes
    double innerR_Endcap = 3550.0; //2828.4;//3550.0 //3125.83
    double l1 = innerR_Endcap;
    //double l1 = innerR_Endcap/cos(thetaB);
    double l2 = (l1+length)/cos(0.5*thetaE);
    //double l = l1-l1*cos(0.5*thetaE)+length;
    double lE = (l2-l1)*cos(0.5*thetaE);
    double rEC = 0.5*(l1+l2)*cos(0.5*thetaE);
    G4Trd* phiBarrel = new G4Trd("phiBarrel",(innerR)*tan(0.5*phi_unit),(innerR+length)*tan(0.5*phi_unit),(innerR)*tan(thetaB),(innerR+length)*tan(thetaB),0.5*length);
    G4Trd* phiBarrel2 = new G4Trd("phiBarrel2",(innerR)*tan(0.5*phi_unit),(innerR+length)*tan(0.5*phi_unit),(innerR)*tan(thetaB),(innerR+length)*tan(thetaB),0.5*length);
    //G4Cons* phiDivision = new G4Cons("phiDivision",innerR-100,innerR+tower_height+100,innerR-100,innerR+tower_height+100,5000,-0.5*phi_unit,phi_unit);
    G4LogicalVolume* phiBLog = new G4LogicalVolume(phiBarrel,Air,"phiBLog");
    G4LogicalVolume* phiB2Log = new G4LogicalVolume(phiBarrel2,Air,"phiB2Log");
    //G4Trap* phiER = new G4Trap("phiER",0.5*lE,0,0,
    //         l1*sin(0.5*thetaE),l1*cos(thetaB+thetaE)*tan(0.5*phi_unit),l1*cos(thetaB)*tan(0.5*phi_unit),0,
    //         l2*sin(0.5*thetaE),l2*cos(thetaB+thetaE)*tan(0.5*phi_unit),l2*cos(thetaB)*tan(0.5*phi_unit),0);
    G4Trap* phiER = new G4Trap("phiER",0.5*lE,0,0,
    	l1*sin(0.5*thetaE),l1*cos(thetaB+thetaE)*tan(0.5*phi_unit),l1*cos(thetaB)*tan(0.5*phi_unit),0,
    	l2*sin(0.5*thetaE),l2*cos(thetaB+thetaE)*tan(0.5*phi_unit),l2*cos(thetaB)*tan(0.5*phi_unit),0);
    G4Trap* phiEL = new G4Trap("phiEL",0.5*(l2-l1)*cos(0.5*thetaE),0,0,
    	l1*sin(0.5*thetaE),l1*cos(thetaB)*tan(0.5*phi_unit),l1*cos(thetaB+thetaE)*tan(0.5*phi_unit),0,
    	l2*sin(0.5*thetaE),l2*cos(thetaB)*tan(0.5*phi_unit),l2*cos(thetaB+thetaE)*tan(0.5*phi_unit),0);
    G4LogicalVolume* phiERLog = new G4LogicalVolume(phiER,Air,"phiERLog");
    G4LogicalVolume* phiELLog = new G4LogicalVolume(phiEL,Air,"phiELLog");
    
    for(int j=0;j<NbOfZRot;j++){ //j<NbOfZRot
        //for(int j=0;j<1;j++){
        //if(j>NbOfZRot/2)break; //half for draw viewer
    	G4RotationMatrix* rmB = new G4RotationMatrix();
    	rmB->rotateZ(M_PI/2.);
    	rmB->rotateZ(-j*phi_unit);
    	rmB->rotateX(M_PI/2.);
        //place physical spacing of barrel
    	new G4PVPlacement(rmB,G4ThreeVector((innerR+0.5*length)*cos(j*phi_unit),(innerR+0.5*length)*sin(j*phi_unit),0),phiBLog,"phiDivPhys",worldLV,false,j,false);
        /*G4RotationMatrix* rmB2 = new G4RotationMatrix();
         rmB2->rotateZ(-M_PI/2.);
         rmB2->rotateZ(-j*phi_unit);
         rmB2->rotateX(0);
         rmB2->rotateY(0);*/
        //new G4PVPlacement(rmB2,G4ThreeVector((innerR+0.5*length)*cos(j*phi_unit),(innerR+0.5*length)*sin(j*phi_unit)-3515,+3515),phiB2Log,"phiDivPhys",worldLV,false,j,false);
        /*G4RotationMatrix* rmB3 = new G4RotationMatrix();
         rmB3->rotateZ(M_PI/2.);
         rmB3->rotateZ(-j*phi_unit);
         rmB3->rotateX(0);
         rmB3->rotateY(0);*/
        //new G4PVPlacement(rmB3,G4ThreeVector((innerR+0.5*length)*cos(j*phi_unit),(innerR+0.5*length)*sin(j*phi_unit)-3515,+3515),phiB2Log,"phiDivPhys",worldLV,false,j,false);
    	G4RotationMatrix* rmER = new G4RotationMatrix();
    	rmER->rotateZ(M_PI/2.);
    	rmER->rotateZ(-j*phi_unit);
    	rmER->rotateX(0.5*M_PI-thetaB-0.5*thetaE);
        /*G4RotationMatrix* rmER = new G4RotationMatrix();
         rmER->rotateZ(0);
         rmER->rotateZ(0);
         rmER->rotateX(0);}*/
        //endcap right
        //new G4PVPlacement(rmER,G4ThreeVector(rEC*cos(thetaB+0.5*thetaE)*cos(j*phi_unit),rEC*cos(thetaB+0.5*thetaE)*sin(j*phi_unit),rEC*sin(thetaB+0.5*thetaE)),
        //    phiERLog,"phiERPhys",worldLV,false,j,fCheckOverlaps);
    	new G4PVPlacement(rmER,G4ThreeVector(rEC*cos(thetaB+0.5*thetaE)*cos(j*phi_unit),rEC*cos(thetaB+0.5*thetaE)*sin(j*phi_unit),rEC*sin(thetaB+0.5*thetaE)),
    		phiERLog,"phiERPhys",worldLV,false,j,fCheckOverlaps);
    	G4RotationMatrix* rmEL = new G4RotationMatrix();
    	rmEL->rotateZ(M_PI/2.);
    	rmEL->rotateZ(-j*phi_unit);
    	rmEL->rotateX(0.5*M_PI+thetaB+0.5*thetaE);
    	new G4PVPlacement(rmEL,G4ThreeVector(rEC*cos(thetaB+0.5*thetaE)*cos(j*phi_unit),rEC*cos(thetaB+0.5*thetaE)*sin(j*phi_unit),-rEC*sin(thetaB+0.5*thetaE)),
    		phiELLog,"phiELPhys",worldLV,false,j,fCheckOverlaps);
    }
    G4VisAttributes* phiVisAttr = new G4VisAttributes(G4Colour(0.8,0.8,0.8,0.3));
    phiVisAttr->SetVisibility(true);
    phiBLog->SetVisAttributes(phiVisAttr);
    phiERLog->SetVisAttributes(phiVisAttr);
    phiELLog->SetVisAttributes(phiVisAttr);
    phiB2Log->SetVisAttributes(phiVisAttr);
    
    G4VisAttributes* towerVisAttr = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
    towerVisAttr->SetVisibility(true);
    towerVisAttr->SetDaughtersInvisible(false);
    towerVisAttr->SetForceWireframe(true);
    
    G4VisAttributes* PMTVisAttr = new G4VisAttributes(G4Colour(0.3,0.6,0.0));
    PMTVisAttr->SetVisibility(true);
    PMTVisAttr->SetDaughtersInvisible(false);
    PMTVisAttr->SetForceWireframe(true);
    
    //BARREL
    dimB = new dimensionB();
    dimB->SetInnerR(innerR);
    dimB->SetTower_height(tower_height);
    dimB->SetNumZRot(NbOfZRot);
    dimB->SetPMTT(PMTT);
    
    // barrel R
    G4cout << "Barrel R..." << G4endl;
    dimB->Rbool(1);
    
    for(int i=0;i<NbOfBarrel;i++){    //i<NbOfBarrel
        //for(int i=0;i<1;i++){
        //if(volnum==0){//
    	thetaofcenter=fulltheta+deltatheta_barrel[i]/2.;
        //for(int k=0;k<NbOfBarrel;k++) thetaofcenter += deltatheta_barrel[k];
        //G4cout << i << " : " << thetaofcenter << G4endl;
    	dimB->SetDeltaTheta(deltatheta_barrel[i]);
    	dimB->SetThetaOfCenter(thetaofcenter);
    	dimB->CalBasic();
    	dimB->Getpt(pt);
        //G4cout << i << " : " << pt[0] << pt[1] << pt[2] << pt[3] << G4endl;
    	sprintf(name,"tower%d",volnum);
    	tower = new G4Trap("TowerBR",pt);
    	towerLogicalBR[i] = new G4LogicalVolume(tower,cu,name);
        //towerLogicalBR[i] = new G4LogicalVolume(tower,mix,name);
    	towerLogicalBR[i]->SetVisAttributes(towerVisAttr);
    	dimB->Getpt_PMTG(pt);
        //pmtg = new G4Trap("PMTGBR",pt);
        //PMTGLogicalBR[i] = new G4LogicalVolume(pmtg,Glass,name);
        //PMTGLogicalBR[i]->SetVisAttributes(PMTVisAttr);
        //G4cout << dimB->GetRM(0) << dimB->GetOrigin(0) << G4endl;
    	G4RotationMatrix* rm = new G4RotationMatrix();
    	rm->rotateX(-thetaofcenter);
    	G4ThreeVector c = dimB->GetOrigin(0);
        //G4cout << tower->GetZHalfLength() << " " << tower->GetYHalfLength1() << " " << tower->GetXHalfLength1() << " " << tower->GetXHalfLength2() << " "
        //       << tower->GetTanAlpha1() << " " << tower->GetYHalfLength2() << " " << tower->GetXHalfLength3() << " " << tower->GetXHalfLength4() << " "
        //       << tower->GetTanAlpha2() << G4endl;
        //G4cout << c << G4endl;
    	G4ThreeVector c_new(c.getY(),-c.getZ(),c.getX()-(innerR+0.5*length));
        //placing towers in barrel R
    	new G4PVPlacement(rm,c_new,towerLogicalBR[i],name,phiBLog,false,volnum,fCheckOverlaps);
        //made by me for new barrel structure
        //if (i>4) new G4PVPlacement(rm,c_new,towerLogicalBR[i],name,phiB2Log,false,volnum,fCheckOverlaps);
    	sprintf(name,"PMT%d",volnum);
        //new G4PVPlacement(dimB->GetRM(0),dimB->GetOrigin_PMTG(0),PMTGLogicalBR[i],name,worldLV,false,i,fCheckOverlaps);
        //for(int j=0;j<NbOfZRot;j++){
        //for(int j=0;j<1;j++){
        //new G4PVPlacement(dimB->GetRM(j),dimB->GetOrigin(j),towerLogicalBR[i],name,worldLogical,false,j,checkOverlaps);
        //new G4PVPlacement(dimB->GetRM(j),dimB->GetOrigin_PMTG(j),PMTGLogicalBR[i],name,worldLogical,false,j,checkOverlaps);
        //}
    	dimB->Getpt(pt);
    	sprintf(name,"fiber%d",volnum);
        //VERY IMPORTANT TO PLACE FIBERS
        //fiberBR(i,deltatheta_barrel[i]); //where fiberBR is defined?
    	G4cout<<"onefiberfromcenter2"<<G4endl;
        //    dimB->Getpt_PMTCath(pt);
        //pmtcath = new G4Trap("PMTCathBR",pt);
        //    PMTCathLogicalBR[i] = new G4LogicalVolume(pmtcath,Al,name);
        //new G4PVPlacement(0,G4ThreeVector(0,0,1.5*PMTT/2.-PMTT/4.),PMTCathLogicalBR[i],name,PMTGLogicalBR[i],false,0,fCheckOverlaps);
        //new G4LogicalSkinSurface("Photocath_surf",PMTCathLogicalBR[i],photocath_opsurf);
        //}//
    	fulltheta = fulltheta+deltatheta_barrel[i];
    	volnum++;
    }
    
    // barrel L
    G4cout << "Barrel L..." << G4endl;
    dimB->Rbool(0);
    thetaofcenter=0;
    fulltheta=0;
    
    for(int i=0;i<NbOfBarrel;i++){  //i<NbOfBarrel
        //if(volnum==99999999){//
    	thetaofcenter=fulltheta+deltatheta_barrel[i]/2.;
    	dimB->SetDeltaTheta(deltatheta_barrel[i]);
    	dimB->SetThetaOfCenter(thetaofcenter);
    	dimB->CalBasic();
    	dimB->Getpt(pt);

    	sprintf(name,"tower%d",volnum);
    	tower = new G4Trap("TowerBL",pt);
    	towerLogicalBL[i] = new G4LogicalVolume(tower,cu,name);
        //towerLogicalBL[i] = new G4LogicalVolume(tower,mix,name);
    	towerLogicalBL[i]->SetVisAttributes(towerVisAttr);
        //dimB->Getpt_PMTG(pt);
        //pmtg = new G4Trap("PMTGBL",pt);
        //PMTGLogicalBL[i] = new G4LogicalVolume(pmtg,Glass,name);
        //PMTGLogicalBL[i]->SetVisAttributes(visAttr);
    	G4RotationMatrix* rm = new G4RotationMatrix();
    	rm->rotateX(thetaofcenter);
        //G4RotationMatrix rm;
        //rm.rotateX(-thetaofcenter);
    	G4ThreeVector c = dimB->GetOrigin(0);
    	G4ThreeVector c_new(c.getY(),-c.getZ(),c.getX()-(innerR+0.5*length));
    	new G4PVPlacement(rm,c_new,towerLogicalBL[i],name,phiBLog,false,volnum,fCheckOverlaps);
        //new G4PVPlacement(dimB->GetRM(0),dimB->GetOrigin_PMTG(0),PMTGLogicalBL[i],name,phiDivLog,false,i,checkOverlaps);
    	//for(int j=0;j<NbOfZRot;j++){
            //for(int j=0;j<1;j++){
            //  new G4PVPlacement(dimB->GetRM(j),dimB->GetOrigin(j),towerLogicalBL[i],name,worldLogical,false,j,checkOverlaps);
            //new G4PVPlacement(dimB->GetRM(j),dimB->GetOrigin_PMTG(j),PMTGLogicalBL[i],name,worldLogical,false,j,checkOverlaps);
    	//}
    	dimB->Getpt(pt);
    	sprintf(name,"fiber%d",volnum);
        //VERY IMPORTANT TO PLACE FIBERS
    	//fiberBL(i,deltatheta_barrel[i]);

        //dimB->Getpt_PMTCath(pt);
        //pmtcath = new G4Trap("PMTCathBL",pt);
        //PMTCathLogicalBL[i] = new G4LogicalVolume(pmtcath,Al,name);
        //new G4PVPlacement(0,G4ThreeVector(0,0,1.5*PMTT/2.-PMTT/4.),PMTCathLogicalBL[i],name,PMTGLogicalBL[i],false,0,checkOverlaps);
        //new G4LogicalSkinSurface("Photocath_surf",PMTCathLogicalBL[i],photocath_opsurf);
        //}//
    	fulltheta = fulltheta+deltatheta_barrel[i];
    	volnum++;
    }
    
    //    cout<<"B trns last"<<dimB->GetTrns_Length()<<endl;
    cout<<"InnerR_New last"<<dimB->GetInnerR_new()<<endl;
    cout<<"fulltheta last"<<fulltheta<<endl;
    //    cout<<"deltatheta last"<<deltatheta_endcap<<endl;
    
    ///endcap///////////////////////////////////////////////////////////
    lastdeltatheta = deltatheta_endcap;
    dimE = new dimensionE();
    dimE->SetInnerR_new(innerR_Endcap);
    //    dimE->SetInnerR_new(3125.8);
    dimE->SetTower_height(tower_height);
    dimE->SetNumZRot(NbOfZRot);
    dimE->SetDeltaTheta(lastdeltatheta);
    dimE->SetPMTT(PMTT);
    
    ///////////////////////////////
    // endcap R
    
    //for(int i=0;i<NbOfBarrel;i++) thetaB += deltatheta_barrel[i];
    
    G4cout << "Endcap R..." << G4endl;
    fulltheta = thetaB;
    //    fulltheta = 0.95717;
    //    fulltheta = 0.957164;
    dimE->Rbool(1);
    
    for(int i=0;i<NbOfEndcap;i++){ //i<NbOfEndcap
        //if(volnum==99999999){//
    	thetaofcenter=fulltheta+lastdeltatheta/2.;
    	dimE->SetThetaOfCenter(thetaofcenter);
    	dimE->CalBasic();
    	dimE->Getpt(pt);
        //G4cout << i << " : " << pt[0] << pt[1] << pt[2] << pt[3] << G4endl;
    	sprintf(name,"tower%d",volnum);
    	G4Trap* tower = new G4Trap("TowerER",pt);

    	towerLogicalER[i] = new G4LogicalVolume(tower,cu,name);
        //towerLogicalER[i] = new G4LogicalVolume(tower,mix,name);
    	towerLogicalER[i]->SetVisAttributes(towerVisAttr);
        //dimE->Getpt_PMTG(pt);
        //pmtg = new G4Trap("PMTGER",pt);
        //PMTGLogicalER[i] = new G4LogicalVolume(pmtg,Glass,name);
        //PMTGLogicalER[i]->SetVisAttributes(visAttr);
    	G4RotationMatrix* rm = new G4RotationMatrix();
    	rm->rotateX(-thetaofcenter+thetaB+0.5*thetaE);
    	G4ThreeVector c = dimE->GetOrigin(0);
    	double r = c.mag();
    	double x = r*sin(thetaB+0.5*thetaE-thetaofcenter);
    	double z = r*cos(thetaB+0.5*thetaE-thetaofcenter)-rEC;
    	G4ThreeVector c_new(c.getY(),x,z);
        new G4PVPlacement(rm,c_new,towerLogicalER[i],name,phiERLog,false,volnum,fCheckOverlaps);
    	//for(int j=0;j<NbOfZRot;j++){
            //for(int j=0;j<1;j++){
            //new G4PVPlacement(dimE->GetRM(j),dimE->GetOrigin(j),towerLogicalER[i],name,worldLogical,false,j,checkOverlaps);
            //new G4PVPlacement(dimE->GetRM(j),dimE->GetOrigin_PMTG(j),PMTGLogicalER[i],name,worldLogical,false,j,checkOverlaps);
    	//}
    	dimE->Getpt(pt);
        //VERY IMPORTANT FOR PLACING FIBERS
        //fiberER(i,lastdeltatheta);

        //dimE->Getpt_PMTCath(pt);
        //pmtcath = new G4Trap("PMTCathER",pt);
        //PMTCathLogicalER[i] = new G4LogicalVolume(pmtcath,Al,name);
        //new G4PVPlacement(0,G4ThreeVector(0,0,1.5*PMTT/2.-PMTT/4.),PMTCathLogicalER[i],name,PMTGLogicalER[i],false,0,checkOverlaps);
        //new G4LogicalSkinSurface("Photocath_surf",PMTCathLogicalER[i],photocath_opsurf);
        //}//
    	fulltheta = fulltheta+lastdeltatheta;
    	volnum++;
    }
    
    // endcap L
    G4cout << "Endcap L..." << G4endl;
    fulltheta = thetaB;
    //    fulltheta = 0.95717;
    //    fulltheta = 0.957164;
    dimE->Rbool(0);
    
    for(int i=0;i<10;i++){
        //if(volnum==99999999){//
    	thetaofcenter=fulltheta+lastdeltatheta/2.;
    	dimE->SetThetaOfCenter(thetaofcenter);
    	dimE->CalBasic();
    	dimE->Getpt(pt);

    	sprintf(name,"tower%d",volnum);
    	G4Trap* tower = new G4Trap("TowerEL",pt);
    	towerLogicalEL[i] = new G4LogicalVolume(tower,cu,name);
        //towerLogicalEL[i] = new G4LogicalVolume(tower,mix,name);
    	towerLogicalEL[i]->SetVisAttributes(towerVisAttr);
        //dimE->Getpt_PMTG(pt);
        //pmtg = new G4Trap("PMTGEL",pt);
        //PMTGLogicalEL[i] = new G4LogicalVolume(pmtg,Glass,name);
        //PMTGLogicalEL[i]->SetVisAttributes(visAttr);
    	G4RotationMatrix* rm = new G4RotationMatrix();
    	rm->rotateX(thetaofcenter-thetaB-0.5*thetaE);
    	G4ThreeVector c = dimE->GetOrigin(0);
    	double r = c.mag();
    	double x = r*sin(thetaofcenter-thetaB-0.5*thetaE);
    	double z = r*cos(thetaofcenter-thetaB-0.5*thetaE)-rEC;
    	G4ThreeVector c_new(c.getY(),x,z);
        if (i==9) new G4PVPlacement(rm,c_new,towerLogicalEL[i],name,phiELLog,false,volnum,fCheckOverlaps);
    	//for(int j=0;j<NbOfZRot;j++){
    	//	for(int j=0;j<1;j++){}
            //new G4PVPlacement(dimE->GetRM(j),dimE->GetOrigin(j),towerLogicalEL[i],name,worldLogical,false,j,fCheckOverlaps);}
            //new G4PVPlacement(dimE->GetRM(j),dimE->GetOrigin_PMTG(j),PMTGLogicalEL[i],name,worldLogical,false,j,checkOverlaps);}
    	//}

    dimE->Getpt(pt);
    fiberEL(i,lastdeltatheta);

        //dimE->Getpt_PMTCath(pt);
        //pmtcath = new G4Trap("PMTCathEL",pt);
        //PMTCathLogicalEL[i] = new G4LogicalVolume(pmtcath,Al,name);
        //new G4PVPlacement(0,G4ThreeVector(0,0,1.5*PMTT/2.-PMTT/4.),PMTCathLogicalEL[i],name,PMTGLogicalEL[i],false,0,checkOverlaps);
        //new G4LogicalSkinSurface("Photocath_surf",PMTCathLogicalEL[i],photocath_opsurf);
        //}//
    fulltheta = fulltheta+lastdeltatheta;
    volnum++;
}

return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DreamDetectorConstruction::ConstructSDandField()
{
    // Create global magnetic field messenger,
    // Uniform magnetic field is then created automatically if
    // the field value is not zero
	G4ThreeVector fieldValue = G4ThreeVector();
	fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
	fMagFieldMessenger->SetVerboseLevel(1);

    // Register the field messenger for deleting
	G4AutoDelete::Register(fMagFieldMessenger);
}

void DreamDetectorConstruction::fiberBR(G4int i,G4double deltatheta_){

	vector<G4double> temp_x;
    vector<G4double> temp_y;// vector for grid calculation
    
    temp_x.clear();
    temp_y.clear();
    temp_x.push_back(0.*mm);
    temp_y.push_back(0.*mm);
    
    for(int j = 0 ; j<1000;j++){
    	temp_x.push_back(0.*mm);
    	temp_y.push_back(1.5*(j+1)*mm);
    	temp_x.push_back(0.*mm);
    	temp_y.push_back(-1.5*(j+1)*mm);
    }
    for(int i=0; i<1000;i++){
    	temp_x.push_back(1.5*(i+1)*mm);
    	temp_y.push_back(0.*mm);
    	temp_x.push_back(-1.5*(i+1)*mm);
    	temp_y.push_back(0.*mm);
    }
    for(int j = 0;j<1000;j++){
    	for(int i = 0;i<1000;i++){
    		temp_x.push_back(1.5*(i+1)*mm);
    		temp_y.push_back(1.5*(j+1)*mm);
    		temp_x.push_back(-1.5*(i+1)*mm);
    		temp_y.push_back(-1.5*(j+1)*mm);
    		temp_x.push_back(-1.5*(i+1)*mm);
    		temp_y.push_back(1.5*(j+1)*mm);
    		temp_x.push_back(1.5*(i+1)*mm);
    		temp_y.push_back(-1.5*(j+1)*mm);
    	}
    }
    
    G4double center_x;
    G4double center_y;
    
    int hi1;
    int hi2;
    int hi1reminder; // 0 -> even 1-> odd
    int hi2reminder;
    
    int numx;
    int numy;
    int reminderx;
    int remindery;
    
    G4ThreeVector v1 = dimB->GetV1();
    G4ThreeVector v2 = dimB->GetV2();
    G4ThreeVector v3 = dimB->GetV3();
    G4ThreeVector v4 = dimB->GetV4();
    
    G4double innerSide_half = dimB->GetInnerR_new()*tan(deltatheta_/2.);
    G4double outerSide_half = (dimB->GetInnerR_new()+tower_height)*tan(deltatheta_/2.);
    
    G4double theta_bc=0;
    
    
    //1. 4 types of grid coordinate(x,y) =  (o,o),(e,o),(o,e),(e,e) { e= 0, o=1 }
    int type_x_BR;
    int type_y_BR;
    
    numx = (int)(((v4.getX()*tan(phi_unit/2.)*2)-1.*mm)/1.5*mm);
    numy = (int)((outerSide_half*2-1.*mm)/(1.5*mm));
    reminderx = numx%2;
    remindery = numy%2;
    if(reminderx == 1) type_x_BR=0;
    if(reminderx == 0) type_x_BR=1;
    if(remindery == 1) type_y_BR=0;
    if(remindery == 0) type_y_BR=1;
    
    ////2. aplying the boundary conditions, reject the fibre which near the boundary of tower surface, get coord of cetres of fibre
    //prepare the grid vector : v2 ~ 114 ( #76 ) , outerSide ~ 132 (#88)
    
    ////////////////////////////////////////////
    vector<G4double> center_x_BR;
    vector<G4double> center_y_BR;
    vector<G4int> bool_cfiber_BR;// 0 -> c fiber 1 s fiber
    
    // select the type of grid coord
    theta_bc=atan(2*outerSide_half/((v2.getX()-v4.getX())*tan(phi_unit/2.)));
    
    if(type_x_BR==1&&type_y_BR==1){//(o,o)
    	for(int num = 0;num<temp_x.size();num++){
    		center_x = temp_x.at(num);
    		center_y = temp_y.at(num);
    		hi1 = fabs(temp_x.at(num))/(1.5*mm);
    		hi2 = fabs(temp_y.at(num))/(1.5*mm);
            hi1reminder = hi1%2; // 0 -> even 1-> odd
            hi2reminder = hi2%2;
            
            
            if(((center_y-cos(theta_bc)*clad_C_rMax)>(2*2*outerSide_half/(2*v2.getX()*tan(phi_unit/2.)-2*v4.getX()*tan(phi_unit/2.))*(center_x+sin(theta_bc)*clad_C_rMax-v2.getX()*tan(phi_unit/2.))+outerSide_half))
            	&&
            	((center_y+clad_C_rMax)<(outerSide_half))
            	&&
            	((center_y-clad_C_rMax)>(-outerSide_half))
            	&&
            	((center_y-cos(theta_bc)*clad_C_rMax)>(2*2*outerSide_half/(2*v4.getX()*tan(phi_unit/2.)-2*v2.getX()*tan(phi_unit/2.))*(center_x-sin(theta_bc)*clad_C_rMax+v2.getX()*tan(phi_unit/2.))+outerSide_half
            		)))
            {
            	center_x_BR.push_back(center_x);
            	center_y_BR.push_back(center_y);
            	if(((hi1reminder==0)&&(hi2reminder==0))||((hi1reminder==1)&&(hi2reminder==1)))
            		bool_cfiber_BR.push_back(0);
            	else bool_cfiber_BR.push_back(1);
            }
            
        }
    }
    
    if(type_x_BR==0&&type_y_BR==0){//(e,e)

    	for(int num = 0;num<temp_x.size();num++){
    		center_x = temp_x.at(num)+0.75*mm;
    		center_y = temp_y.at(num)-0.75*mm;

    		hi1 = fabs(temp_x.at(num))/(1.5*mm);
    		hi2 = fabs(temp_y.at(num))/(1.5*mm);
            hi1reminder = hi1%2; // 0 -> even 1-> odd
            hi2reminder = hi2%2;
            
            if(((center_y-cos(theta_bc)*clad_C_rMax)>(2*2*outerSide_half/(2*v2.getX()*tan(phi_unit/2.)-2*v4.getX()*tan(phi_unit/2.))*(center_x+sin(theta_bc)*clad_C_rMax-v2.getX()*tan(phi_unit/2.))+outerSide_half))
            	&&
            	((center_y+clad_C_rMax)<(outerSide_half))
            	&&
            	((center_y-clad_C_rMax)>(-outerSide_half))
            	&&
            	((center_y-cos(theta_bc)*clad_C_rMax)>(2*2*outerSide_half/(2*v4.getX()*tan(phi_unit/2.)-2*v2.getX()*tan(phi_unit/2.))*(center_x-sin(theta_bc)*clad_C_rMax+v2.getX()*tan(phi_unit/2.))+outerSide_half
            		)))
            {
            	center_x_BR.push_back(center_x);
            	center_y_BR.push_back(center_y);
            	if(((hi1reminder==0)&&(hi2reminder==0))||((hi1reminder==1)&&(hi2reminder==1)))
            		bool_cfiber_BR.push_back(0);
            	else bool_cfiber_BR.push_back(1);
            }
        }
    }
    
    
    if(type_x_BR==0&&type_y_BR==1){//(e,o)

    	for(int num = 0;num<temp_x.size();num++){
    		center_x = temp_x.at(num)+0.75*mm;
    		center_y = temp_y.at(num);

    		hi1 = fabs(temp_x.at(num))/(1.5*mm);
    		hi2 = fabs(temp_y.at(num))/(1.5*mm);
            hi1reminder = hi1%2; // 0 -> even 1-> odd
            hi2reminder = hi2%2;
            
            if(((center_y-cos(theta_bc)*clad_C_rMax)>(2*2*outerSide_half/(2*v2.getX()*tan(phi_unit/2.)-2*v4.getX()*tan(phi_unit/2.))*(center_x+sin(theta_bc)*clad_C_rMax-v2.getX()*tan(phi_unit/2.))+outerSide_half))
            	&&
            	((center_y+clad_C_rMax)<(outerSide_half))
            	&&
            	((center_y-clad_C_rMax)>(-outerSide_half))
            	&&
            	((center_y-cos(theta_bc)*clad_C_rMax)>(2*2*outerSide_half/(2*v4.getX()*tan(phi_unit/2.)-2*v2.getX()*tan(phi_unit/2.))*(center_x-sin(theta_bc)*clad_C_rMax+v2.getX()*tan(phi_unit/2.))+outerSide_half
            		)))
            {
            	center_x_BR.push_back(center_x);
            	center_y_BR.push_back(center_y);
            	if(((hi1reminder==0)&&(hi2reminder==0))||((hi1reminder==1)&&(hi2reminder==1)))
            		bool_cfiber_BR.push_back(0);
            	else bool_cfiber_BR.push_back(1);
            }
            
        }
    }
    
    if(type_x_BR==1&&type_y_BR==0){//(o,e)

    	for(int num = 0;num<temp_x.size();num++){
    		center_x = temp_x.at(num);
    		center_y = temp_y.at(num)-0.75*mm;

    		hi1 = fabs(temp_x.at(num))/(1.5*mm);
    		hi2 = fabs(temp_y.at(num))/(1.5*mm);
            hi1reminder = hi1%2; // 0 -> even 1-> odd
            hi2reminder = hi2%2;
            
            if(((center_y-cos(theta_bc)*clad_C_rMax)>(2*2*outerSide_half/(2*v2.getX()*tan(phi_unit/2.)-2*v4.getX()*tan(phi_unit/2.))*(center_x+sin(theta_bc)*clad_C_rMax-v2.getX()*tan(phi_unit/2.))+outerSide_half))
            	&&
            	((center_y+clad_C_rMax)<(outerSide_half))
            	&&
            	((center_y-clad_C_rMax)>(-outerSide_half))
            	&&
            	((center_y-cos(theta_bc)*clad_C_rMax)>(2*2*outerSide_half/(2*v4.getX()*tan(phi_unit/2.)-2*v2.getX()*tan(phi_unit/2.))*(center_x-sin(theta_bc)*clad_C_rMax+v2.getX()*tan(phi_unit/2.))+outerSide_half
            		)))
            {
            	center_x_BR.push_back(center_x);
            	center_y_BR.push_back(center_y);
            	if(((hi1reminder==0)&&(hi2reminder==0))||((hi1reminder==1)&&(hi2reminder==1)))
            		bool_cfiber_BR.push_back(0);
            	else bool_cfiber_BR.push_back(1);
            }
            
        }
    }
    for(int j = 0; j<center_x_BR.size();j++){
        // determine z value for center of fibre
        //sprintf(name,"fiber%d",j);
        //if(fabs(center_x_BR.at(j))>1) continue;
    	double z = tower->GetZHalfLength();
        //G4cout << i << "," << j << ": " << G4ThreeVector(center_x_BR.at(j),center_y_BR.at(j),0) << " " << z << G4endl;
    	bool outside = false;
    	for(int ip = 0;ip<4;ip++){
    		TrapSidePlane plane = tower->GetSidePlane(ip);
    		double zpoint = (-plane.a*center_x_BR.at(j)-plane.b*center_y_BR.at(j)-plane.d)/plane.c;
    		outside = (tower->Inside(G4ThreeVector(center_x_BR.at(j),center_y_BR.at(j),zpoint))==kOutside);
    		G4cout << ip << ": " << plane.a << " " << plane.b << " " << plane.c << " " << plane.d << " " << outside << G4endl;
    		if(!outside){
    			G4ThreeVector normal = tower->SurfaceNormal(G4ThreeVector(center_x_BR.at(j),center_y_BR.at(j),zpoint));
    			double angle = normal.angle(G4ThreeVector(0,0,-1));
    			double shift = fabs(clad_C_rMax/tan(0.5*M_PI-angle));
    			int length = z - zpoint - shift;
    			G4cout << length << G4endl;
    			if(length>=1&&length<=2*z){
    				new G4PVPlacement(0,G4ThreeVector(center_x_BR.at(j),center_y_BR.at(j),z-0.5*length),
    					(bool_cfiber_BR.at(j)==0)?fiberCLog[length]:fiberSLog[length],name,towerLogicalBR[i],false,j,false);

    			}
    			G4cout << G4ThreeVector(center_x_BR.at(j),center_y_BR.at(j),zpoint) << " " << length << " " << shift << " " << angle << G4endl;
    			break;
    		}
    	}
    	if(outside) {new G4PVPlacement(0,G4ThreeVector(center_x_BR.at(j),center_y_BR.at(j),0),
    		(bool_cfiber_BR.at(j)==0)?fiberCLog[2*z]:fiberSLog[2*z],name,towerLogicalBR[i],false,j,false);
				//G4cout<<"djskdfjkslfjdlkfjdlkfjdlk"<<G4endl;
    }

        if(bool_cfiber_BR.at(j)==0){//c fibre
            /*
             intersect = new G4IntersectionSolid("fiber_",fiber,tower,0,G4ThreeVector(-center_x_BR.at(j),-center_y_BR.at(j),0.));
             fiberLogical_BR[i].push_back(new G4LogicalVolume(intersect,clad_C_Material,name));
             new G4PVPlacement(0,G4ThreeVector(center_x_BR.at(j),center_y_BR.at(j),0),fiberLogical_BR[i].at(j),name,towerLogicalBR[i],false,j,fCheckOverlaps);
             
             intersect_ = new G4IntersectionSolid("fiber_",fiberC,tower,0,G4ThreeVector(-center_x_BR.at(j),-center_y_BR.at(j),0.));
             fiberLogical_BR_[i].push_back(new G4LogicalVolume(intersect_,core_C_Material,name));
             //new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),fiberLogical_BR_[i].at(j),name,fiberLogical_BR[i].at(j),false,j,checkOverlaps);
             new G4PVPlacement(0,G4ThreeVector(center_x_BR.at(j),center_y_BR.at(j),0),fiberCladCLog,name,towerLogicalBR[i],false,j,fCheckOverlaps);
             
             visAttr = new G4VisAttributes(G4Colour(0.5,0.5,0.5,0.5));
             //visAttr->SetForceSolid(true);
             visAttr->SetVisibility(true);
             fiberLogical_BR[i].at(j)->SetVisAttributes(visAttr);
             fVisAttributes.push_back(visAttr);
             
             visAttr = new G4VisAttributes(G4Colour(0.,0.,1.0,0.5));
             //visAttr->SetForceSolid(true);
             visAttr->SetVisibility(true);
             fiberLogical_BR_[i].at(j)->SetVisAttributes(visAttr);
             fVisAttributes.push_back(visAttr);
             */
        }
        else // s fibre
        {
            /*
             intersect = new G4IntersectionSolid("fiber_",fiber,tower,0,G4ThreeVector(-center_x_BR.at(j),-center_y_BR.at(j),0.));
             fiberLogical_BR[i].push_back(new G4LogicalVolume(intersect,clad_S_Material,name));
             new G4PVPlacement(0,G4ThreeVector(center_x_BR.at(j),center_y_BR.at(j),0),fiberLogical_BR[i].at(j),name,towerLogicalBR[i],false,j,checkOverlaps);
             
             intersect_ = new G4IntersectionSolid("fiber_",fiberS,tower,0,G4ThreeVector(-center_x_BR.at(j),-center_y_BR.at(j),0.));
             fiberLogical_BR_[i].push_back(new G4LogicalVolume(intersect_,core_S_Material,name));
             new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),fiberLogical_BR_[i].at(j),name,fiberLogical_BR[i].at(j),false,j,checkOverlaps);
             //new G4PVPlacement(0,G4ThreeVector(center_x_BR.at(j),center_y_BR.at(j),0),fiberCladSLog,name,towerLogicalBR[i],false,j,checkOverlaps);
             
             visAttr = new G4VisAttributes(G4Colour(0.5,0.5,0.5,0.5));
             //visAttr->SetForceSolid(true);
             visAttr->SetVisibility(true);
             fiberLogical_BR[i].at(j)->SetVisAttributes(visAttr);
             fVisAttributes.push_back(visAttr);
             
             visAttr = new G4VisAttributes(G4Colour(1.0,1.0,0.,0.5));
             //visAttr->SetForceSolid(true);
             visAttr->SetVisibility(true);
             fiberLogical_BR_[i].at(j)->SetVisAttributes(visAttr);
             fVisAttributes.push_back(visAttr);
             */
        }
        
    }
    
}

void DreamDetectorConstruction::fiberBL(G4int i, G4double deltatheta_){
	vector<G4double> temp_x;
    vector<G4double> temp_y;// vector for grid calculation
    
    temp_x.clear();
    temp_y.clear();
    temp_x.push_back(0.*mm);
    temp_y.push_back(0.*mm);
    
    for(int j = 0 ; j<1000;j++){
    	temp_x.push_back(0.*mm);
    	temp_y.push_back(1.5*(j+1)*mm);
    	temp_x.push_back(0.*mm);
    	temp_y.push_back(-1.5*(j+1)*mm);
    }
    for(int i=0; i<1000;i++){
    	temp_x.push_back(1.5*(i+1)*mm);
    	temp_y.push_back(0.*mm);
    	temp_x.push_back(-1.5*(i+1)*mm);
    	temp_y.push_back(0.*mm);
    }
    for(int j = 0;j<1000;j++){
    	for(int i = 0;i<1000;i++){
    		temp_x.push_back(1.5*(i+1)*mm);
    		temp_y.push_back(1.5*(j+1)*mm);
    		temp_x.push_back(-1.5*(i+1)*mm);
    		temp_y.push_back(-1.5*(j+1)*mm);
    		temp_x.push_back(-1.5*(i+1)*mm);
    		temp_y.push_back(1.5*(j+1)*mm);
    		temp_x.push_back(1.5*(i+1)*mm);
    		temp_y.push_back(-1.5*(j+1)*mm);
    	}
    }
    
    G4double center_x;
    G4double center_y;
    int hi1;
    int hi2;
    int hi1reminder; // 0 -> even 1-> odd
    int hi2reminder;
    
    int numx;
    int numy;
    int reminderx;
    int remindery;
    
    G4ThreeVector v1 = dimB->GetV1();
    G4ThreeVector v2 = dimB->GetV2();
    G4ThreeVector v3 = dimB->GetV3();
    G4ThreeVector v4 = dimB->GetV4();
    
    G4double innerSide_half = dimB->GetInnerR_new()*tan(deltatheta_/2.);
    G4double outerSide_half = (dimB->GetInnerR_new()+tower_height)*tan(deltatheta_/2.);
    
    G4double theta_bc=0;
    /////////////////////
    //1. 4 types of grid coordinate(x,y) =  (o,o),(e,o),(o,e),(e,e) { e= 0, o=1 }
    
    
    
    int type_x_BL;
    int type_y_BL;
    
    numx = (int)(((v4.getX()*tan(phi_unit/2.)*2)-1.*mm)/1.5*mm);
    numy = (int)((outerSide_half*2-1.*mm)/(1.5*mm));
    reminderx = numx%2;
    remindery = numy%2;
    
    if(reminderx == 1) type_x_BL=0;
    if(reminderx == 0) type_x_BL=1;
    if(remindery == 1) type_y_BL=0;
    if(remindery == 0) type_y_BL=1;
    
    
    ////2. aplying the boundary conditions, reject the fibre which near the boundary of tower surface, get coord of cetres of fibre
    
    
    vector<G4double> center_x_BL;
    vector<G4double> center_y_BL;
    vector<G4int> bool_cfiber_BL;// 0 -> c fiber 1 s fiber
    
    
    // select the type of grid coord
    
    theta_bc=atan(2*outerSide_half/((v2.getX()-v4.getX())*tan(phi_unit/2.)));
    if(type_x_BL==1&&type_y_BL==1){//(o,o)
    	for(int num = 0;num<temp_x.size();num++){
    		center_x = temp_x.at(num);
    		center_y = temp_y.at(num);
    		hi1 = fabs(temp_x.at(num))/(1.5*mm);
    		hi2 = fabs(temp_y.at(num))/(1.5*mm);
            hi1reminder = hi1%2; // 0 -> even 1-> odd
            hi2reminder = hi2%2;
            
            if(((center_y+cos(theta_bc)*clad_C_rMax)<(2*2*outerSide_half/(2*v4.getX()*tan(phi_unit/2.)-2*v2.getX()*tan(phi_unit/2.))*(center_x+sin(theta_bc)*clad_C_rMax-v4.getX()*tan(phi_unit/2.))+outerSide_half))
            	&&
            	((center_y+clad_C_rMax)<(outerSide_half))
            	&&
            	((center_y-clad_C_rMax)>(-outerSide_half))
            	&&
            	((center_y+cos(theta_bc)*clad_C_rMax)<(2*2*outerSide_half/(2*v2.getX()*tan(phi_unit/2.)-2*v4.getX()*tan(phi_unit/2.))*(center_x-sin(theta_bc)*clad_C_rMax+v4.getX()*tan(phi_unit/2.))+outerSide_half
            		)))
            {
            	center_x_BL.push_back(center_x);
            	center_y_BL.push_back(center_y);
            	if(((hi1reminder==0)&&(hi2reminder==0))||((hi1reminder==1)&&(hi2reminder==1)))
            		bool_cfiber_BL.push_back(0);
            	else bool_cfiber_BL.push_back(1);
            }
            
        }
    }
    
    if(type_x_BL==0&&type_y_BL==0){//(e,e)

    	for(int num = 0;num<temp_x.size();num++){
    		center_x = temp_x.at(num)+0.75*mm;
    		center_y = temp_y.at(num)-0.75*mm;

    		hi1 = fabs(temp_x.at(num))/(1.5*mm);
    		hi2 = fabs(temp_y.at(num))/(1.5*mm);
            hi1reminder = hi1%2; // 0 -> even 1-> odd
            hi2reminder = hi2%2;
            if(((center_y+cos(theta_bc)*clad_C_rMax)<(2*2*outerSide_half/(2*v4.getX()*tan(phi_unit/2.)-2*v2.getX()*tan(phi_unit/2.))*(center_x+sin(theta_bc)*clad_C_rMax-v4.getX()*tan(phi_unit/2.))+outerSide_half))
            	&&
            	((center_y+clad_C_rMax)<(outerSide_half))
            	&&
            	((center_y-clad_C_rMax)>(-outerSide_half))
            	&&
            	((center_y+cos(theta_bc)*clad_C_rMax)<(2*2*outerSide_half/(2*v2.getX()*tan(phi_unit/2.)-2*v4.getX()*tan(phi_unit/2.))*(center_x-sin(theta_bc)*clad_C_rMax+v4.getX()*tan(phi_unit/2.))+outerSide_half
            		)))
            {
            	center_x_BL.push_back(center_x);
            	center_y_BL.push_back(center_y);
            	if(((hi1reminder==0)&&(hi2reminder==0))||((hi1reminder==1)&&(hi2reminder==1)))
            		bool_cfiber_BL.push_back(0);
            	else bool_cfiber_BL.push_back(1);
            }
        }
    }
    
    
    if(type_x_BL==0&&type_y_BL==1){//(e,o)

    	for(int num = 0;num<temp_x.size();num++){
    		center_x = temp_x.at(num)+0.75*mm;
    		center_y = temp_y.at(num);

    		hi1 = fabs(temp_x.at(num))/(1.5*mm);
    		hi2 = fabs(temp_y.at(num))/(1.5*mm);
            hi1reminder = hi1%2; // 0 -> even 1-> odd
            hi2reminder = hi2%2;
            
            if(((center_y+cos(theta_bc)*clad_C_rMax)<(2*2*outerSide_half/(2*v4.getX()*tan(phi_unit/2.)-2*v2.getX()*tan(phi_unit/2.))*(center_x+sin(theta_bc)*clad_C_rMax-v4.getX()*tan(phi_unit/2.))+outerSide_half))
            	&&
            	((center_y+clad_C_rMax)<(outerSide_half))
            	&&
            	((center_y-clad_C_rMax)>(-outerSide_half))
            	&&
            	((center_y+cos(theta_bc)*clad_C_rMax)<(2*2*outerSide_half/(2*v2.getX()*tan(phi_unit/2.)-2*v4.getX()*tan(phi_unit/2.))*(center_x-sin(theta_bc)*clad_C_rMax+v4.getX()*tan(phi_unit/2.))+outerSide_half
            		)))
            {
            	center_x_BL.push_back(center_x);
            	center_y_BL.push_back(center_y);
            	if(((hi1reminder==0)&&(hi2reminder==0))||((hi1reminder==1)&&(hi2reminder==1)))
            		bool_cfiber_BL.push_back(0);
            	else bool_cfiber_BL.push_back(1);
            }
            
        }
    }
    if(type_x_BL==1&&type_y_BL==0){//(o,e)

    	for(int num = 0;num<temp_x.size();num++){
    		center_x = temp_x.at(num);
    		center_y = temp_y.at(num)-0.75*mm;

    		hi1 = fabs(temp_x.at(num))/(1.5*mm);
    		hi2 = fabs(temp_y.at(num))/(1.5*mm);
            hi1reminder = hi1%2; // 0 -> even 1-> odd
            hi2reminder = hi2%2;
            
            if(((center_y+cos(theta_bc)*clad_C_rMax)<(2*2*outerSide_half/(2*v4.getX()*tan(phi_unit/2.)-2*v2.getX()*tan(phi_unit/2.))*(center_x+sin(theta_bc)*clad_C_rMax-v4.getX()*tan(phi_unit/2.))+outerSide_half))
            	&&
            	((center_y+clad_C_rMax)<(outerSide_half))
            	&&
            	((center_y-clad_C_rMax)>(-outerSide_half))
            	&&
            	((center_y+cos(theta_bc)*clad_C_rMax)<(2*2*outerSide_half/(2*v2.getX()*tan(phi_unit/2.)-2*v4.getX()*tan(phi_unit/2.))*(center_x-sin(theta_bc)*clad_C_rMax+v4.getX()*tan(phi_unit/2.))+outerSide_half
            		)))
            {
            	center_x_BL.push_back(center_x);
            	center_y_BL.push_back(center_y);
            	if(((hi1reminder==0)&&(hi2reminder==0))||((hi1reminder==1)&&(hi2reminder==1)))
            		bool_cfiber_BL.push_back(0);
            	else bool_cfiber_BL.push_back(1);
            }
            
        }
    }
    
    for(int j = 0; j<center_x_BL.size();j++){
        // determine z value for center of fibre
    	double z = tower->GetZHalfLength();
    	bool outside = false;
    	for(int ip = 0;ip<4;ip++){
    		TrapSidePlane plane = tower->GetSidePlane(ip);
    		double zpoint = (-plane.a*center_x_BL.at(j)-plane.b*center_y_BL.at(j)-plane.d)/plane.c;
    		outside = (tower->Inside(G4ThreeVector(center_x_BL.at(j),center_y_BL.at(j),zpoint))==kOutside);
    		if(!outside){
    			G4ThreeVector normal = tower->SurfaceNormal(G4ThreeVector(center_x_BL.at(j),center_y_BL.at(j),zpoint));
    			double angle = normal.angle(G4ThreeVector(0,0,-1));
    			double shift = fabs(clad_C_rMax/tan(0.5*M_PI-angle));
    			int length = z - zpoint - shift;
    			if(length>=1&&length<=2*z){
    				new G4PVPlacement(0,G4ThreeVector(center_x_BL.at(j),center_y_BL.at(j),z-0.5*length),
    					(bool_cfiber_BL.at(j)==0)?fiberCLog[length]:fiberSLog[length],name,towerLogicalBL[i],false,j,false);
    			}
    			break;
    		}
    	}
    	if(outside) new G4PVPlacement(0,G4ThreeVector(center_x_BL.at(j),center_y_BL.at(j),0),
    		(bool_cfiber_BL.at(j)==0)?fiberCLog[2*z]:fiberSLog[2*z],name,towerLogicalBL[i],false,j,false);
        if(bool_cfiber_BL.at(j)==0){//c fibre

            //intersect = new G4IntersectionSolid("fiber_",fiber,tower,0,G4ThreeVector(-center_x_BL.at(j),-center_y_BL.at(j),0.));
            //  fiberLogical_BL[i].push_back(new G4LogicalVolume(intersect,clad_C_Material,name));
            //  new G4PVPlacement(0,G4ThreeVector(center_x_BL.at(j),center_y_BL.at(j),0.),fiberLogical_BL[i].at(j),name,towerLogicalBL[i],false,j,checkOverlaps);

            //  intersect_ = new G4IntersectionSolid("fiber_",fiberC,tower,0,G4ThreeVector(-center_x_BL.at(j),-center_y_BL.at(j),0.));
            //  fiberLogical_BL_[i].push_back(new G4LogicalVolume(intersect_,core_C_Material,name));
            //  new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),fiberLogical_BL_[i].at(j),name,fiberLogical_BL[i].at(j),false,j,checkOverlaps);

            //  visAttr = new G4VisAttributes(G4Colour(0.5,0.5,0.5,0.5));
            //visAttr->SetForceSolid(true);
            //  visAttr->SetVisibility(true);
            //  fiberLogical_BL[i].at(j)->SetVisAttributes(visAttr);
            //  fVisAttributes.push_back(visAttr);

            //  visAttr = new G4VisAttributes(G4Colour(0,0,1,0.5));
            //visAttr->SetForceSolid(true);
            //  visAttr->SetVisibility(true);
            //  fiberLogical_BL_[i].at(j)->SetVisAttributes(visAttr);
            //  fVisAttributes.push_back(visAttr);
            //new G4PVPlacement(0,G4ThreeVector(center_x_BL.at(j),center_y_BL.at(j),0),fiberCladCLog,name,towerLogicalBL[i],false,j,checkOverlaps);
        }
        
        else // s fibre
        {

            /*
             intersect = new G4IntersectionSolid("fiber_",fiber,tower,0,G4ThreeVector(-center_x_BL.at(j),-center_y_BL.at(j),0.));
             fiberLogical_BL[i].push_back(new G4LogicalVolume(intersect,clad_S_Material,name));
             new G4PVPlacement(0,G4ThreeVector(center_x_BL.at(j),center_y_BL.at(j),0.),fiberLogical_BL[i].at(j),name,towerLogicalBL[i],false,j,checkOverlaps);
             
             intersect_ = new G4IntersectionSolid("fiber_",fiberS,tower,0,G4ThreeVector(-center_x_BL.at(j),-center_y_BL.at(j),0.));
             fiberLogical_BL_[i].push_back(new G4LogicalVolume(intersect_,core_S_Material,name));
             new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),fiberLogical_BL_[i].at(j),name,fiberLogical_BL[i].at(j),false,j,checkOverlaps);
             
             visAttr = new G4VisAttributes(G4Colour(0.5,0.5,0.5,0.5));
             //visAttr->SetForceSolid(true);
             visAttr->SetVisibility(true);
             fiberLogical_BL[i].at(j)->SetVisAttributes(visAttr);
             fVisAttributes.push_back(visAttr);
             
             visAttr = new G4VisAttributes(G4Colour(1,1,0,0.5));
             //visAttr->SetForceSolid(true);
             visAttr->SetVisibility(true);
             fiberLogical_BL_[i].at(j)->SetVisAttributes(visAttr);
             fVisAttributes.push_back(visAttr);
             */
            //new G4PVPlacement(0,G4ThreeVector(center_x_BL.at(j),center_y_BL.at(j),0),fiberCladSLog,name,towerLogicalBL[i],false,j,checkOverlaps);
        }
    }
}

void DreamDetectorConstruction::fiberER(G4int i,G4double deltatheta_){
	vector<G4double> temp_x;
    vector<G4double> temp_y;// vector for grid calculation
    
    temp_x.clear();
    temp_y.clear();
    temp_x.push_back(0.*mm);
    temp_y.push_back(0.*mm);
    
    for(int j = 0 ; j<1000;j++){
    	temp_x.push_back(0.*mm);
    	temp_y.push_back(1.5*(j+1)*mm);
    	temp_x.push_back(0.*mm);
    	temp_y.push_back(-1.5*(j+1)*mm);
    }
    for(int i=0; i<1000;i++){
    	temp_x.push_back(1.5*(i+1)*mm);
    	temp_y.push_back(0.*mm);
    	temp_x.push_back(-1.5*(i+1)*mm);
    	temp_y.push_back(0.*mm);
    }
    for(int j = 0;j<1000;j++){
    	for(int i = 0;i<1000;i++){
    		temp_x.push_back(1.5*(i+1)*mm);
    		temp_y.push_back(1.5*(j+1)*mm);
    		temp_x.push_back(-1.5*(i+1)*mm);
    		temp_y.push_back(-1.5*(j+1)*mm);
    		temp_x.push_back(-1.5*(i+1)*mm);
    		temp_y.push_back(1.5*(j+1)*mm);
    		temp_x.push_back(1.5*(i+1)*mm);
    		temp_y.push_back(-1.5*(j+1)*mm);
    	}
    }
    
    G4double center_x;
    G4double center_y;
    int hi1;
    int hi2;
    int hi1reminder; // 0 -> even 1-> odd
    int hi2reminder;
    
    int numx;
    int numy;
    int reminderx;
    int remindery;
    
    G4ThreeVector v1_ = dimE->GetV1();
    G4ThreeVector v2_ = dimE->GetV2();
    G4ThreeVector v3_ = dimE->GetV3();
    G4ThreeVector v4_ = dimE->GetV4();
    
    G4double innerSide_half_ = dimE->GetInnerR_new()*tan(deltatheta_/2.);
    G4double outerSide_half_ = (dimE->GetInnerR_new()+tower_height)*tan(deltatheta_/2.);
    
    G4double theta_bc=0;
    /////////////////////
    //1. 4 types of grid coordinate(x,y) =  (o,o),(e,o),(o,e),(e,e) { e= 0, o=1 }
    
    
    int type_x_ER;
    int type_y_ER;
    
    numx = (int)(((v4_.getX()*tan(phi_unit/2.)*2)-1.*mm)/1.5*mm);
    numy = (int)((outerSide_half_*2-1.*mm)/(1.5*mm));
    reminderx = numx%2;
    remindery = numy%2;
    
    if(reminderx == 1) type_x_ER=0;
    if(reminderx == 0) type_x_ER=1;
    if(remindery == 1) type_y_ER=0;
    if(remindery == 0) type_y_ER=1;
    
    ////2. aplying the boundary conditions, reject the fibre which near the boundary of tower surface, get coord of cetres of fibre
    
    vector<G4double> center_x_ER;
    vector<G4double> center_y_ER;
    vector<G4int> bool_cfiber_ER;// 0 -> c fiber 1 s fiber
    
    theta_bc=atan(2*outerSide_half_/((v2_.getX()-v4_.getX())*tan(phi_unit/2.)));
    
    if(type_x_ER==1&&type_y_ER==1){//(o,o)
    	for(int num = 0;num<temp_x.size();num++){
    		center_x = temp_x.at(num);
    		center_y = temp_y.at(num);
    		hi1 = fabs(temp_x.at(num))/(1.5*mm);
    		hi2 = fabs(temp_y.at(num))/(1.5*mm);
            hi1reminder = hi1%2; // 0 -> even 1-> odd
            hi2reminder = hi2%2;
            
            if(((center_y-cos(theta_bc)*clad_C_rMax)>(2*2*outerSide_half_/(2*v2_.getX()*tan(phi_unit/2.)-2*v4_.getX()*tan(phi_unit/2.))*(center_x+sin(theta_bc)*clad_C_rMax-v2_.getX()*tan(phi_unit/2.))+outerSide_half_))
            	&&
            	((center_y+clad_C_rMax)<(outerSide_half_))
            	&&
            	((center_y-clad_C_rMax)>(-outerSide_half_))
            	&&
            	((center_y-cos(theta_bc)*clad_C_rMax)>(2*2*outerSide_half_/(2*v4_.getX()*tan(phi_unit/2.)-2*v2_.getX()*tan(phi_unit/2.))*(center_x-sin(theta_bc)*clad_C_rMax+v2_.getX()*tan(phi_unit/2.))+outerSide_half_
            		)))
            {
            	center_x_ER.push_back(center_x);
            	center_y_ER.push_back(center_y);
            	if(((hi1reminder==0)&&(hi2reminder==0))||((hi1reminder==1)&&(hi2reminder==1)))
            		bool_cfiber_ER.push_back(0);
            	else bool_cfiber_ER.push_back(1);
            }
            
        }
    }
    
    if(type_x_ER==0&&type_y_ER==0){//(e,e)

    	for(int num = 0;num<temp_x.size();num++){
    		center_x = temp_x.at(num)+0.75*mm;
    		center_y = temp_y.at(num)-0.75*mm;

    		hi1 = fabs(temp_x.at(num))/(1.5*mm);
    		hi2 = fabs(temp_y.at(num))/(1.5*mm);
            hi1reminder = hi1%2; // 0 -> even 1-> odd
            hi2reminder = hi2%2;
            
            if(((center_y-cos(theta_bc)*clad_C_rMax)>(2*2*outerSide_half_/(2*v2_.getX()*tan(phi_unit/2.)-2*v4_.getX()*tan(phi_unit/2.))*(center_x+sin(theta_bc)*clad_C_rMax-v2_.getX()*tan(phi_unit/2.))+outerSide_half_))
            	&&
            	((center_y+clad_C_rMax)<(outerSide_half_))
            	&&
            	((center_y-clad_C_rMax)>(-outerSide_half_))
            	&&
            	((center_y-cos(theta_bc)*clad_C_rMax)>(2*2*outerSide_half_/(2*v4_.getX()*tan(phi_unit/2.)-2*v2_.getX()*tan(phi_unit/2.))*(center_x-sin(theta_bc)*clad_C_rMax+v2_.getX()*tan(phi_unit/2.))+outerSide_half_
            		)))
            {
            	center_x_ER.push_back(center_x);
            	center_y_ER.push_back(center_y);
            	if(((hi1reminder==0)&&(hi2reminder==0))||((hi1reminder==1)&&(hi2reminder==1)))
            		bool_cfiber_ER.push_back(0);
            	else bool_cfiber_ER.push_back(1);
            }
        }
    }
    
    
    if(type_x_ER==0&&type_y_ER==1){//(e,o)

    	for(int num = 0;num<temp_x.size();num++){
    		center_x = temp_x.at(num)+0.75*mm;
    		center_y = temp_y.at(num);

    		hi1 = fabs(temp_x.at(num))/(1.5*mm);
    		hi2 = fabs(temp_y.at(num))/(1.5*mm);
            hi1reminder = hi1%2; // 0 -> even 1-> odd
            hi2reminder = hi2%2;
            
            if(((center_y-cos(theta_bc)*clad_C_rMax)>(2*2*outerSide_half_/(2*v2_.getX()*tan(phi_unit/2.)-2*v4_.getX()*tan(phi_unit/2.))*(center_x+sin(theta_bc)*clad_C_rMax-v2_.getX()*tan(phi_unit/2.))+outerSide_half_))

            	&&
            	((center_y+clad_C_rMax)<(outerSide_half_))
            	&&
            	((center_y-clad_C_rMax)>(-outerSide_half_))
            	&&
            	((center_y-cos(theta_bc)*clad_C_rMax)>(2*2*outerSide_half_/(2*v4_.getX()*tan(phi_unit/2.)-2*v2_.getX()*tan(phi_unit/2.))*(center_x-sin(theta_bc)*clad_C_rMax+v2_.getX()*tan(phi_unit/2.))+outerSide_half_
            		)))
            {
            	center_x_ER.push_back(center_x);
            	center_y_ER.push_back(center_y);
            	if(((hi1reminder==0)&&(hi2reminder==0))||((hi1reminder==1)&&(hi2reminder==1)))
            		bool_cfiber_ER.push_back(0);
            	else bool_cfiber_ER.push_back(1);
            }
            
        }
    }
    if(type_x_ER==1&&type_y_ER==0){//(o,e)

    	for(int num = 0;num<temp_x.size();num++){
    		center_x = temp_x.at(num);
    		center_y = temp_y.at(num)-0.75*mm;

    		hi1 = fabs(temp_x.at(num))/(1.5*mm);
    		hi2 = fabs(temp_y.at(num))/(1.5*mm);
            hi1reminder = hi1%2; // 0 -> even 1-> odd
            hi2reminder = hi2%2;
            
            if(((center_y-cos(theta_bc)*clad_C_rMax)>(2*2*outerSide_half_/(2*v2_.getX()*tan(phi_unit/2.)-2*v4_.getX()*tan(phi_unit/2.))*(center_x+sin(theta_bc)*clad_C_rMax-v2_.getX()*tan(phi_unit/2.))+outerSide_half_))
            	&&
            	((center_y+clad_C_rMax)<(outerSide_half_))
            	&&
            	((center_y-clad_C_rMax)>(-outerSide_half_))
            	&&
            	((center_y-cos(theta_bc)*clad_C_rMax)>(2*2*outerSide_half_/(2*v4_.getX()*tan(phi_unit/2.)-2*v2_.getX()*tan(phi_unit/2.))*(center_x-sin(theta_bc)*clad_C_rMax+v2_.getX()*tan(phi_unit/2.))+outerSide_half_
            		)))
            {
            	center_x_ER.push_back(center_x);
            	center_y_ER.push_back(center_y);
            	if(((hi1reminder==0)&&(hi2reminder==0))||((hi1reminder==1)&&(hi2reminder==1)))
            		bool_cfiber_ER.push_back(0);
            	else bool_cfiber_ER.push_back(1);
            }
            
        }
    }
    
    for(int j=0;j<center_x_ER.size();j++){
    	double z = tower->GetZHalfLength();
    	bool outside = false;
    	for(int ip = 0;ip<4;ip++){
    		TrapSidePlane plane = tower->GetSidePlane(ip);
    		double zpoint = (-plane.a*center_x_ER.at(j)-plane.b*center_y_ER.at(j)-plane.d)/plane.c;
    		outside = (tower->Inside(G4ThreeVector(center_x_ER.at(j),center_y_ER.at(j),zpoint))==kOutside);
    		if(!outside){
    			G4ThreeVector normal = tower->SurfaceNormal(G4ThreeVector(center_x_ER.at(j),center_y_ER.at(j),zpoint));
    			double angle = normal.angle(G4ThreeVector(0,0,-1));
    			double shift = fabs(clad_C_rMax/tan(0.5*M_PI-angle));
    			int length = z - zpoint - shift;
    			if(length>=1&&length<=2*z){
    				new G4PVPlacement(0,G4ThreeVector(center_x_ER.at(j),center_y_ER.at(j),z-0.5*length),
    					(bool_cfiber_ER.at(j)==0)?fiberCLog[length]:fiberSLog[length],name,towerLogicalER[i],false,j,false);
    			}
    			break;
    		}
    	}
    	if(outside) new G4PVPlacement(0,G4ThreeVector(center_x_ER.at(j),center_y_ER.at(j),0),
    		(bool_cfiber_ER.at(j)==0)?fiberCLog[2*z]:fiberSLog[2*z],name,towerLogicalER[i],false,j,false);
        if(bool_cfiber_ER.at(j)==0){//c fibre
            /*
             intersect = new G4IntersectionSolid("fiber_",fiber,tower,0,G4ThreeVector(-center_x_ER.at(j),-center_y_ER.at(j),0.));
             fiberLogical_ER[i].push_back(new G4LogicalVolume(intersect,clad_C_Material,name));
             new G4PVPlacement(0,G4ThreeVector(center_x_ER.at(j),center_y_ER.at(j),0.),fiberLogical_ER[i].at(j),name,towerLogicalER[i],false,j,checkOverlaps);
             
             intersect_= new G4IntersectionSolid("fiber_",fiberC,tower,0,G4ThreeVector(-center_x_ER.at(j),-center_y_ER.at(j),0.));
             fiberLogical_ER_[i].push_back(new G4LogicalVolume(intersect_,core_C_Material,name));
             new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),fiberLogical_ER_[i].at(j),name,fiberLogical_ER[i].at(j),false,j,checkOverlaps);
             
             
             visAttr = new G4VisAttributes(G4Colour(0.5,0.5,0.5,0.5));
             //visAttr->SetForceSolid(true);
             visAttr->SetVisibility(true);
             fiberLogical_ER[i].at(j)->SetVisAttributes(visAttr);
             fVisAttributes.push_back(visAttr);
             
             visAttr = new G4VisAttributes(G4Colour(0.,0.,1.0,0.5));
             //visAttr->SetForceSolid(true);
             visAttr->SetVisibility(true);
             fiberLogical_ER_[i].at(j)->SetVisAttributes(visAttr);
             fVisAttributes.push_back(visAttr);
             */
            //new G4PVPlacement(0,G4ThreeVector(center_x_ER.at(j),center_y_ER.at(j),0),fiberCladCLog,name,towerLogicalER[i],false,j,checkOverlaps);
        }
        
        else // s fibre
        {
            /*
             intersect = new G4IntersectionSolid("fiber_",fiber,tower,0,G4ThreeVector(-center_x_ER.at(j),-center_y_ER.at(j),0.));
             fiberLogical_ER[i].push_back(new G4LogicalVolume(intersect,clad_S_Material,name));
             new G4PVPlacement(0,G4ThreeVector(center_x_ER.at(j),center_y_ER.at(j),0.),fiberLogical_ER[i].at(j),name,towerLogicalER[i],false,j,checkOverlaps);
             
             intersect_ = new G4IntersectionSolid("fiber_",fiberS,tower,0,G4ThreeVector(-center_x_ER.at(j),-center_y_ER.at(j),0.));
             fiberLogical_ER_[i].push_back(new G4LogicalVolume(intersect_,core_S_Material,name));
             new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),fiberLogical_ER_[i].at(j),name,fiberLogical_ER[i].at(j),false,j,checkOverlaps);
             
             visAttr = new G4VisAttributes(G4Colour(0.5,0.5,0.5,0.5));
             //visAttr->SetForceSolid(true);
             visAttr->SetVisibility(true);
             fiberLogical_ER[i].at(j)->SetVisAttributes(visAttr);
             fVisAttributes.push_back(visAttr);
             
             visAttr = new G4VisAttributes(G4Colour(1.0,1.,0.,0.5));
             //visAttr->SetForceSolid(true);
             visAttr->SetVisibility(true);
             fiberLogical_ER_[i].at(j)->SetVisAttributes(visAttr);
             fVisAttributes.push_back(visAttr);
             */
            //new G4PVPlacement(0,G4ThreeVector(center_x_ER.at(j),center_y_ER.at(j),0),fiberCladSLog,name,towerLogicalER[i],false,j,checkOverlaps);
        }
    }
}

void DreamDetectorConstruction::fiberEL(G4int i,G4double deltatheta_){
	vector<G4double> temp_x;
    vector<G4double> temp_y;// vector for grid calculation
    
    temp_x.clear();
    temp_y.clear();
    temp_x.push_back(0.*mm);
    temp_y.push_back(0.*mm);
    
    for(int j = 0 ; j<1000;j++){
    	temp_x.push_back(0.*mm);
    	temp_y.push_back(1.5*(j+1)*mm);
    	temp_x.push_back(0.*mm);
    	temp_y.push_back(-1.5*(j+1)*mm);
    }
    for(int i=0; i<1000;i++){
    	temp_x.push_back(1.5*(i+1)*mm);
    	temp_y.push_back(0.*mm);
    	temp_x.push_back(-1.5*(i+1)*mm);
    	temp_y.push_back(0.*mm);
    }
    for(int j = 0;j<1000;j++){
    	for(int i = 0;i<50;i++){
    		temp_x.push_back(1.5*(i+1)*mm);
    		temp_y.push_back(1.5*(j+1)*mm);
    		temp_x.push_back(-1.5*(i+1)*mm);
    		temp_y.push_back(-1.5*(j+1)*mm);
    		temp_x.push_back(-1.5*(i+1)*mm);
    		temp_y.push_back(1.5*(j+1)*mm);
    		temp_x.push_back(1.5*(i+1)*mm);
    		temp_y.push_back(-1.5*(j+1)*mm);
    	}
    }
    
    G4double center_x;
    G4double center_y;
    int hi1;
    int hi2;
    int hi1reminder; // 0 -> even 1-> odd
    int hi2reminder;
    
    int numx;
    int numy;
    int reminderx;
    int remindery;
    
    G4ThreeVector v1_ = dimE->GetV1();
    G4ThreeVector v2_ = dimE->GetV2();
    G4ThreeVector v3_ = dimE->GetV3();
    G4ThreeVector v4_ = dimE->GetV4();
    G4double innerSide_half_ = dimE->GetInnerR_new()*tan(deltatheta_/2.);
    G4double outerSide_half_ = (dimE->GetInnerR_new()+tower_height)*tan(deltatheta_/2.);
    
    G4double theta_bc=0;
    /////////////////////
    //1. 4 types of grid coordinate(x,y) =  (o,o),(e,o),(o,e),(e,e) { e= 0, o=1 }
    
    int type_x_EL;
    int type_y_EL;
    
    numx = (int)(((v4_.getX()*tan(phi_unit/2.)*2)-1.*mm)/1.5*mm);
    numy = (int)((outerSide_half_*2-1.*mm)/(1.5*mm));
    reminderx = numx%2;
    remindery = numy%2;
    
    if(reminderx == 1) type_x_EL=0;
    if(reminderx == 0) type_x_EL=1;
    if(remindery == 1) type_y_EL=0;
    if(remindery == 0) type_y_EL=1;
    
    ////2. aplying the boundary conditions, reject the fibre which near the boundary of tower surface, get coord of cetres of fibre
    
    
    vector<G4double> center_x_EL;
    vector<G4double> center_y_EL;
    vector<G4int> bool_cfiber_EL;// 0 -> c fiber 1 s fiber
    
    theta_bc=atan(2*outerSide_half_/((v2_.getX()-v4_.getX())*tan(phi_unit/2.)));
    
    if(type_x_EL==1&&type_y_EL==1){//(o,o)
    	for(int num = 0;num<temp_x.size();num++){
    		center_x = temp_x.at(num);
    		center_y = temp_y.at(num);
    		hi1 = fabs(temp_x.at(num))/(1.5*mm);
    		hi2 = fabs(temp_y.at(num))/(1.5*mm);
            hi1reminder = hi1%2; // 0 -> even 1-> odd
            hi2reminder = hi2%2;
            
            if(((center_y+cos(theta_bc)*clad_C_rMax)<(2*2*outerSide_half_/(2*v4_.getX()*tan(phi_unit/2.)-2*v2_.getX()*tan(phi_unit/2.))*(center_x+sin(theta_bc)*clad_C_rMax-v4_.getX()*tan(phi_unit/2.))+outerSide_half_))
            	&&
            	((center_y+clad_C_rMax)<(outerSide_half_))
            	&&
            	((center_y-clad_C_rMax)>(-outerSide_half_))
            	&&
            	((center_y+cos(theta_bc)*clad_C_rMax)<(2*2*outerSide_half_/(2*v2_.getX()*tan(phi_unit/2.)-2*v4_.getX()*tan(phi_unit/2.))*(center_x-sin(theta_bc)*clad_C_rMax+v4_.getX()*tan(phi_unit/2.))+outerSide_half_
            		)))
            {
            	center_x_EL.push_back(center_x);
            	center_y_EL.push_back(center_y);
            	if(((hi1reminder==0)&&(hi2reminder==0))||((hi1reminder==1)&&(hi2reminder==1)))
            		bool_cfiber_EL.push_back(0);
            	else bool_cfiber_EL.push_back(1);
            }
            
        }
    }
    if(type_x_EL==0&&type_y_EL==0){//(e,e)

    	for(int num = 0;num<temp_x.size();num++){
    		center_x = temp_x.at(num)+0.75*mm;
    		center_y = temp_y.at(num)-0.75*mm;

    		hi1 = fabs(temp_x.at(num))/(1.5*mm);
    		hi2 = fabs(temp_y.at(num))/(1.5*mm);
            hi1reminder = hi1%2; // 0 -> even 1-> odd
            hi2reminder = hi2%2;
            
            if(((center_y+cos(theta_bc)*clad_C_rMax)<(2*2*outerSide_half_/(2*v4_.getX()*tan(phi_unit/2.)-2*v2_.getX()*tan(phi_unit/2.))*(center_x+sin(theta_bc)*clad_C_rMax-v4_.getX()*tan(phi_unit/2.))+outerSide_half_))
            	&&
            	((center_y+clad_C_rMax)<(outerSide_half_))
            	&&
            	((center_y-clad_C_rMax)>(-outerSide_half_))
            	&&
            	((center_y+cos(theta_bc)*clad_C_rMax)<(2*2*outerSide_half_/(2*v2_.getX()*tan(phi_unit/2.)-2*v4_.getX()*tan(phi_unit/2.))*(center_x-sin(theta_bc)*clad_C_rMax+v4_.getX()*tan(phi_unit/2.))+outerSide_half_
            		)))
            {
            	center_x_EL.push_back(center_x);
            	center_y_EL.push_back(center_y);
            	if(((hi1reminder==0)&&(hi2reminder==0))||((hi1reminder==1)&&(hi2reminder==1)))
            		bool_cfiber_EL.push_back(0);
            	else bool_cfiber_EL.push_back(1);
            }
        }
    }
    if(type_x_EL==0&&type_y_EL==1){//(e,o)

    	for(int num = 0;num<temp_x.size();num++){
    		center_x = temp_x.at(num)+0.75*mm;
    		center_y = temp_y.at(num);

    		hi1 = fabs(temp_x.at(num))/(1.5*mm);
    		hi2 = fabs(temp_y.at(num))/(1.5*mm);
            hi1reminder = hi1%2; // 0 -> even 1-> odd
            hi2reminder = hi2%2;
            
            if(((center_y+cos(theta_bc)*clad_C_rMax)<(2*2*outerSide_half_/(2*v4_.getX()*tan(phi_unit/2.)-2*v2_.getX()*tan(phi_unit/2.))*(center_x+sin(theta_bc)*clad_C_rMax-v4_.getX()*tan(phi_unit/2.))+outerSide_half_))
            	&&
            	((center_y+clad_C_rMax)<(outerSide_half_))
            	&&
            	((center_y-clad_C_rMax)>(-outerSide_half_))
            	&&
            	((center_y+cos(theta_bc)*clad_C_rMax)<(2*2*outerSide_half_/(2*v2_.getX()*tan(phi_unit/2.)-2*v4_.getX()*tan(phi_unit/2.))*(center_x-sin(theta_bc)*clad_C_rMax+v4_.getX()*tan(phi_unit/2.))+outerSide_half_
            		)))
            {
            	center_x_EL.push_back(center_x);
            	center_y_EL.push_back(center_y);
            	if(((hi1reminder==0)&&(hi2reminder==0))||((hi1reminder==1)&&(hi2reminder==1)))
            		bool_cfiber_EL.push_back(0);
            	else bool_cfiber_EL.push_back(1);
            }
            
        }
    }
    if(type_x_EL==1&&type_y_EL==0){//(o,e)

    	for(int num = 0;num<temp_x.size();num++){
    		center_x = temp_x.at(num);
    		center_y = temp_y.at(num)-0.75*mm;

    		hi1 = fabs(temp_x.at(num))/(1.5*mm);
    		hi2 = fabs(temp_y.at(num))/(1.5*mm);
            hi1reminder = hi1%2; // 0 -> even 1-> odd
            hi2reminder = hi2%2;
            
            if(((center_y+cos(theta_bc)*clad_C_rMax)<(2*2*outerSide_half_/(2*v4_.getX()*tan(phi_unit/2.)-2*v2_.getX()*tan(phi_unit/2.))*(center_x+sin(theta_bc)*clad_C_rMax-v4_.getX()*tan(phi_unit/2.))+outerSide_half_))
            	&&
            	((center_y+clad_C_rMax)<(outerSide_half_))
            	&&
            	((center_y-clad_C_rMax)>(-outerSide_half_))
            	&&
            	((center_y+cos(theta_bc)*clad_C_rMax)<(2*2*outerSide_half_/(2*v2_.getX()*tan(phi_unit/2.)-2*v4_.getX()*tan(phi_unit/2.))*(center_x-sin(theta_bc)*clad_C_rMax+v4_.getX()*tan(phi_unit/2.))+outerSide_half_
            		)))
            {
            	center_x_EL.push_back(center_x);
            	center_y_EL.push_back(center_y);
            	if(((hi1reminder==0)&&(hi2reminder==0))||((hi1reminder==1)&&(hi2reminder==1)))
            		bool_cfiber_EL.push_back(0);
            	else bool_cfiber_EL.push_back(1);
            }
            
        }
    }
    
    for(int j = 0; j<center_x_EL.size();j++){
    	double z = tower->GetZHalfLength();
    	bool outside = false;
    	for(int ip = 0;ip<4;ip++){
    		TrapSidePlane plane = tower->GetSidePlane(ip);
    		double zpoint = (-plane.a*center_x_EL.at(j)-plane.b*center_y_EL.at(j)-plane.d)/plane.c;
    		outside = (tower->Inside(G4ThreeVector(center_x_EL.at(j),center_y_EL.at(j),zpoint))==kOutside);
    		if(!outside){
    			G4ThreeVector normal = tower->SurfaceNormal(G4ThreeVector(center_x_EL.at(j),center_y_EL.at(j),zpoint));
    			double angle = normal.angle(G4ThreeVector(0,0,-1));
    			double shift = fabs(clad_C_rMax/tan(0.5*M_PI-angle));
    			int length = z - zpoint - shift;
    			if(length>=1&&length<=2*z){
    				new G4PVPlacement(0,G4ThreeVector(center_x_EL.at(j),center_y_EL.at(j),z-0.5*length),
    					(bool_cfiber_EL.at(j)==0)?fiberCLog[length]:fiberSLog[length],name,towerLogicalEL[i],false,j,false);
    			}
    			break;
    		}
    	}
    	if(outside) new G4PVPlacement(0,G4ThreeVector(center_x_EL.at(j),center_y_EL.at(j),0),
    		(bool_cfiber_EL.at(j)==0)?fiberCLog[2*z]:fiberSLog[2*z],name,towerLogicalEL[i],false,j,false);
        if(bool_cfiber_EL.at(j)==0){//c fibre
            /*
             intersect = new G4IntersectionSolid("fiber_",fiber,tower,0,G4ThreeVector(-center_x_EL.at(j),-center_y_EL.at(j),0.));
             fiberLogical_EL[i].push_back(new G4LogicalVolume(intersect,clad_C_Material,name));
             new G4PVPlacement(0,G4ThreeVector(center_x_EL.at(j),center_y_EL.at(j),0),fiberLogical_EL[i].at(j),name,towerLogicalEL[i],false,j,checkOverlaps);
             
             intersect_ = new G4IntersectionSolid("fiber_",fiberC,tower,0,G4ThreeVector(-center_x_EL.at(j),-center_y_EL.at(j),0.));
             fiberLogical_EL_[i].push_back(new G4LogicalVolume(intersect_,core_C_Material,name));
             new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),fiberLogical_EL_[i].at(j),name,fiberLogical_EL[i].at(j),false,j,checkOverlaps);
             
             visAttr = new G4VisAttributes(G4Colour(0.5,0.5,0.5,0.5));
             // visAttr->SetForceSolid(true);
             visAttr->SetVisibility(true);
             fiberLogical_EL[i].at(j)->SetVisAttributes(visAttr);
             fVisAttributes.push_back(visAttr);
             
             visAttr = new G4VisAttributes(G4Colour(0.,0.,1.0,0.5));
             // visAttr->SetForceSolid(true);
             visAttr->SetVisibility(true);
             fiberLogical_EL_[i].at(j)->SetVisAttributes(visAttr);
             fVisAttributes.push_back(visAttr);
             */
            //new G4PVPlacement(0,G4ThreeVector(center_x_EL.at(j),center_y_EL.at(j),0),fiberCladCLog,name,towerLogicalEL[i],false,j,checkOverlaps);
        }
        
        else // s fibre
        {/*
          intersect = new G4IntersectionSolid("fiber_",fiber,tower,0,G4ThreeVector(-center_x_EL.at(j),-center_y_EL.at(j),0.));
          fiberLogical_EL[i].push_back(new G4LogicalVolume(intersect,clad_S_Material,name));
          new G4PVPlacement(0,G4ThreeVector(center_x_EL.at(j),center_y_EL.at(j),0.),fiberLogical_EL[i].at(j),name,towerLogicalEL[i],false,j,checkOverlaps);
          
          intersect_ = new G4IntersectionSolid("fiber_",fiberS,tower,0,G4ThreeVector(-center_x_EL.at(j),-center_y_EL.at(j),0.));
          fiberLogical_EL_[i].push_back(new G4LogicalVolume(intersect_,core_S_Material,name));
          new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),fiberLogical_EL_[i].at(j),name,fiberLogical_EL[i].at(j),false,j,checkOverlaps);
          visAttr = new G4VisAttributes(G4Colour(0.5,0.5,0.5,0.5));
          // visAttr->SetForceSolid(true);
          visAttr->SetVisibility(true);
          fiberLogical_EL[i].at(j)->SetVisAttributes(visAttr);
          fVisAttributes.push_back(visAttr);
          
          visAttr = new G4VisAttributes(G4Colour(1.0,1.0,0.0,0.5));
          // visAttr->SetForceSolid(true);
          visAttr->SetVisibility(true);
          fiberLogical_EL_[i].at(j)->SetVisAttributes(visAttr);
          fVisAttributes.push_back(visAttr);
          */
            //new G4PVPlacement(0,G4ThreeVector(center_x_EL.at(j),center_y_EL.at(j),0),fiberCladSLog,name,towerLogicalEL[i],false,j,checkOverlaps);
    }
}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

