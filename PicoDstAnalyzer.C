/**
 * \brief Example of how to read a file (list of files) using StFemtoEvent classes
 *
 * RunFemtoDstAnalyzer.C is an example of reading STAR FemtoDst format.
 * One can use either FemtoDst file or a list of femtoDst files (inFile.lis or
 * inFile.list) as an input, and preform physics analysis
 *
 * \author Grigory Nigmatkulov
 * \date May 29, 2018
 */

// This is needed for calling standalone classes (not needed on RACF)
#define _VANILLA_ROOT_

// C++ headers
#include <iostream>
#include <ctime>
#include <cmath>

// ROOT headers
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TMath.h"
#include "TString.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom.h"


// FemtoDst headers
#include "StRoot/StPicoEvent/StPicoDstReader.h"
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoEpdHit.h"
/*
//We use our own Epd EP reconstruction
#include "StRoot/StEpdUtil/StEpdEpFinder.h"
#include "StRoot/StEpdUtil/StEpdEpInfo.h"
#include "StRoot/StEpdUtil/StEpdGeom.h"
*/
//#include "StRoot/StFemtoEvent/BBCTile.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

// Load libraries (for ROOT_VERSTION_CODE >= 393215)
/*
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
R__LOAD_LIBRARY(libStPicoDst.so)
#endif
*/
// Forward declarations
bool EventCut(StPicoEvent *event);
Bool_t isGoodTrack(const TVector3& vect, const Int_t& nHits,
                   const Int_t& nHitsPoss);

const Float_t electron_mass = 0.0005485799;
const Float_t pion_mass = 0.13957061;
const Float_t kaon_mass = 0.493677;
const Float_t proton_mass = 0.9382720813;

const Float_t electron_mass_sqr = 0.000000301;
const Float_t pion_mass_sqr = 0.019479955;
const Float_t kaon_mass_sqr = 0.24371698;
const Float_t proton_mass_sqr = 0.880354499;

const Double_t lambda_mass = 1.115683;
const Double_t kaon0S_mass = 0.497611;

const Double_t alpha_Lam = 0.642;
const Double_t alpha_LamBar = -0.642;

//Func
double GetPsi(int iOrd, double Qx, double Qy);
int CentralityBin(int Multiplicity);
float GetBBCTilePhi(const Int_t e_w, const Int_t iTile);
TVector3 GetRandomPointOnTile(int position, int tile, int row, int ew);
bool ScipZeroWeight(double arr[], int size);
const int nSub = 3;


// inFile - is a name of name.FemtoDst.root file or a name
//          of a name.lis(t) files that contains a list of
//          name1.FemtoDst.root files
//_________________
void PicoDstAnalyzer(const Char_t *inFile, const Char_t *outputFile,
                      const string mode) {
//    gSystem->Load("StEpdUtil");
//    gSystem->Load("StRefMultCorr");
//    gSystem->Load("libStPicoDst");
    //mods
    /*
    1)6sem //To make profiles for recentring
    2)Centred //To make profiles for flatenning
    3)Flatt //To calculate Pz

    //Does not matter at Pz calculation
    4)Weighted
    5)WeightedCentred
    6)WeightedFlatt
    */

    //Histograms with corrections
    TProfile *v1_average[9][2];
    TH1F *Coef_A_n_TH_Psi1[10][nSub];
    TH1F *Coef_B_n_TH_Psi1[10][nSub];
    TH1F *Coef_A_n_TH_Psi2[10][nSub];
    TH1F *Coef_B_n_TH_Psi2[10][nSub];
    TH1F *Qvec1Prof_TH[nSub*2];
    TH1F *Qvec2Prof_TH[nSub*2]; 
   

    if(mode == "Centred" || mode == "Flatt"){
        //Recentering
        TFile *input = new TFile("/star/u/mmorozov/14p5SecondRP/correctionsEP_output/output6sem_14p5_new.root", "read");
        for(int iSub=0; iSub!=2*nSub; iSub++){
            Qvec1Prof_TH[iSub] = (TH1F*)input->Get(Form("Qvec1Prof_%i", iSub));
            Qvec2Prof_TH[iSub] = (TH1F*)input->Get(Form("Qvec2Prof_%i", iSub));
        }

        if(mode == "Flatt"){

            //Flattening
            TFile *inputCentred = new TFile("/star/u/mmorozov/MaximTest/correctionsEP_output/outputCentred_PolarZ_70.root", "read");
            
            for(int iProf=0; iProf!=10; iProf++){
                for(int iSub=0; iSub!=nSub; iSub++){
                    Coef_A_n_TH_Psi1[iProf][iSub] = (TH1F*)inputCentred->Get(Form("Coef_A_n_Psi1_%i_%i", iProf, iSub));
                    Coef_B_n_TH_Psi1[iProf][iSub] = (TH1F*)inputCentred->Get(Form("Coef_B_n_Psi1_%i_%i", iProf, iSub));
                    Coef_A_n_TH_Psi2[iProf][iSub] = (TH1F*)inputCentred->Get(Form("Coef_A_n_Psi2_%i_%i", iProf, iSub));
                    Coef_B_n_TH_Psi2[iProf][iSub] = (TH1F*)inputCentred->Get(Form("Coef_B_n_Psi2_%i_%i", iProf, iSub));
                }
            }      
        }
    }

    TFile *inputV1 = new TFile("/star/u/mmorozov/14p5MakeCorrections/correctionsEP_output/outputV1.root", "read");
    
    for(int iCent=0; iCent!=9; iCent++){
        for(int iSub=0; iSub!=2; iSub++){
            v1_average[iCent][iSub] = (TProfile*)inputV1->Get(Form("v1_average_%i_%i", iCent, iSub));
        }
        
    }
    
    TFile *output = new TFile(outputFile, "RECREATE");
  
    //Reaction Plane  
    //0-West TPC
    //1-East TPC 
    TH1F *Qvec1Hist[9][2*nSub];
    TH1F *Psi1Hist[9][nSub];
    TH1F *Qvec2Hist[9][2*nSub];
    TH1F *Psi2Hist[9][nSub];


    for(int iHist=0; iHist!=9; iHist++){
        for(int iSub=0; iSub!=nSub; iSub++){
            Qvec1Hist[iHist][2*iSub] = new TH1F(Form("Qvec1Hist_%i_%i", iHist, 2*iSub), "", 500, -1, 1);
            Qvec1Hist[iHist][2*iSub+1] = new TH1F(Form("Qvec1Hist_%i_%i", iHist, 2*iSub+1), "", 500, -1, 1);
            Qvec2Hist[iHist][2*iSub] = new TH1F(Form("Qvec2Hist_%i_%i", iHist, 2*iSub), "", 500, -1, 1);
            Qvec2Hist[iHist][2*iSub+1] = new TH1F(Form("Qvec2Hist_%i_%i", iHist, 2*iSub+1), "", 500, -1, 1);
            Psi1Hist[iHist][iSub] = new TH1F(Form("Psi1Hist_%i_%i", iHist, iSub), "", 140, 0, 7);
            Psi2Hist[iHist][iSub] = new TH1F(Form("Psi2Hist_%i_%i", iHist, iSub), "", 70, 0, 3.5);
        }
    }

    //Profiles for centering
    TProfile *Qvec1Prof[2*nSub];
    TProfile *Qvec2Prof[2*nSub];
    for(int iSub=0; iSub!=2*nSub; iSub++){
        Qvec1Prof[iSub] = new TProfile(Form("Qvec1Prof_%i", iSub), "", 9, 0, 9);
        Qvec2Prof[iSub] = new TProfile(Form("Qvec2Prof_%i", iSub), "", 9, 0, 9);
    }
   
    
    //Profiles for flattening
    //A_n = (-2/n)*<sin(nPsi)>
    //B_n = (2/n)*<cos(nPsi)>
    TProfile *Coef_A_n_Psi1[10][nSub];
    TProfile *Coef_B_n_Psi1[10][nSub];
    TProfile *Coef_A_n_Psi2[10][nSub];
    TProfile *Coef_B_n_Psi2[10][nSub];


    for(int iProf=0; iProf!=10; iProf++){
        for(int iSub=0; iSub!=nSub; iSub++){
            Coef_A_n_Psi1[iProf][iSub] = new TProfile(Form("Coef_A_n_Psi1_%i_%i", iProf, iSub), "", 9, 0, 9);
            Coef_B_n_Psi1[iProf][iSub] = new TProfile(Form("Coef_B_n_Psi1_%i_%i", iProf, iSub), "", 9, 0, 9);
            Coef_A_n_Psi2[iProf][iSub] = new TProfile(Form("Coef_A_n_Psi2_%i_%i", iProf, iSub), "", 9, 0, 9);
            Coef_B_n_Psi2[iProf][iSub] = new TProfile(Form("Coef_B_n_Psi1_%i_%i", iProf, iSub), "", 9, 0, 9);
	    }
    }

    //Resolution
    TProfile *CosOfDiff_1 = new TProfile("CosOfDiff_1", "CosOfDiff for Psi_1", 9, 0, 9);
    TProfile *CosOfDiff_2 = new TProfile("CosOfDiff_2", "CosOfDiff for Psi_2", 9, 0, 9);
        

    //Variables
    //Reaction Plane
    double QWeight_1[nSub] = {}, QWeight_2[nSub] = {};
    double Qvec_1[2*nSub] = {}, Qvec_2[2*nSub] = {};
    double Psi1[nSub] = {}, Psi2[nSub] = {};
    double deltaPsi1[nSub] = {}, deltaPsi2[nSub] = {};
    double phi, w;
    int PP, TT, EW, iSide = 0, row;


    //corners EPD
    double phiCenter[12][31][2];
    double deltaPhi = (30.0/180.0)*TMath::Pi(); 
    int ew = 0;//east

    for(int pp=1; pp<13; pp++){
        double phiPpCenter = TMath::Pi()/2.0 - (pp-0.5)*deltaPhi;
    	if(phiPpCenter<0.0) phiPpCenter += 2.0*TMath::Pi();
    	phiCenter[pp-1][0][ew] = phiPpCenter;

    	for(int tt=2; tt<32; tt+=2){
      	    phiCenter[pp-1][tt-1][ew] = phiPpCenter - deltaPhi/4.0;
    	}

    	for(int tt=3; tt<32; tt+=2){
      	    phiCenter[pp-1][tt-1][ew] = phiPpCenter + deltaPhi/4.0;
    	}
    }
  
    ew = 1;//west 5.89049
	

    for(int pp=1; pp<13; pp++){
    	double phiPpCenter = TMath::Pi()/2.0 + (pp-0.5)*deltaPhi;
    	if(phiPpCenter>2.0*TMath::Pi()) phiPpCenter -= 2.0*TMath::Pi();
    	phiCenter[pp-1][0][ew] = phiPpCenter;

    	for(int tt=2; tt<32; tt+=2){
     	    phiCenter[pp-1][tt-1][ew] = phiPpCenter + deltaPhi/4.0;
    	}
    	for(int tt=3; tt<32; tt+=2){
      	    phiCenter[pp-1][tt-1][ew] = phiPpCenter - deltaPhi/4.0;
    	}

    }
   

    
    //-------------------------------Start analizing--------------------------------------------------
    std::cout << "Hi! Lets do some physics, Master!" << std::endl;

    Int_t runIdBins = 2000000;
    Int_t runIdRange[2] = { 11000000, 13000000 };
/*
    gSystem->Load("./libStFemtoDst.so");
    #if ROOT_VERSION_CODE < ROOT_VERSION(6,0,0)
        gSystem->Load("./libStFemtoDst.so");
    #endif
*/
    StPicoDstReader* femtoReader = new StPicoDstReader(inFile);
    femtoReader->Init();

    //Long64_t events2read = femtoReader->chain()->GetEntries();

    // This is a way if you want to spead up IO
    std::cout << "Explicit read status for some branches" << std::endl;
    femtoReader->SetStatus("*",0);
    femtoReader->SetStatus("Event",1);
    femtoReader->SetStatus("Track",1);
    //femtoReader->SetStatus("KFP", 1);
    femtoReader->SetStatus("EpdHit", 1);
    std::cout << "Status has been set" << std::endl;

    std::cout << "Now I know what to read, Master!" << std::endl;

    if( !femtoReader->chain() ) {
        std::cout << "No chain has been found." << std::endl;
    }
    Long64_t eventsInTree = femtoReader->tree()->GetEntries();
    std::cout << "eventsInTree: "  << eventsInTree << std::endl;
    Long64_t events2read = femtoReader->chain()->GetEntries();

    std::cout << "Number of events to read: " << events2read << std::endl;  

    
    // Loop over events
    for(Long64_t iEvent=0; iEvent<events2read; iEvent++) {

        std::cout << "Working on event #[" << (iEvent+1)
                << "/" << events2read << "]" << mode << std::endl;


        Bool_t readEvent = femtoReader->readPicoEvent(iEvent);
        if( !readEvent ) {
            std::cout << "Something went wrong, Master! Nothing to analyze..." << std::endl;
            break;
        }

        // Retrieve femtoDst
        StPicoDst *dst = femtoReader->picoDst();

        // Retrieve event information
        StPicoEvent *event = dst->event();
        if( !event ) {
            std::cout << "Something went wrong, Master! Event is hiding from me..." << std::endl;
            break;
        }


	    // Simple event cut
        if (EventCut(event) == false ) continue;
	
        
        //Centrality
        //TODO: Заменить поиск центральности и плохих событий на поиск по StRefMultCorr (см. /star/u/alpatov/EPCOnstructor/PicoFourteenAnalyzer.C)
        StRefMultCorr* refmultCorrUtil = CentralityMaker::instance()->getRefMultCorr() ;
        refmultCorrUtil -> init(event->runId());
        refmultCorrUtil -> initEvent(event->refMult(), event->primaryVertex().z(), event->ZDCx());
        Bool_t isBadRun = refmultCorrUtil-> isBadRun(event->runId()); //reject bad runs
        Bool_t isPileUpEvt = !refmultCorrUtil->passnTofMatchRefmultCut(1.*event->refMult(), 1.*event->nBTOFMatch()); //reject pileup events

        int cent = refmultCorrUtil->getCentralityBin9() ;
        if (cent < 0 || isBadRun || isPileUpEvt) continue; 

        if(cent < 0) continue;

        for(int iSub=0; iSub!=nSub; iSub++){
            Qvec_1[2*iSub] = 0.;
            Qvec_1[2*iSub+1] = 0.;
            QWeight_1[iSub] = 0.;
            Psi1[iSub] = 0.;
            Qvec_2[2*iSub] = 0.;
            Qvec_2[2*iSub+1] = 0.;
            QWeight_2[iSub] = 0.;
            Psi2[iSub] = 0.;
        }

        
        
        //-------------------------------Track analysis--------------------------------------------------
        Int_t nTracks = dst->numberOfTracks();

        // Track loop
        for(Int_t iTrk=0; iTrk<nTracks; iTrk++) {
	   
        

            // Retrieve i-th femto track
            StPicoTrack *femtoTrack = dst->track(iTrk);

            if (!femtoTrack) continue;
            //std::cout << "Track #[" << (iTrk+1) << "/" << nTracks << "]"  << std::endl;
	        double eta = femtoTrack->pMom().Eta();

	        if (femtoTrack->nHits() < 15 || femtoTrack->pPt() < 0.1 || femtoTrack->pPt() > 2. || fabs(femtoTrack->pMom().Eta()) > 1.5 || fabs(eta)<0.1) continue;
            if ((Float_t)femtoTrack->nHits()/femtoTrack->nHitsPoss() < 0.52 ) continue;

            double wEff = femtoTrack->pPt();
            double phi = femtoTrack->pMom().Phi(); //femtoTrack->isPrimary() ? femtoTrack->pMom().Phi() : 0.;
            if(eta>0) iSide = 0;//West
            else if(eta<0) iSide = 1;//East

            Qvec_2[2*iSide] += wEff*TMath::Cos(2*phi);
            Qvec_2[2*iSide+1] += wEff*TMath::Sin(2*phi);
            QWeight_2[iSide] += wEff;

            
        }//for(Int_t iTrk=0; iTrk<nTracks; iTrk++)

		//.........................................................EpdHit analysis..................................................................
        int nEpdHits = dst->numberOfEpdHits();

        //EpdHit Loop
        for(int iEpd=0; iEpd<nEpdHits; iEpd++){
            //Retrieve i-th EpdHit
            StPicoEpdHit *femtoEpdHit = dst->epdHit(iEpd);

            if(!femtoEpdHit->isGood()) continue;
    
            PP = femtoEpdHit->position();
            TT = femtoEpdHit->tile();
            EW = femtoEpdHit->side();
            w = femtoEpdHit->nMIP();
            row = femtoEpdHit->row();

            double wEff_EPD = w;
            if(wEff_EPD<0.3) continue;
            if(wEff_EPD>3) wEff_EPD = 3;

            if(EW == 1){
                wEff_EPD *= fabs(v1_average[cent][0]->GetBinContent(row));
                Qvec_1[0] += wEff_EPD*TMath::Cos(phiCenter[PP-1][TT-1][EW]);
                Qvec_1[1] += wEff_EPD*TMath::Sin(phiCenter[PP-1][TT-1][EW]);
                QWeight_1[0] += wEff_EPD;
            }

            if(EW==-1){
                wEff_EPD *= fabs(v1_average[cent][1]->GetBinContent(row));
                Qvec_1[2] += wEff_EPD*TMath::Cos(phiCenter[PP-1][TT-1][0]);
                Qvec_1[3] += wEff_EPD*TMath::Sin(phiCenter[PP-1][TT-1][0]);
                QWeight_1[1] += wEff_EPD;
            }
        }

            
        //-------------------------------reaction plane--------------------------------------------------
        if(ScipZeroWeight(QWeight_1, nSub) == false) continue;
        if(ScipZeroWeight(QWeight_2, nSub) == false) continue;
            
        //Get Q vectors
        bool check = true;            
        for(int iSub=0; iSub!=nSub; iSub++){
            Qvec_1[2*iSub] = Qvec_1[2*iSub]/QWeight_1[iSub];
            Qvec_1[2*iSub+1] = Qvec_1[2*iSub+1]/QWeight_1[iSub];
            if(fabs(Qvec_1[2*iSub])>999 || fabs(Qvec_1[2*iSub+1])>999) check = false;

            Qvec_2[2*iSub] = Qvec_2[2*iSub]/QWeight_2[iSub];
            Qvec_2[2*iSub+1] = Qvec_2[2*iSub+1]/QWeight_2[iSub];
            if(fabs(Qvec_2[2*iSub])>999 || fabs(Qvec_2[2*iSub+1])>999) check = false;
        }
        if(!check) continue;
        Qvec_1[4] = Qvec_1[0] + Qvec_1[2];
        Qvec_1[5] = Qvec_1[1] + Qvec_1[3];
        Qvec_2[4] = Qvec_2[0] + Qvec_2[2];
        Qvec_2[5] = Qvec_2[1] + Qvec_2[3];

        
        //Centring collibration
        if(mode == "Centred" || mode == "Flatt"){
            //Recentering       
            for(int iSub = 0; iSub!=2*nSub; iSub++){
                Qvec_1[iSub] -= Qvec1Prof_TH[iSub]->GetBinContent(cent+1);
                Qvec_2[iSub] -= Qvec2Prof_TH[iSub]->GetBinContent(cent+1);
            }
        }

        
        //Get Psi
        for(int iSub=0; iSub!=nSub; iSub++){
            Psi1[iSub] = GetPsi(1, Qvec_1[2*iSub], Qvec_1[2*iSub+1]);
            if(Psi1[iSub] > 7) check = false;

            Psi2[iSub] = GetPsi(2, Qvec_2[2*iSub], Qvec_2[2*iSub+1]);
            if(Psi2[iSub] > 3.5) check = false;
        }
        if(!check) continue;
                            

        //Flattening collibration
        if(mode == "Flatt"){
            //Flattaning
            for(int iProf=0; iProf!=10; iProf++){
                for(int iSub=0; iSub!=nSub; iSub++){
                    deltaPsi1[iSub] += Coef_A_n_TH_Psi1[iProf][iSub]->GetBinContent(cent+1)*TMath::Cos((iProf+1)*Psi1[iSub]) + \
                                        Coef_B_n_TH_Psi1[iProf][iSub]->GetBinContent(cent+1)*TMath::Sin((iProf+1)*Psi1[iSub]);
                    deltaPsi2[iSub] += Coef_A_n_TH_Psi2[iProf][iSub]->GetBinContent(cent+1)*TMath::Cos((iProf+1)*2*Psi2[iSub]) + \
                                        Coef_B_n_TH_Psi2[iProf][iSub]->GetBinContent(cent+1)*TMath::Sin((iProf+1)*2*Psi2[iSub]);
                }
            }
        
            for(int iSub=0; iSub!=nSub; iSub++){
                Psi1[iSub] += deltaPsi1[iSub];
                while(Psi1[iSub]>2*TMath::Pi()) Psi2[iSub] -= 2*TMath::Pi();
                while(Psi1[iSub]<0.0) Psi2[iSub] += 2*TMath::Pi();
                deltaPsi1[iSub] = 0.0;
                
                Psi2[iSub] += deltaPsi2[iSub]/2.;
                while(Psi2[iSub]>TMath::Pi()) Psi2[iSub] -= TMath::Pi();
                while(Psi2[iSub]<0.0) Psi2[iSub] += TMath::Pi();
                deltaPsi2[iSub] = 0.0;
            }
        }    
                        
                                    
        //Fill histograms
        for(int iSub = 0; iSub!=nSub; iSub++){
            Qvec1Hist[cent][2*iSub]->Fill(Qvec_1[2*iSub]);
            Qvec1Hist[cent][2*iSub+1]->Fill(Qvec_1[2*iSub+1]);
            Psi1Hist[cent][iSub]->Fill(Psi1[iSub]);

            Qvec2Hist[cent][2*iSub]->Fill(Qvec_2[2*iSub]);
            Qvec2Hist[cent][2*iSub+1]->Fill(Qvec_2[2*iSub+1]);
            Psi2Hist[cent][iSub]->Fill(Psi2[iSub]);
        }            


        if(mode == "6sem"){
            //Profiles for recentering
            for(int iSub=0; iSub!=nSub*2; iSub++){
                Qvec1Prof[iSub]->Fill(cent, Qvec_1[iSub]);
                Qvec2Prof[iSub]->Fill(cent, Qvec_2[iSub]);
            }
        }    


        if(mode == "Centred"){
            //Profiles for flattening
            for(int iProf=0; iProf!=10; iProf++){
                for(int iSub=0; iSub!=nSub; iSub++){
                    Coef_A_n_Psi1[iProf][iSub]->Fill(cent, (-2.0/(iProf+1))*TMath::Sin((iProf+1)*Psi1[iSub]));
                    Coef_B_n_Psi1[iProf][iSub]->Fill(cent, (2.0/(iProf+1))*TMath::Cos((iProf+1)*Psi1[iSub]));

                    Coef_A_n_Psi2[iProf][iSub]->Fill(cent, (-2.0/(iProf+1))*TMath::Sin((iProf+1)*2*Psi2[iSub]));
                    Coef_B_n_Psi2[iProf][iSub]->Fill(cent, (2.0/(iProf+1))*TMath::Cos((iProf+1)*2*Psi2[iSub]));
                }
            }        
        }

        CosOfDiff_1->Fill(cent, TMath::Cos(Psi1[1] - Psi1[0]));
        CosOfDiff_2->Fill(cent, TMath::Cos(2*(Psi2[1] - Psi2[0])));
        
        

    }//for(Long64_t iEvent=0; iEvent<events2read; iEvent++)


    //EP Distributions
    for(int iCent=0; iCent!=9; iCent++){
	    Qvec1Hist[iCent][0]->SetTitle(Form("Qx West EPD first harm (%i)", iCent));
        Qvec1Hist[iCent][1]->SetTitle(Form("Qy West EPD first harm (%i)", iCent));
        Qvec1Hist[iCent][2]->SetTitle(Form("Qx East EPD first harm (%i)", iCent));
        Qvec1Hist[iCent][3]->SetTitle(Form("Qy East EPD first harm (%i)", iCent));
        Qvec1Hist[iCent][4]->SetTitle(Form("Qx Comb EPD first harm (%i)", iCent));
        Qvec1Hist[iCent][5]->SetTitle(Form("Qy Comb EPD first harm (%i)", iCent));

        Qvec2Hist[iCent][0]->SetTitle(Form("Qx West TPC second harm (%i)", iCent));
        Qvec2Hist[iCent][1]->SetTitle(Form("Qy West TPC second harm (%i)", iCent));
        Qvec2Hist[iCent][2]->SetTitle(Form("Qx East TPC second harm (%i)", iCent));
        Qvec2Hist[iCent][3]->SetTitle(Form("Qy East TPC second harm (%i)", iCent));
        Qvec2Hist[iCent][4]->SetTitle(Form("Qx Comb TPC second harm (%i)", iCent));
        Qvec2Hist[iCent][5]->SetTitle(Form("Qy Comb TPC second harm (%i)", iCent));

        Psi1Hist[iCent][0]->SetTitle(Form("West Psi_1 EPD (%i)", iCent));
        Psi1Hist[iCent][1]->SetTitle(Form("East Psi_1 EPD (%i)", iCent));
        Psi1Hist[iCent][2]->SetTitle(Form("Comb Psi_1 EPD (%i)", iCent));
        
        Psi2Hist[iCent][0]->SetTitle(Form("West Psi_2 TPC (%i)", iCent));
        Psi2Hist[iCent][1]->SetTitle(Form("East Psi_2 TPC (%i)", iCent));
        Psi2Hist[iCent][2]->SetTitle(Form("Comb Psi_2 TPC (%i)", iCent));

    }
    Qvec1Prof[0]->SetTitle("Profile Qx West EPD first harm");
    Qvec1Prof[1]->SetTitle("Profile Qy West EPD first harm");
    Qvec1Prof[2]->SetTitle("Profile Qx East EPD first harm");
    Qvec1Prof[3]->SetTitle("Profile Qy East EPD first harm");
    Qvec1Prof[4]->SetTitle("Profile Qx Comb EPD first harm");
    Qvec1Prof[5]->SetTitle("Profile Qy Comb EPD first harm");

    Qvec2Prof[0]->SetTitle("Profile Qx West TPC second harm");
    Qvec2Prof[1]->SetTitle("Profile Qy West TPC second harm");
    Qvec2Prof[2]->SetTitle("Profile Qx East TPC second harm");
    Qvec2Prof[3]->SetTitle("Profile Qy East TPC second harm");
    Qvec2Prof[4]->SetTitle("Profile Qx Comb TPC second harm");
    Qvec2Prof[5]->SetTitle("Profile Qy Comb TPC second harm");

    for(int i=0; i!=10; i++){
	    Coef_A_n_Psi1[i][0]->SetTitle(Form("A_n_%i West Psi_1 EPD", i));
	    Coef_B_n_Psi1[i][0]->SetTitle(Form("B_n_%i West Psi_1 EPD", i));
	    Coef_A_n_Psi1[i][1]->SetTitle(Form("A_n_%i East Psi_1 EPD", i));
        Coef_B_n_Psi1[i][1]->SetTitle(Form("B_n_%i East Psi_1 EPD", i));
	    Coef_A_n_Psi1[i][2]->SetTitle(Form("A_n_%i Comb Psi_1 EPD", i));
        Coef_B_n_Psi1[i][2]->SetTitle(Form("B_n_%i Comb Psi_1 EPD", i));
        
        Coef_A_n_Psi2[i][0]->SetTitle(Form("A_n_%i West Psi_2 TPC", i));
	    Coef_B_n_Psi2[i][0]->SetTitle(Form("B_n_%i West Psi_2 TPC", i));
	    Coef_A_n_Psi2[i][1]->SetTitle(Form("A_n_%i East Psi_2 TPC", i));
        Coef_B_n_Psi2[i][1]->SetTitle(Form("B_n_%i East Psi_2 TPC", i));
	    Coef_A_n_Psi2[i][2]->SetTitle(Form("A_n_%i Comb Psi_2 TPC", i));
        Coef_B_n_Psi2[i][2]->SetTitle(Form("B_n_%i Comb Psi_2 TPC", i));
    }

    femtoReader->Finish();

    output->Write();

    std::cout << "I'm done with analysis. We'll have a Nobel Prize, Master!"
        << std::endl;

    
}

bool EventCut(StPicoEvent *event)
{
  bool cut = true;
  if (!event->isTrigger(650000) && !event->isTrigger(650001) && !event->isTrigger(650002) && !event->isTrigger(650003) &&
      !event->isTrigger(650007) && !event->isTrigger(650004) && !event->isTrigger(650005) && !event->isTrigger(650006) &&
      !event->isTrigger(650009) ) cut = false;
  if (event->primaryVertex().Z() < -145. || event->primaryVertex().Z() > 145. || ( pow(event->primaryVertex().X(),2) + pow(event->primaryVertex().Y(),2) > 4)) cut = false;
//  if (!GoodRun(event))  cut = false;
  return cut;


}

//________________
Bool_t isGoodTrack(const TVector3& mom, const Int_t& nHits,
                   const Int_t& nHitsPoss) {
  return ( mom.Mag() >= 0.15 &&
           mom.Mag() <= 3 &&
           nHits >= 14 &&
           (Float_t)nHits/nHitsPoss >= 0.51 );
}

//Get corner EPD
double GetPsi(int iOrd, double Qx, double Qy){
    double Psi;
    Psi = atan2(Qy, Qx)/iOrd;
    if(Psi<0.) Psi += 2*TMath::Pi() / iOrd;
    return Psi;

}

/*double GetPsiWest(double Qx, double Qy){
  double PsiWest;
  if(Qx>0){
    PsiWest = TMath::ATan(Qy/Qx);
    if(PsiWest<0.0) PsiWest += 2.0*TMath::Pi();
  }else{
    PsiWest = TMath::Pi() + TMath::ATan(Qy/Qx);
  }
  return PsiWest;    
}*/

//Centrality
int CentralityBin(int Multiplicity) {
    int CentId;
    if (Multiplicity >= 290) CentId = 8; //0-5%
    else if (Multiplicity >= 233) CentId = 7; //5-10%
    else if (Multiplicity >= 150) CentId = 6; //10-20%
    else if (Multiplicity >= 94) CentId = 5; //20-30%
    else if (Multiplicity >= 57) CentId = 4; //30-40%
    else if (Multiplicity >= 32) CentId = 3; //40-50%
    else if (Multiplicity >= 17) CentId = 2; //50-60%
    else if (Multiplicity >= 9) CentId = 1; //60-70%
    else if (Multiplicity >= 4) CentId = 0; //70-80%
    else CentId = -1;
    return CentId;
}

float GetBBCTilePhi(const Int_t e_w, const Int_t iTile){
    Double_t pi = TMath::Pi();
    Double_t phi_div  = pi/6.0;
    float bbc_phi = phi_div;
    switch(iTile) {
    case 0: bbc_phi=3*phi_div;
    break;
    case 1: bbc_phi=phi_div;
    break;
    case 2: bbc_phi=-1*phi_div;
    break;
    case 3: bbc_phi=-3*phi_div;
    break;
    case 4: bbc_phi=-5*phi_div;
    break;
    case 5: bbc_phi=5*phi_div;
    break;
    case 6: bbc_phi= 3*phi_div;
    break;
    case 7: bbc_phi=2*phi_div;
    break;
    case 8: bbc_phi=phi_div;
    break;
    case 9: bbc_phi=0.;
    break;
    case 10: bbc_phi=-phi_div;
    break;
    case 11: bbc_phi=-2*phi_div;
    break;
    case 12: bbc_phi=-3*phi_div;
    break;
    case 13: bbc_phi=-4*phi_div;
    break;
    case 14: bbc_phi=-5*phi_div;
    break;
    case 15: bbc_phi=pi;
    break;
    case 16: bbc_phi=5*phi_div;
    break;
    case 17: bbc_phi=4*phi_div;
    break;
    }
    if(e_w==0){if (bbc_phi > -0.001){ bbc_phi = pi-bbc_phi;}
                    else {bbc_phi= -pi-bbc_phi;}
                }
    if(bbc_phi<0.0) bbc_phi +=2*pi;
    if(bbc_phi>2*pi) bbc_phi -=2*pi;
    return bbc_phi;
}

//Get random point on a tile
TVector3 GetRandomPointOnTile(int position, int tile, int row, int ew){
    double phiCenter[12][31][2];
    double deltaPhi = (30.0/180.0)*TMath::Pi(); 
    int EW = 0;//east

    for(int pp=1; pp<13; pp++){
        double phiPpCenter = TMath::Pi()/2.0 - (pp-0.5)*deltaPhi;
        if(phiPpCenter<0.0) phiPpCenter += 2.0*TMath::Pi();
        phiCenter[pp-1][0][EW] = phiPpCenter;

        for(int tt=2; tt<32; tt+=2){
        phiCenter[pp-1][tt-1][EW] = phiPpCenter - deltaPhi/4.0;      
        }

        for(int tt=3; tt<32; tt+=2){
        phiCenter[pp-1][tt-1][EW] = phiPpCenter + deltaPhi/4.0;
        }
        
    }
    EW = 1;//west 5.89049

    for(int pp=1; pp<13; pp++){
        double phiPpCenter = TMath::Pi()/2.0 + (pp-0.5)*deltaPhi;
        if(phiPpCenter>2.0*TMath::Pi()) phiPpCenter -= 2.0*TMath::Pi();
        phiCenter[pp-1][0][EW] = phiPpCenter;

        for(int tt=2; tt<32; tt+=2){
        phiCenter[pp-1][tt-1][EW] = phiPpCenter + deltaPhi/4.0;
        }
        for(int tt=3; tt<32; tt+=2){
        phiCenter[pp-1][tt-1][EW] = phiPpCenter - deltaPhi/4.0;
        }

    }

    double Row[17] = {4.6, 9.0, 13.4, 17.8, 23.33, 28.86, 34.39, 39.92, 45.45, 50.98, 56.51, 62.05, 67.58, 73.11, 78.64, 84.17, 89.70};
    TRandom *rand = new TRandom();
    //Random value (0,1)
    double randValue = rand->Rndm();

    double randRow = (Row[row] - Row[row-1])*randValue + Row[row-1];
    if(ew == -1) ew = 0;
    double randDeltaPhi = (deltaPhi/2.0)*randValue;
    if(randDeltaPhi > deltaPhi/4.0) randDeltaPhi = -randDeltaPhi + deltaPhi/4.0;
    double randPhi = phiCenter[position - 1][tile - 1][ew] + randDeltaPhi;

    double x = randRow*TMath::Cos(randPhi);
    double y = randRow*TMath::Sin(randPhi);
    double z = 375;
    if(ew == 0) z = -375;
    TVector3 RandPoint(x,y,z);

    return RandPoint;
}


bool ScipZeroWeight(double arr[], int size) {
    for (size_t i = 0; i < size; ++i) {
        if (arr[i] == 0) {
            return false;
        }
    }
    return true;
}
