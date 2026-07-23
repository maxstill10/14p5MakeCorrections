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
double GetPsi(double Qx, double Qy);
int CentralityBin(int Multiplicity);
float GetBBCTilePhi(const Int_t e_w, const Int_t iTile);
TVector3 GetRandomPointOnTile(int position, int tile, int row, int ew);
bool ScipZeroWeight(double arr[], int size);
const int nSub = 7;


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
    TH1F *Coef_A_n_TH[10][nSub];
    TH1F *Coef_B_n_TH[10][nSub];
    TH1F *QvecProf_TH[nSub*2]; 
   

    if(mode == "Centred" || mode == "Flatt"){
        //Recentering
        TFile *input = new TFile("/star/u/mmorozov/14p5SecondRP/correctionsEP_output/output6sem_14p5_new.root", "read");
        for(int iSub=0; iSub!=2*nSub; iSub++){
	    QvecProf_TH[iSub] = (TH1F*)input->Get(Form("QvecProf_%i", iSub));
	}

        if(mode == "Flatt"){

            //Flattening
            TFile *inputCentred = new TFile("/star/u/mmorozov/MaximTest/correctionsEP_output/outputCentred_PolarZ_70.root", "read");
            
            for(int iProf=0; iProf!=10; iProf++){
		for(int iSub=0; iSub!=nSub; iSub++){
		    Coef_A_n_TH[iProf][iSub] = (TH1F*)inputCentred->Get(Form("Coef_A_n_%i_%i", iProf, iSub));
                    Coef_B_n_TH[iProf][iSub] = (TH1F*)inputCentred->Get(Form("Coef_B_n_%i_%i", iProf, iSub));
		}
            }      
        }
    }


    
    TFile *output = new TFile(outputFile, "RECREATE");
  
  
    //General distribution
    TH2F *DEdx = new TH2F("DEdx", "dEdx", 300, -3, 3, 300, -2, 15); 
    TH2F *DEdxProt = new TH2F("DEdxProt", "dEdx of p", 300, -3, 3, 300, -2, 10);
    TH2F *DEdxPion = new TH2F("DEdxPion", "dEdx of #pi", 300, -3, 3, 300, -2, 10); 
    TH2F *NsigmakaonMagn = new TH2F("NsigmakaonMagn", "NsigmakaonMagn", 300, -3, 3, 300, -10, 10);
    TH2F *NsigmapionMagn = new TH2F("NsigmapionMagn", "NsigmapionMagn", 600, -3, 3, 600, -10, 10);
    TH2F *NsigmaprotMagn = new TH2F("NsigmaprotMagn", "NsigmaprotMagn", 600, -3, 3, 600, -10, 10);
    TH2F *NsigmakaonTrans = new TH2F("NsigmakaonTrans", "NsigmakaonTrans", 300, -3, 3, 300, -10, 10);
    TH2F *NsigmapionTrans = new TH2F("NsigmapionTrans", "NsigmapionTrans", 300, -3, 3, 300, -10, 10);
    TH2F *NsigmaprotTrans = new TH2F("NsigmaprotTrans", "NsigmaprotTrans", 300, -3, 3, 300, -10, 10);
    TH2F *PrimeVertXvsY = new TH2F("PrimeVertXvsY", "Primery vertex X Vs Y", 500, -5, 5, 500, -5, 5);
    TH1F *PrimeVertZ = new TH1F("PrimeVertZ", "Primerty vertex Z", 100, -100, 100);
    TH1F *NHits = new TH1F("NHits", "nHits", 100, 10, 90);
    TH1F *Pseudorap = new TH1F("Pseudorap", "Pseudorepidity", 100, -10, 10);
    TH1F *Mom = new TH1F("Mom", "Momentum", 100, 0, 5);
    //TH2F *BettaVsMom = new TH2F("BettaVsMom", "Betta vs  Momentum", 100, -1, 1, 100, 0, 5);



    //Reaction Plane  
    //0-West TPC
    //1-East TPC
    //2-Comb TPC
    //3-E+W TPC
    //4-West EPD
    //5-East EPD
    //6-E+W EPD   
    TH1F *QvecHist[9][2*nSub];
    TH1F *Psi2Hist[9][nSub];
    TH1F *Psi2Full_EW[9][6];


    for(int iHist=0; iHist!=9; iHist++){
	for(int iSub=0; iSub!=nSub; iSub++){
	    QvecHist[iHist][2*iSub] = new TH1F(Form("QvecHist_%i_%i", iHist, 2*iSub), "", 500, -1, 1);
	    QvecHist[iHist][2*iSub+1] = new TH1F(Form("QvecHist_%i_%i", iHist, 2*iSub+1), "", 500, -1, 1);
	    Psi2Hist[iHist][iSub] = new TH1F(Form("Psi2Hist_%i_%i", iHist, iSub), "", 70, 0, 3.5);
	}
	//Psi[i] - Psi[j]
	for(int i=0; i!=6; i++){
	    Psi2Full_EW[iHist][i] = new TH1F(Form("Psi2Full_EW_%i_%i", iHist, i), "", 140, -3.5, 3.5);
	}
    }

    //Profiles for centering
    TProfile *QvecProf[2*nSub];
    for(int iSub=0; iSub!=2*nSub; iSub++){
        QvecProf[iSub] = new TProfile(Form("QvecProf_%i", iSub), "", 9, 0, 9);
    }
   
    
    //Profiles for flattening
    //A_n = (-2/n)*<sin(nPsi)>
    //B_n = (2/n)*<cos(nPsi)>
    TProfile *Coef_A_n[10][nSub];
    TProfile *Coef_B_n[10][nSub];


    for(int iProf=0; iProf!=10; iProf++){
        for(int iSub=0; iSub!=nSub; iSub++){
	    Coef_A_n[iProf][iSub] = new TProfile(Form("Coef_A_n_%i_%i", iProf, iSub), "", 9, 0, 9);
	    Coef_B_n[iProf][iSub] = new TProfile(Form("Coef_B_n_%i_%i", iProf, iSub), "", 9, 0, 9);
	}
    }

    //Resolution
    TProfile *CosOfDiff = new TProfile("CosOfDiff", "CosOfDiff", 9, 0, 9);
    TProfile *CosOfDiff_EPD = new TProfile("CosOfDiff_EPD", "Comb - EPD_West", 9, 0, 9);
    TProfile *CosOfDiff_Comb = new TProfile("CosOfDiff_Comb", "Comb - EPD_East", 9, 0, 9);
        

    //Variables
    //Reaction Plane
    double QWeight[nSub] = {};
    double Qvec[2*nSub] = {};
    double Psi2[nSub] = {};
    double deltaPsi2[nSub] = {};
    double phi, w;
    int PP, TT, EW, iSide = 0;


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

        TVector3 pVtx = event->primaryVertex();
        
        PrimeVertXvsY->Fill(pVtx.X(), pVtx.Y());
        PrimeVertZ->Fill(pVtx.Z());

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
	    Qvec[2*iSub] = 0.;
	    Qvec[2*iSub+1] = 0.;
	    QWeight[iSub] = 0.;
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

            DEdx->Fill(femtoTrack->pPtot()*femtoTrack->charge(), femtoTrack->dEdx());
            NsigmakaonMagn->Fill(femtoTrack->pPtot()*femtoTrack->charge(), femtoTrack->nSigmaKaon());
            NsigmapionMagn->Fill(femtoTrack->pPtot()*femtoTrack->charge(), femtoTrack->nSigmaPion());
            NsigmaprotMagn->Fill(femtoTrack->pPtot()*femtoTrack->charge(), femtoTrack->nSigmaProton());
            NsigmakaonTrans->Fill(femtoTrack->pPt()*femtoTrack->charge(), femtoTrack->nSigmaKaon());
            NsigmapionTrans->Fill(femtoTrack->pPt()*femtoTrack->charge(), femtoTrack->nSigmaPion());
            NsigmaprotTrans->Fill(femtoTrack->pPt()*femtoTrack->charge(), femtoTrack->nSigmaProton());
            NHits->Fill(femtoTrack->nHits());
            Pseudorap->Fill(eta);
            Mom->Fill(femtoTrack->pPtot());
            //BettaVsMom->Fill(femtoTrack->pPtot()*femtoTrack->charge(), femtoTrack->beta());

	    if (femtoTrack->nHits() < 15 || femtoTrack->pPt() < 0.1 || femtoTrack->pPt() > 2. || fabs(femtoTrack->pMom().Eta()) > 1.5 || fabs(femtoTrack->pMom().Eta())<0.1) continue;
            if ((Float_t)femtoTrack->nHits()/femtoTrack->nHitsPoss() < 0.52 ) continue;

	    /*if( !isGoodTrack( femtoTrack->pMom(), femtoTrack->nHits(),
                                femtoTrack->nHitsPoss() ) ) {
                continue;
            }*/

	    
            //if(TMath::Abs(eta)>1.5 || TMath::Abs(eta)<0.1) continue;
                                   
 
            double wEff = femtoTrack->pPt();
            double phi = femtoTrack->pMom().Phi(); //femtoTrack->isPrimary() ? femtoTrack->pMom().Phi() : 0.;
	    if(eta>0) iSide = 0;//West
	    else if(eta<0) iSide = 1;//East

	    Qvec[2*iSide] += wEff*TMath::Cos(2*phi);
            Qvec[2*iSide+1] += wEff*TMath::Sin(2*phi);
            QWeight[iSide] += wEff;

	    Qvec[4] += wEff*TMath::Cos(2*phi);
            Qvec[5] += wEff*TMath::Sin(2*phi);
            QWeight[2] += wEff;
	    QWeight[3] += wEff;

            
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

	    double wEff_EPD = w;
      	    if(wEff_EPD<0.3) continue;
      	    if(wEff_EPD>3) wEff_EPD = 3;

	    if(EW == 1){
        
        	Qvec[8] += wEff_EPD*TMath::Cos(2*phiCenter[PP-1][TT-1][EW]);
        	Qvec[9] += wEff_EPD*TMath::Sin(2*phiCenter[PP-1][TT-1][EW]);
        	QWeight[4] += wEff_EPD;
      	    }

      	    if(EW==-1){
        	EW = 0;
        	Qvec[10] += wEff_EPD*TMath::Cos(2*phiCenter[PP-1][TT-1][0]);
        	Qvec[11] += wEff_EPD*TMath::Sin(2*phiCenter[PP-1][TT-1][0]);
        	QWeight[5] += wEff_EPD;
      	    }
	}
	QWeight[6] = 1;

        
        //-------------------------------reaction plane--------------------------------------------------
        if(ScipZeroWeight(QWeight, nSub) == false) continue;
        
	//Get Q vectors
	int check = 0;            
        for(int iSub=0; iSub!=nSub; iSub++){
	    Qvec[2*iSub] = Qvec[2*iSub]/QWeight[iSub];
	    Qvec[2*iSub+1] = Qvec[2*iSub+1]/QWeight[iSub];
	    if(fabs(Qvec[2*iSub])>999 || fabs(Qvec[2*iSub+1])>999) check = 1;
	}
	if(check == 1) continue;
	Qvec[6] = Qvec[0] + Qvec[2];
	Qvec[7] = Qvec[1] + Qvec[3];
	Qvec[12] = Qvec[8] + Qvec[10];
	Qvec[13] = Qvec[9] + Qvec[11];

	
	//Centring collibration
	if(mode == "Centred" || mode == "Flatt"){
            //Recentering       
	    for(int iSub = 0; iSub!=2*nSub; iSub++){
	        Qvec[iSub] -= QvecProf_TH[iSub]->GetBinContent(cent+1);
	    }
        }

	
	//Get Psi2
	for(int iSub=0; iSub!=nSub; iSub++){
	    Psi2[iSub] = GetPsi(Qvec[2*iSub], Qvec[2*iSub+1]);
	    if(Psi2[iSub] > 3.5) check = 1;
	}
	if(check == 1) continue;
			            

	//Flattening collibration
        if(mode == "Flatt"){
            //Flattaning
            for(int iProf=0; iProf!=10; iProf++){
		for(int iSub=0; iSub!=nSub; iSub++){
		    deltaPsi2[iSub] += Coef_A_n_TH[iProf][iSub]->GetBinContent(cent+1)*TMath::Cos((iProf+1)*2*Psi2[iSub]) + \
                                       Coef_B_n_TH[iProf][iSub]->GetBinContent(cent+1)*TMath::Sin((iProf+1)*2*Psi2[iSub]);
		}	
            }
	    
	    for(int iSub=0; iSub!=nSub; iSub++){
		Psi2[iSub] += deltaPsi2[iSub]/2.;
		while(Psi2[iSub]>TMath::Pi()) Psi2[iSub] -= TMath::Pi();
                while(Psi2[iSub]<0.0) Psi2[iSub] += TMath::Pi();
                deltaPsi2[iSub] = 0.0;
	    }
        }    
                     
                                
        //Fill histograms
	for(int iSub = 0; iSub!=nSub; iSub++){
	    QvecHist[cent][2*iSub]->Fill(Qvec[2*iSub]);
	    QvecHist[cent][2*iSub+1]->Fill(Qvec[2*iSub+1]);
	    Psi2Hist[cent][iSub]->Fill(Psi2[iSub]);
	}            


        if(mode == "6sem"){
            //Profiles for recentering
            for(int iSub=0; iSub!=nSub*2; iSub++){
	        QvecProf[iSub]->Fill(cent, Qvec[iSub]);
	    }
        }    


        if(mode == "Centred"){
            //Profiles for flattening
            for(int iProf=0; iProf!=10; iProf++){
	        for(int iSub=0; iSub!=nSub; iSub++){
	      	    Coef_A_n[iProf][iSub]->Fill(cent, (-2.0/(iProf+1))*TMath::Sin((iProf+1)*2*Psi2[iSub]));
		    Coef_B_n[iProf][iSub]->Fill(cent, (2.0/(iProf+1))*TMath::Cos((iProf+1)*2*Psi2[iSub]));
		}
            }        
        }

	Psi2Full_EW[cent][0]->Fill(Psi2[2]-Psi2[0]);
	Psi2Full_EW[cent][1]->Fill(Psi2[2]-Psi2[1]);
	Psi2Full_EW[cent][2]->Fill(Psi2[3]-Psi2[0]);
	Psi2Full_EW[cent][3]->Fill(Psi2[3]-Psi2[1]);
	Psi2Full_EW[cent][4]->Fill(Psi2[2]-Psi2[3]);
	Psi2Full_EW[cent][5]->Fill(Psi2[0]-Psi2[1]);
        CosOfDiff->Fill(cent, TMath::Cos(2*(Psi2[1] - Psi2[0])));
	CosOfDiff_EPD->Fill(cent, TMath::Cos(2*(Psi2[2] - Psi2[4])));
	CosOfDiff_Comb->Fill(cent, TMath::Cos(2*(Psi2[2] - Psi2[5])));
        

    }//for(Long64_t iEvent=0; iEvent<events2read; iEvent++)

    //General distribution
    DEdx->GetXaxis()->SetTitle("p [GeV/c]");
    DEdx->GetYaxis()->SetTitle("dE/dx [GeV/m]");
    DEdx->GetXaxis()->SetLabelSize(0.05);
    DEdx->GetXaxis()->SetTitleSize(0.05);
    DEdx->GetYaxis()->SetLabelSize(0.05);
    DEdx->GetYaxis()->SetTitleSize(0.05);
       
    PrimeVertXvsY->GetXaxis()->SetTitle("x [sm]");
    PrimeVertXvsY->GetYaxis()->SetTitle("y [sm]");
    PrimeVertXvsY->GetXaxis()->SetLabelSize(0.05);
    PrimeVertXvsY->GetXaxis()->SetTitleSize(0.05);
    PrimeVertXvsY->GetYaxis()->SetLabelSize(0.05);
    PrimeVertXvsY->GetYaxis()->SetTitleSize(0.05);
    PrimeVertZ->GetXaxis()->SetTitle("z [sm]");
    PrimeVertZ->GetXaxis()->SetLabelSize(0.05);
    PrimeVertZ->GetXaxis()->SetTitleSize(0.05);

    //EP Distributions
    for(int iCent=0; iCent!=9; iCent++){
	QvecHist[iCent][0]->SetTitle(Form("Qx West TPC (%i)", iCent));
        QvecHist[iCent][1]->SetTitle(Form("Qy West TPC (%i)", iCent));
        QvecHist[iCent][2]->SetTitle(Form("Qx East TPC (%i)", iCent));
        QvecHist[iCent][3]->SetTitle(Form("Qy East TPC (%i)", iCent));
        QvecHist[iCent][4]->SetTitle(Form("Qx Comb TPC (%i)", iCent));
        QvecHist[iCent][5]->SetTitle(Form("Qy Comb TPC (%i)", iCent));
        QvecHist[iCent][6]->SetTitle(Form("Qx E+W TPC (%i)", iCent));
        QvecHist[iCent][7]->SetTitle(Form("Qy E+W TPC (%i)", iCent));
        QvecHist[iCent][8]->SetTitle(Form("Qx West EPD (%i)", iCent));
        QvecHist[iCent][9]->SetTitle(Form("Qy West EPD (%i)", iCent));
        QvecHist[iCent][10]->SetTitle(Form("Qx East EPD (%i)", iCent));
        QvecHist[iCent][11]->SetTitle(Form("Qy East EPD (%i)", iCent));
        QvecHist[iCent][12]->SetTitle(Form("Qx E+W EPD (%i)", iCent));
        QvecHist[iCent][13]->SetTitle(Form("Qy E+W EPD (%i)", iCent));

        Psi2Hist[iCent][0]->SetTitle(Form("West TPC (%i)", iCent));
        Psi2Hist[iCent][1]->SetTitle(Form("East TPC (%i)", iCent));
        Psi2Hist[iCent][2]->SetTitle(Form("Comb TPC (%i)", iCent));
        Psi2Hist[iCent][3]->SetTitle(Form("E+W TPC (%i)", iCent));
        Psi2Hist[iCent][4]->SetTitle(Form("West EPD (%i)", iCent));
        Psi2Hist[iCent][5]->SetTitle(Form("East EPD (%i)", iCent));
        Psi2Hist[iCent][6]->SetTitle(Form("E+W EPD (%i)", iCent));

	//Psi_i - Psi_j
	Psi2Full_EW[iCent][0]->SetTitle(Form("Comb - West (%i)", iCent));
	Psi2Full_EW[iCent][1]->SetTitle(Form("Comb - East (%i)", iCent));
	Psi2Full_EW[iCent][2]->SetTitle(Form("E+W - West (%i)", iCent));
	Psi2Full_EW[iCent][3]->SetTitle(Form("E+W - East (%i)", iCent));
	Psi2Full_EW[iCent][4]->SetTitle(Form("Comb - E+W (%i)", iCent));
	Psi2Full_EW[iCent][5]->SetTitle(Form("West - East (%i)", iCent));

    }
    QvecProf[0]->SetTitle("Profile Qx West TPC");
    QvecProf[1]->SetTitle("Profile Qy West TPC");
    QvecProf[2]->SetTitle("Profile Qx East TPC");
    QvecProf[3]->SetTitle("Profile Qy East TPC");
    QvecProf[4]->SetTitle("Profile Qx Comb TPC");
    QvecProf[5]->SetTitle("Profile Qy Comb TPC");
    QvecProf[6]->SetTitle("Profile Qx E+W TPC");
    QvecProf[7]->SetTitle("Profile Qy E+W TPC");
    QvecProf[8]->SetTitle("Profile Qx West EPD");
    QvecProf[9]->SetTitle("Profile Qy West EPD");
    QvecProf[10]->SetTitle("Profile Qx East EPD");
    QvecProf[11]->SetTitle("Profile Qy East EPD");
    QvecProf[12]->SetTitle("Profile Qx E+W EPD");
    QvecProf[13]->SetTitle("Profile Qy E+W EPD");

    for(int i=0; i!=10; i++){
	Coef_A_n[i][0]->SetTitle(Form("A_n_%i West TPC", i));
	Coef_B_n[i][0]->SetTitle(Form("B_n_%i West TPC", i));
	Coef_A_n[i][1]->SetTitle(Form("A_n_%i East TPC", i));
        Coef_B_n[i][1]->SetTitle(Form("B_n_%i East TPC", i));
	Coef_A_n[i][2]->SetTitle(Form("A_n_%i Comb TPC", i));
        Coef_B_n[i][2]->SetTitle(Form("B_n_%i Comb TPC", i));
	Coef_A_n[i][3]->SetTitle(Form("A_n_%i E+W TPC", i));
        Coef_B_n[i][3]->SetTitle(Form("B_n_%i E+W TPC", i));
	Coef_A_n[i][4]->SetTitle(Form("A_n_%i West EPD", i));
        Coef_B_n[i][4]->SetTitle(Form("B_n_%i West EPD", i));
	Coef_A_n[i][5]->SetTitle(Form("A_n_%i East EPD", i));
        Coef_B_n[i][5]->SetTitle(Form("B_n_%i East EPD", i));
	Coef_A_n[i][6]->SetTitle(Form("A_n_%i E+W EPD", i));
        Coef_B_n[i][6]->SetTitle(Form("B_n_%i E+W EPD", i));
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
double GetPsi(double Qx, double Qy){
    double Psi;
    Psi = atan2(Qy, Qx)/2.;
    if(Psi<0.) Psi += TMath::Pi();
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
