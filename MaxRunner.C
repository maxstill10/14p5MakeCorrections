/**
 *  Author: Grigory Nigmatkulov
 *  Date:   July 5, 2018
 *
 *  Description:
 *  This macros takes inFileName argument with a picoDst.root file
 *  or with a list of files (name.lis or name.list). It sets _VANILLA_ROOT_
 *  (necessary for standalone and can be skipped on RACF), loads pre compiled
 *  libStPicoDst.so (from StPicoEvent), compiles and executes a text
 *  PicoDstAnalyzer.C macro with passing inFileName to it, and
 *  cleans up the directory from the compilation products at the end.
 *
 *  Some details:
 *    inFileName - is a name of name.picoDst.root file or a name
 *                 of a name.lis(t) files that contains a list of
 *                 name1.picoDst.root files.
 *    NOTE: inFileName should contain either /absolutePath/inFileName
 *          or /relative2currentDir/inFileName
 *  It is assumed that PicoDstAnalyzer.C is placed in the same
 *  directory where the RunAnalyzer.C is stored.
 **/

#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"

//_________________
void MaxRunner(const Char_t *inFileName = "root://xrdstar.rcf.bnl.gov:1095//home/starlib/home/starreco/reco/production_14p5GeV_2019/ReversedFullField/P21ic/2019/140/20140001/st_physics_20140001_raw_0000002.picoDst.root", 
                 const Char_t *outFileName = "/star/u/mmorozov/MaximTest/out.root",
		 const Char_t *mode = "6sem") {
  // Next line is not needed if you are not running in a standalone mode
  cout << "Started Runner!" << endl;
  gSystem->Load("libStPicoDst.so");//libStPicoDst
  cout << "Loaded StPicoDst" << endl;
//  gSystem->Load("StEpdUtil");
  gSystem->Load("StRefMultCorr");
  cout << "Loaded StRefMultCorr" << endl;
//  gSystem->Load("StEpConstructor");
  TString str;
  str = ".x PicoDstAnalyzer.C+(\"";
  str += inFileName;
  str += "\",\"";
  str += outFileName;
  str += "\",\"";
  str += mode;
  str += "\")";
  gROOT->ProcessLine( str.Data() );
  // Next line should be commented if you run in a batch mode
  }
