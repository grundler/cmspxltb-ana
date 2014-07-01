#include "tbAna.hh"
#include "treeCorrelator.hh"
#include "plotter.hh"

#include <cstring>
#include <iostream>
#include <sstream>

#include <sys/types.h>
#include <sys/stat.h>

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TEntryList.h>
#include <TGraphAsymmErrors.h>

using namespace std;

tbAna::tbAna(int dutID, string board, int spill0, int spill1, int algo) 
   : _testBoard(board),
     _DUTID(dutID),
     _firstSpill(spill0),
     _finalSpill(spill1),
     _nSpills(1+_finalSpill-_firstSpill),
     _algo(algo),
     maxFluxRatio(5.),
     use_correctTriggerPhase(true),
     use_isFiducial(true),
     use_slopeCut(true)
{
   //Set TCuts
   chi2_50 = "Chi2<50.";
   chi2_600 = "Chi2<600.";

   char tmptxt[256];
   sprintf(tmptxt,"((dut%s-fit%s_%d)<%f && (dut%s-fit%s_%d)>%f) && ((dut%s-fit%s_%d)<%f && (dut%s-fit%s_%d)>%f)",
           D[0], D[0], _DUTID, resHi[0], D[0], D[0], _DUTID, resLo[0],
           D[1], D[1], _DUTID, resHi[1], D[1], D[1], _DUTID, resLo[1]);
   nearTrack = tmptxt;

   ostringstream stream;
   stream << outDir << "/" << _testBoard << "_DUT" << _DUTID << "_" << _firstSpill << "-" << _finalSpill;
   _outDir = stream.str();
   mkdir(_outDir.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);

   _plotter = new plotter(_firstSpill,_finalSpill,_outDir);

   cout << "initialized tbAna\n";
}

tbAna::~tbAna() {
   delete _trackTree;
   
}

void tbAna::analyze(TCut myCut) {

   cout << "Running analysis of " << _testBoard << " spills " << _firstSpill << " to " << _finalSpill << endl;
   cout << "\tDUT: " << _DUTID << endl;
   cout << "\tQIE correlation algorithm: " << _algo << endl;
   cout << "\tRequire correct trigger phase: " << use_correctTriggerPhase << endl;
   cout << "\tRequire track to be fiducial: " << use_isFiducial << endl;
   cout << "\tRequire low slope: " << use_slopeCut << endl;
   cout << "\tMaximum flux ratio: " << maxFluxRatio << endl;
   cout << endl;

   float maxFlux = 0.;
   int maxSpill = 0, maxEvt = 0;

   for(int iSpill=_firstSpill;  iSpill<=_finalSpill; iSpill++) {
      if(!initSpill(iSpill)) {
         cout << "Initializing spill " << iSpill << " failed. Skipping.\n";
         continue;
      }

      // TCut aCut = chi2_50 && slopeCut && isFiducial && myCut;
      TCut aCut = chi2_50 && myCut;
      if(use_slopeCut)
         aCut += slopeCut;
      if(use_isFiducial)
         aCut += isFiducial;
      //aCut.Print();
      TCut bCut = aCut && nearTrack;
      //bCut.Print();
      _trackTree->Draw(">>goodtracks", aCut, "entrylist");
      TEntryList *goodTracks = (TEntryList*) gDirectory->Get("goodtracks");
      _trackTree->SetEntryList(goodTracks);
      _trackTree->Draw(">>goodhits", bCut, "entrylist");
      TEntryList *goodHits = (TEntryList*) gDirectory->Get("goodhits");

      for(int iD=0; iD<nD; iD++) {
         char hh[100];
         char var[100];

         sprintf(hh,"hres%s",D[iD]);
         sprintf(var,"dut%s-fit%s_%d >> %s", D[iD], D[iD], _DUTID, hh);         
         _trackTree->Draw(var,bCut, "goff");
         TH1F *h = (TH1F*) gDirectory->Get(hh);

         _plotter->h_resSpill[_wbcBin][iD][0]->SetBinContent(1+iSpill-_firstSpill, h->GetMean());
         _plotter->h_resSpill[_wbcBin][iD][0]->SetBinError(1+iSpill-_firstSpill, h->GetMeanError());
         _plotter->h_resSpill[_wbcBin][iD][1]->SetBinContent(1+iSpill-_firstSpill, h->GetRMS());
         _plotter->h_resSpill[_wbcBin][iD][1]->SetBinError(1+iSpill-_firstSpill, h->GetRMSError());
      }

      cout << "\tgoodTracks: " << goodTracks->GetN() << endl;
      //goodTracks->Print("all");
      cout << "\tgoodHits: " << goodHits->GetN() << endl;
      //goodHits->Print("all");

      for(int iEvt=0; iEvt<(int)goodTracks->GetN(); iEvt++) {
         if(iEvt%5000==0)
            cout << "\t\tProcessing event " << iEvt << endl;

         int entry = goodTracks->GetEntry(iEvt);
         loadTrackEntry(entry);

         int TriggerPhase = _tc->getTriggerPhase(EvtNr);
         if(use_correctTriggerPhase && TriggerPhase != correctTriggerPhase)
            continue;

         int qieevt = _tc->getQieEvent(EvtNr);
         float flux = _tc->getFlux(qieevt);
         float ratio = _tc->getFluxRatio(qieevt);

         // cout << "\tiEvt " << iEvt << " qieevt " << qieevt << " flux " << flux << endl;
         if(flux < 0.) continue;
         if(ratio > maxFluxRatio) continue;

         //Track has passed, fill track-level information
         _plotter->h_effFlux[_wbcBin][0]->Fill(flux);
         _plotter->h_effNHits[_wbcBin][0]->Fill(nhits_4);
         _plotter->h_effSpill[_wbcBin][0]->Fill(iSpill);
         _plotter->h_nHitsFlux[_wbcBin]->Fill(flux,nhits_4);

         float dutFitX=-999., dutFitY=-999.;
         getDUTFitPosition(dutFitX, dutFitY);

         _plotter->h_effMap[_wbcBin][0]->Fill(dutFitX, dutFitY);
         _plotter->h_effMapWide[_wbcBin][0]->Fill(dutFitX, dutFitY);

         if(flux > maxFlux) { 
            maxFlux = flux;
            maxSpill = iSpill;
            maxEvt = EvtNr;
         }

         //If hit on track, fill hit-level information
         if(goodHits->Contains(entry)) {
            // cout << "\t\t\tgoodHit\n";
            _plotter->h_effFlux[_wbcBin][1]->Fill(flux);
            _plotter->h_effNHits[_wbcBin][1]->Fill(nhits_4);
            _plotter->h_effSpill[_wbcBin][1]->Fill(iSpill);

            _plotter->h_effMap[_wbcBin][1]->Fill(dutFitX, dutFitY);
            _plotter->h_effMapWide[_wbcBin][1]->Fill(dutFitX, dutFitY);
         }

      }

      delete _trackTree;
      delete _tc;
   }

   cout << "max flux seen: " << maxFlux << ", at spill/event " << maxSpill << "/" << maxEvt << endl;

   ostringstream stream;
   stream << _outDir << "/" << _testBoard << "_DUT" << _DUTID << "_" << _firstSpill << "-" << _finalSpill << ".root";
   string filename = stream.str();
   _plotter->writeFile(filename);
}

void tbAna::loadTrackEntry(int entry) {

   _trackTree->LoadTree(entry);
   _trackTree->GetEntry(entry);
 
}

void tbAna::getDUTFitPosition(float &posX, float &posY) {
   if(_DUTID == 0) {
      posX = fitX_0;      posY = fitY_0;
   }
   else if(_DUTID == 1) {
      posX = fitX_1;      posY = fitY_1;
   }
   else if(_DUTID == 2) {
      posX = fitX_2;      posY = fitY_2;
   }
   else if(_DUTID == 3) {
      posX = fitX_3;      posY = fitY_3;
   }
   else if(_DUTID == 4) {
      posX = fitX_4;      posY = fitY_4;
   }
   else if(_DUTID == 5) {
      posX = fitX_5;      posY = fitY_5;
   }
   else if(_DUTID == 6) {
      posX = fitX_6;      posY = fitY_6;
   }
   else if(_DUTID == 7) {
      posX = fitX_7;      posY = fitY_7;
   }

}

bool tbAna::initSpill(int spill) {
   cout << "get spill " << spill << endl;

   //Get Tree
   initTrackTree(spill);
   if(_trackTree==NULL) {
      cout << "\tCould not get all information for spill " << spill << ". Skipping\n";
      return false;
   }

   //Set cuts
   if(use_slopeCut && !setSlopeCut())
      return false;

   if(use_isFiducial && !setFiducialCut())
      return false;

   //Get maps to correlate with triggerphase and qie trees
   _tc = new treeCorrelator(spill,_testBoard, _algo);
   if(!_tc->isInitialized())
      return false;

   //get WBC;
   _wbcBin=-1;
   int wbc=_tc->getWBC();
   cout << "  WBC=" << wbc << endl;
   for(int iwbc=0; iwbc<nWBC; iwbc++) {
      if(wbc == WBCvalue[iwbc]) {
         _wbcBin=iwbc;
         break;
      }
   }
   if(_wbcBin < 0) {
      cout << "\tCould not find proper WBC value for spill " << spill << ". Found WBC=" << wbc << ". Skipping\n";
      return false;
   }

   //Set one more cut
   if(use_correctTriggerPhase && !setTriggerPhaseCut(spill))
      return false;

   return true;
}

bool tbAna::setSlopeCut() {

   slopeCut = "";
   double slopecuts[nD][2] = {{-10.,10.},{-10.,10.}}; //starting min/max cuts for x,y

   char hh[100];
   char var[100];

   TH1F *h[nD];

   float hMax[nD];
   bool foundMin[nD];

   for(int iD=0; iD<nD; iD++) {
      sprintf(hh,"hslope%s",D[iD]);
      sprintf(var, "fit%s_0-fit%s_7 >> %s",D[iD],D[iD],hh);

      _trackTree->Draw(var,chi2_50,"goff");
      h[iD] = (TH1F*) gDirectory->Get(hh);

      if(!h[iD]) {
         cout << "\tDon't have histogram. Tree empty?\n";
         return false;
      }
      hMax[iD] = h[iD]->GetMaximum();
      foundMin[iD] = false;

      for(int ibin=1; ibin <= h[iD]->GetNbinsX(); ibin++) {
         if(h[iD]->GetBinContent(ibin) > hMax[iD]*0.5) {
            if(!foundMin[iD]) {
               slopecuts[iD][0] = h[iD]->GetBinCenter(ibin-2);
               foundMin[iD] = true;
            }
            slopecuts[iD][1] = h[iD]->GetBinCenter(ibin+2);
         }
      }
   }

   slopeCut = Form("( ((fitX_0-fitX_7)>%f && (fitX_0-fitX_7)<%f) && ((fitY_0-fitY_7)>%f && (fitY_0-fitY_7)<%f) )",slopecuts[0][0],slopecuts[0][1],slopecuts[1][0],slopecuts[1][1]);
   cout << "\tSet slope cut to " << slopeCut << endl;

   return true;
}

bool tbAna::setFiducialCut() {
   isFiducial = "";

   TH1F *h[nD];
   double fidcuts[nD][2] = {{0.,0.},{0.,0.}}; //starting min/max cuts for x,y

   char hh[100];
   char var[100];
   for(int iD=0; iD<nD; iD++) {
      sprintf(hh,"hfit%s",D[iD]);
      sprintf(var, "fit%s_%d >> %s",D[iD],_DUTID,hh);

      _trackTree->Draw(var,chi2_50, "goff");
      h[iD] = (TH1F*) gDirectory->Get(hh);

      if(!h[iD]) {
         cout << "\tDon't have histogram. Tree empty?\n";
         return false;
      }

      double tmpmin=0., tmpmax=0.;
      for(int ibin=1; ibin <= h[iD]->GetNbinsX(); ibin++) {
         if(h[iD]->GetBinContent(ibin) > 0.) {
            tmpmin = h[iD]->GetBinCenter(ibin);
            break;
         }
      }
      for(int ibin=h[iD]->GetNbinsX(); ibin>=1; ibin--) {
         if(h[iD]->GetBinContent(ibin) > 0.) {
            tmpmax = h[iD]->GetBinCenter(ibin);
            break;
         }
      }

      fidcuts[iD][0] = tmpmin + (tmpmax-tmpmin)*0.05;
      fidcuts[iD][1] = tmpmax - (tmpmax-tmpmin)*0.05;
   }

   isFiducial = Form("((fitX_%d>%f && fitX_%d<%f) && (fitY_%d>%f && fitY_%d<%f))",_DUTID,fidcuts[0][0],_DUTID,fidcuts[0][1],_DUTID,fidcuts[1][0],_DUTID,fidcuts[1][1]);
   cout << "\tSet fiducial cut to " << isFiducial << endl;


   return true;
}

bool tbAna::setTriggerPhaseCut(int spill) {

   correctTriggerPhase = -1;

   TCut aCut = chi2_50;
   TCut bCut = aCut && nearTrack;
   _trackTree->Draw(">>tptracks", aCut, "entrylist");
   TEntryList *tpTracks = (TEntryList*) gDirectory->Get("tptracks");
   _trackTree->SetEntryList(tpTracks);
   _trackTree->Draw(">>tphits", bCut, "entrylist");
   TEntryList *tpHits = (TEntryList*) gDirectory->Get("tphits");

   for(int iEvt=0; iEvt<(int)tpTracks->GetN(); iEvt++) {
      int entry = tpTracks->GetEntry(iEvt);
      loadTrackEntry(entry);      

      int TriggerPhase = _tc->getTriggerPhase(EvtNr);
      // cout << "iEvt " << iEvt << ", event " << EvtNr << ", tp " << TriggerPhase;
      if(tpHits->Contains(entry)) {
         _plotter->h_tpSpill[_wbcBin][1]->Fill(spill, TriggerPhase);
         // cout << ", has hit";
      }
      // cout << endl;

      _plotter->h_tpSpill[_wbcBin][0]->Fill(spill, TriggerPhase);
   }

   TGraphAsymmErrors *tgaeTP = new TGraphAsymmErrors();
   tgaeTP->SetName("tgaeTP");
   tgaeTP->SetTitle("Efficiency vs TriggerPhase");

   int spillBin = 1 + spill - _firstSpill;
   tgaeTP->Divide(_plotter->h_tpSpill[_wbcBin][1]->ProjectionY("_py",spillBin,spillBin,"e"),
                  _plotter->h_tpSpill[_wbcBin][0]->ProjectionY("_py",spillBin,spillBin,"e"),
                  "cl=0.683 b(1,1) mode");

   int tpMax=-1;
   float effMax=0.;
   for(int ibin=0; ibin<tgaeTP->GetN(); ibin++) {
      double tp,eff;
      tgaeTP->GetPoint(ibin, tp, eff);
      cout << " TriggerPhase " << ibin << ", tp,eff=" << tp << "," << eff << endl; 
      if(eff>effMax) {
         tpMax = ibin;
         effMax = eff;
      }
   }

   correctTriggerPhase = tpMax;
   cout << "\tCorrect TriggerPhase is " << correctTriggerPhase << endl;

   if(correctTriggerPhase == -1) {
      cout << "ERROR: Somehow we didn't get a trigger phase\n";
      return false;
   }
   return true;
}

void tbAna::initTrackTree(int spill) {

   //delete _trackTree;

   ostringstream stream;
   stream << subdir << "/" << _testBoard << "/histograms/" << spill << "-tracks.root";
   string filename = stream.str();

   TFile *f = new TFile(filename.c_str());
   if(f->IsZombie()) {
      cout << "\tFile " << filename << " does not exist\n";
      _trackTree = NULL;
      return;
   }

   _trackTree = new TChain("MyEUTelFitTuple/EUFit");

   _trackTree->SetBranchStatus("*",kFALSE);

   _trackTree->SetBranchStatus("RunNr",kTRUE);
   _trackTree->SetBranchStatus("EvtNr",kTRUE);
   _trackTree->SetBranchStatus("fitX_*",kTRUE);
   _trackTree->SetBranchStatus("fitY_*",kTRUE);
   _trackTree->SetBranchStatus("dutX",kTRUE);
   _trackTree->SetBranchStatus("dutY",kTRUE);
   _trackTree->SetBranchStatus("nTrack",kTRUE);
   _trackTree->SetBranchStatus("Ndf",kTRUE);
   _trackTree->SetBranchStatus("Chi2",kTRUE);
   _trackTree->SetBranchStatus("nhits_*",kTRUE);

   _trackTree->SetBranchAddress("RunNr", &RunNr);
   _trackTree->SetBranchAddress("EvtNr", &EvtNr);
   _trackTree->SetBranchAddress("fitX_0", &fitX_0);
   _trackTree->SetBranchAddress("fitY_0", &fitY_0);
   _trackTree->SetBranchAddress("nhits_0", &nhits_0);
   _trackTree->SetBranchAddress("fitX_1", &fitX_1);
   _trackTree->SetBranchAddress("fitY_1", &fitY_1);
   _trackTree->SetBranchAddress("nhits_1", &nhits_1);
   _trackTree->SetBranchAddress("fitX_2", &fitX_2);
   _trackTree->SetBranchAddress("fitY_2", &fitY_2);
   _trackTree->SetBranchAddress("nhits_2", &nhits_2);
   _trackTree->SetBranchAddress("fitX_3", &fitX_3);
   _trackTree->SetBranchAddress("fitY_3", &fitY_3);
   _trackTree->SetBranchAddress("nhits_3", &nhits_3);
   _trackTree->SetBranchAddress("fitX_4", &fitX_4);
   _trackTree->SetBranchAddress("fitY_4", &fitY_4);
   _trackTree->SetBranchAddress("nhits_4", &nhits_4);
   _trackTree->SetBranchAddress("fitX_5", &fitX_5);
   _trackTree->SetBranchAddress("fitY_5", &fitY_5);
   _trackTree->SetBranchAddress("nhits_5", &nhits_5);
   _trackTree->SetBranchAddress("fitX_6", &fitX_6);
   _trackTree->SetBranchAddress("fitY_6", &fitY_6);
   _trackTree->SetBranchAddress("nhits_6", &nhits_6);
   _trackTree->SetBranchAddress("fitX_7", &fitX_7);
   _trackTree->SetBranchAddress("fitY_7", &fitY_7);
   _trackTree->SetBranchAddress("nhits_7", &nhits_7);
   _trackTree->SetBranchAddress("dutX", &dutX);
   _trackTree->SetBranchAddress("dutY", &dutY);
   _trackTree->SetBranchAddress("nTrack", &nTrack);

   _trackTree->SetBranchAddress("Ndf", &Ndf);
   _trackTree->SetBranchAddress("Chi2",&Chi2);
   _trackTree->SetBranchAddress("nhits_7",&nhits_7);

   _trackTree->Add(filename.c_str());

}

