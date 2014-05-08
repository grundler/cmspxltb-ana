#include "tbAna.hh"
#include "treeCorrelator.hh"

#include <cstring>
#include <iostream>
#include <sstream>

#include <TFile.h>
#include <TH1F.h>
#include <TLine.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TEntryList.h>
#include <TGraphAsymmErrors.h>

using namespace std;

tbAna::tbAna(int dutID, string board, int spill0, int spill1) 
   : _testBoard(board),
     _DUTID(dutID),
     _firstSpill(spill0),
     _finalSpill(spill1)
{
   //Set TCuts
   chi2_50 = "Chi2<50.";
   chi2_600 = "Chi2<600.";

   char tmptxt[256];
   sprintf(tmptxt,"((dut%s-fit%s_%d)<%f && (dut%s-fit%s_%d)>%f) && ((dut%s-fit%s_%d)<%f && (dut%s-fit%s_%d)>%f)",
           D[0], D[0], _DUTID, resHi[0], D[0], D[0], _DUTID, resLo[0],
           D[1], D[1], _DUTID, resHi[1], D[1], D[1], _DUTID, resLo[1]);
   nearTrack = tmptxt;

   //Book histograms
   bookHistos();

   cout << "initialized tbAna\n";
}

tbAna::~tbAna() {
   delete _trackTree;

   for(int iwbc=0; iwbc<nWBC; iwbc++)
      for(int i=0; i<2; i++) {
         delete h_effFlux[iwbc][i];
         delete h_effNHits[iwbc][i];
         delete h_effSpill[iwbc][i];
         
         delete h_resSpill[iwbc][i];
      }
}

void tbAna::analyze(TCut myCut) {

   int maxFlux = 0.;

   for(int iSpill=_firstSpill;  iSpill<=_finalSpill; iSpill++) {
      if(!initSpill(iSpill)) {
         cout << "Initializing spill " << iSpill << " failed. Skipping.\n";
         continue;
      }

      TCut aCut = chi2_50 && slopeCut && isFiducial && myCut;
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

         h_resSpill[_wbcBin][iD][0]->SetBinContent(1+iSpill-_firstSpill, h->GetMean());
         h_resSpill[_wbcBin][iD][0]->SetBinError(1+iSpill-_firstSpill, h->GetMeanError());
         h_resSpill[_wbcBin][iD][1]->SetBinContent(1+iSpill-_firstSpill, h->GetRMS());
         h_resSpill[_wbcBin][iD][1]->SetBinError(1+iSpill-_firstSpill, h->GetRMSError());
      }

      cout << "\tgoodTracks: " << goodTracks->GetN() << endl;
      //goodTracks->Print("all");
      cout << "\tgoodHits: " << goodHits->GetN() << endl;
      //goodHits->Print("all");

      for(int iEvt=0; iEvt<(int)goodTracks->GetN(); iEvt++) {
         if(iEvt%1000==0)
            cout << "\t\tProcessing event " << iEvt << endl;

         int entry = goodTracks->GetEntry(iEvt);
         loadTrackEntry(entry);

         int TriggerPhase = _tc->getTriggerPhase(EvtNr);
         if(TriggerPhase != correctTriggerPhase)
            continue;

         int qieevt = _tc->getQieEvent(EvtNr);
         int flux = _tc->getFlux(qieevt);

         // cout << "\tiEvt " << iEvt << " qieevt " << qieevt << " flux " << flux << endl;
         if(flux < 0.)
            continue;

         //Track has passed, fill track-level information
         h_effFlux[_wbcBin][0]->Fill(flux);
         h_effNHits[_wbcBin][0]->Fill(nhits_4);
         h_effSpill[_wbcBin][0]->Fill(iSpill);

         if(flux > maxFlux) 
            maxFlux = flux;

         //If hit on track, fill hit-level information
         if(goodHits->Contains(entry)) {
            // cout << "\t\t\tgoodHit\n";
            h_effFlux[_wbcBin][1]->Fill(flux);
            h_effNHits[_wbcBin][1]->Fill(nhits_4);
            h_effSpill[_wbcBin][1]->Fill(iSpill);
         }

      }

      delete _trackTree;
      delete _tc;
   }

   cout << "max flux seen: " << maxFlux << endl;

   ostringstream stream;
   stream << "output_" << _testBoard << "_DUT" << _DUTID << "_" << _firstSpill << "-" << _finalSpill << ".root";
   string filename = stream.str();

   TFile *f = new TFile(filename.c_str(),"recreate");

   for(int iwbc=0; iwbc<nWBC; iwbc++) {
      g_effFlux[iwbc]->BayesDivide(h_effFlux[iwbc][1],h_effFlux[iwbc][0]);
      g_effNHits[iwbc]->BayesDivide(h_effNHits[iwbc][1],h_effNHits[iwbc][0]);
      g_effSpill[iwbc]->BayesDivide(h_effSpill[iwbc][1],h_effSpill[iwbc][0]);

      g_effFlux[iwbc]->Write();
      g_effNHits[iwbc]->Write();
      g_effSpill[iwbc]->Write();

      for(int i=0; i<2; i++) {
         h_effFlux[iwbc][i]->Write();
         h_effNHits[iwbc][i]->Write();
         h_effSpill[iwbc][i]->Write();
         for(int iD=0; iD<nD; iD++) {
            h_resSpill[iwbc][iD][i]->Write();
         }
      }
   }
   f->Close();
}


void tbAna::loadTrackEntry(int entry) {

   int lflag = _trackTree->LoadTree(entry);
   _trackTree->GetBranch("EvtNr")->GetEntry(lflag);
 
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
   if(!setSlopeCut())
      return false;

   if(!setFiducialCut())
      return false;

   //Get maps to correlate with triggerphase and qie trees
   _tc = new treeCorrelator(spill,_testBoard);
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
   if(!setTriggerPhaseCut())
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

bool tbAna::setTriggerPhaseCut() {

   correctTriggerPhase = -1;

   TCut aCut = chi2_50;
   TCut bCut = aCut && nearTrack;
   _trackTree->Draw(">>tptracks", aCut, "entrylist");
   TEntryList *tpTracks = (TEntryList*) gDirectory->Get("tptracks");
   _trackTree->SetEntryList(tpTracks);
   _trackTree->Draw(">>tphits", bCut, "entrylist");
   TEntryList *tpHits = (TEntryList*) gDirectory->Get("tphits");

   int nPhases = 8;
   char title[128], name[256];
   sprintf(name, "hitsTP");
   sprintf(title, "HitsOnTrack vs TriggerPhase");
   TH1F* hitsTP = new TH1F(name, title, nPhases, -0.5, nPhases-0.5);

   sprintf(name, "tracksTP");
   sprintf(title, "Tracks vs TriggerPhase");
   TH1F* tracksTP = new TH1F(name, title, nPhases, -0.5, nPhases-0.5);

   for(int iEvt=0; iEvt<(int)tpTracks->GetN(); iEvt++) {
      int entry = tpTracks->GetEntry(iEvt);
      loadTrackEntry(entry);      

      int TriggerPhase = _tc->getTriggerPhase(EvtNr);
      if(tpHits->Contains(entry)) {
         hitsTP->Fill(TriggerPhase);
      }

      tracksTP->Fill(TriggerPhase);
   }

   TGraphAsymmErrors *tgaeTP = new TGraphAsymmErrors();
   tgaeTP->SetName("tgaeTP");
   tgaeTP->SetTitle("Efficiency vs TriggerPhase");
   tgaeTP->BayesDivide(hitsTP,tracksTP);

   int tpMax=-1;
   float effMax=0.;
   for(int ibin=0; ibin<tgaeTP->GetN(); ibin++) {
      double tp,eff;
      tgaeTP->GetPoint(ibin, tp, eff);
      // cout << " TriggerPhase " << ibin << ", tp,eff=" << tp << "," << eff << endl; 
      if(eff>effMax) {
         tpMax = ibin;
         effMax = eff;
      }
   }

   correctTriggerPhase = tpMax;
   cout << "\tCorrect TriggerPhase is " << correctTriggerPhase << endl;

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

void tbAna::bookHistos() {

   char title[128], name[256];

   // const int nIntBins = 11;
   // const float intHist[nIntBins+1] = {0,25,50,100,200,400,800,1600,3200,6400,12800,20000};

   int nSpills = 1+_finalSpill-_firstSpill;

   char suffix[16];
   for(int iwbc=0; iwbc<nWBC; iwbc++) {
      for(int i=0; i<2; i++) { //loop over histograms (usually tracks and hits, for efficiency)
         if(0==i) sprintf(suffix, "tracks");
         else     sprintf(suffix, "hits");

         sprintf(name, "h_effFlux_wbc%d_%s",WBCvalue[iwbc],suffix);
         sprintf(title, "%s vs flux (WBC %d)",suffix, WBCvalue[iwbc]);
         //h_effFlux[iwbc][i] = new TH1F(name, title, nIntBins, intHist);
         h_effFlux[iwbc][i] = new TH1F(name, title, 200, 0., 1000.);

         sprintf(name, "h_effNHits_wbc%d_%s",WBCvalue[iwbc],suffix);
         sprintf(title, "%s vs pixel hits (WBC %d)",suffix,WBCvalue[iwbc]);
         h_effNHits[iwbc][i] = new TH1F(name, title, 100, -0.5, 99.5);

         sprintf(name, "h_effSpill_wbc%d_%s",WBCvalue[iwbc],suffix);
         sprintf(title, "%s vs spill (WBC %d)",suffix,WBCvalue[iwbc]);
         h_effSpill[iwbc][i] = new TH1F(name, title, nSpills, _firstSpill-0.5, _finalSpill+0.5);

         if(0==i) sprintf(suffix, "mean");
         else     sprintf(suffix, "sigma");

         for(int iD=0; iD<nD; iD++) {
            sprintf(name, "h_res%sSpill_wbc%d_%s",D[iD],WBCvalue[iwbc],suffix);
            sprintf(title, "%s residual %s vs spill (WBC %d)",D[iD],suffix,WBCvalue[iwbc]);
            h_resSpill[iwbc][iD][i] = new TH1F(name, title, nSpills, _firstSpill-0.5, _finalSpill+0.5);
         }
      }

      sprintf(name, "g_effFlux_wbc%d", WBCvalue[iwbc]);
      sprintf(title, "Efficiency vs flux (WBC %d)", WBCvalue[iwbc]);
      g_effFlux[iwbc] = new TGraphAsymmErrors();
      g_effFlux[iwbc]->SetName(name);
      g_effFlux[iwbc]->SetTitle(title);

      sprintf(name, "g_effNHits_wbc%d", WBCvalue[iwbc]);
      sprintf(title, "Efficiency vs pixel hits (WBC %d)", WBCvalue[iwbc]);
      g_effNHits[iwbc] = new TGraphAsymmErrors();
      g_effNHits[iwbc]->SetName(name);
      g_effNHits[iwbc]->SetTitle(title);

      sprintf(name, "g_effSpill_wbc%d", WBCvalue[iwbc]);
      sprintf(title, "Efficiency vs spill (WBC %d)", WBCvalue[iwbc]);
      g_effSpill[iwbc] = new TGraphAsymmErrors();
      g_effSpill[iwbc]->SetName(name);
      g_effSpill[iwbc]->SetTitle(title);

   }

}
