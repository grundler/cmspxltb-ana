#include "tbAna.hh"
#include "constants.hh"

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

   for(int i=0; i<2; i++) {
      delete h_effFlux[i];
      delete h_effSpill[i];

      delete h_resSpill[i];
   }
}

void tbAna::plotEffVsIntensity(TCut myCut) {

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

         // TCanvas *c = new TCanvas("c","",800,600);
         // h->Draw();
         // c->Print(Form("res%s_spill%d.pdf",D[iD],iSpill));

         // cout << h->GetMean() << " " << h->GetRMS() << endl;
                  
         h_resSpill[iD][0]->Fill(iSpill, h->GetMean());
         h_resSpill[iD][1]->Fill(iSpill, h->GetRMS());
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

         //Track has passed, fill track-level information
         if(flux > -1) {
            h_effFlux[0]->Fill(flux);
         }
         h_effSpill[0]->Fill(iSpill);

         if(flux > maxFlux) 
            maxFlux = flux;

         //If hit on track, fill hit-level information
         if(goodHits->Contains(entry)) {
            // cout << "\t\t\tgoodHit\n";
            if(flux > -1) {
               h_effFlux[1]->Fill(flux);
            }
            h_effSpill[1]->Fill(iSpill);
         }

      }

      delete _trackTree;
      delete _tc;
   }

   cout << "max flux seen: " << maxFlux << endl;
   cout << "h_tracksFlux overflow: " << h_effFlux[0]->GetBinContent(h_effFlux[0]->GetNbinsX()+1) << endl;
   cout << "h_hitsFlux overflow: " << h_effFlux[1]->GetBinContent(h_effFlux[1]->GetNbinsX()+1) << endl;
   TGraphAsymmErrors *effVsFlux = new TGraphAsymmErrors();
   effVsFlux->BayesDivide(h_effFlux[1],h_effFlux[0]);
   TCanvas *c = new TCanvas("c1","", 800,600);
   effVsFlux->Draw();

   TGraphAsymmErrors *effVsSpill = new TGraphAsymmErrors();
   effVsSpill->BayesDivide(h_effSpill[1],h_effSpill[0]);
   c = new TCanvas("c2","", 800,600);
   effVsSpill->Draw();

   for(int iD=0; iD<nD; iD++) {
      for(int i=0; i<2; i++) {
         c = new TCanvas(Form("c%d",3+iD*2+i),"", 800, 600);
         h_resSpill[iD][i]->Draw();
      }
   }
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
   _trackTree->SetBranchStatus("nhits_7",kTRUE);

   _trackTree->SetBranchAddress("RunNr", &RunNr);
   _trackTree->SetBranchAddress("EvtNr", &EvtNr);
   _trackTree->SetBranchAddress("fitX_0", &fitX_0);
   _trackTree->SetBranchAddress("fitY_0", &fitY_0);
   _trackTree->SetBranchAddress("fitX_1", &fitX_1);
   _trackTree->SetBranchAddress("fitY_1", &fitY_1);
   _trackTree->SetBranchAddress("fitX_2", &fitX_2);
   _trackTree->SetBranchAddress("fitY_2", &fitY_2);
   _trackTree->SetBranchAddress("fitX_3", &fitX_3);
   _trackTree->SetBranchAddress("fitY_3", &fitY_3);
   _trackTree->SetBranchAddress("fitX_4", &fitX_4);
   _trackTree->SetBranchAddress("fitY_4", &fitY_4);
   _trackTree->SetBranchAddress("fitX_5", &fitX_5);
   _trackTree->SetBranchAddress("fitY_5", &fitY_5);
   _trackTree->SetBranchAddress("fitX_6", &fitX_6);
   _trackTree->SetBranchAddress("fitY_6", &fitY_6);
   _trackTree->SetBranchAddress("fitX_7", &fitX_7);
   _trackTree->SetBranchAddress("fitY_7", &fitY_7);
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
   for(int i=0; i<2; i++) { //loop over histograms (usually tracks and hits, for efficiency)
      if(0==i) sprintf(suffix, "tracks");
      else     sprintf(suffix, "hits");

      sprintf(name, "h_effFlux_%s",suffix);
      sprintf(title, "%s vs flux",suffix);
      //h_effFlux[i] = new TH1F(name, title, nIntBins, intHist);
      h_effFlux[i] = new TH1F(name, title, 100, 0., 3000.);

      sprintf(name, "h_effSpill_%s",suffix);
      sprintf(title, "%s vs spill",suffix);
      h_effSpill[i] = new TH1F(name, title, nSpills, _firstSpill-0.5, _finalSpill+0.5);

      if(0==i) sprintf(suffix, "mean");
      else     sprintf(suffix, "sigma");

      for(int iD=0; iD<nD; iD++) {
         sprintf(name, "h_res%sSpill_%s",D[iD],suffix);
         sprintf(title, "%s residual %s vs spill",D[iD],suffix);
         h_resSpill[iD][i] = new TH1F(name, title, nSpills, _firstSpill-0.5, _finalSpill+0.5);
      }
   }

}
