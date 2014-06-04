#include "tbAna.hh"
#include "treeCorrelator.hh"
#include "utils.hh"

#include <cstring>
#include <iostream>
#include <sstream>

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
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
     _algo(algo)
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

   for(int iwbc=0; iwbc<nWBC; iwbc++) {
      for(int i=0; i<3; i++) {
         delete h_effFlux[iwbc][i];
         delete h_effNHits[iwbc][i];
         delete h_effSpill[iwbc][i];

         delete h_effMap[iwbc][i];
         delete h_effMapWide[iwbc][i];
         
         if(i<2)
            for(int j=0; j<nD; j++)
               delete h_resSpill[iwbc][j][i];
      }
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
         if(iEvt%5000==0)
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

         //add ratio cut?

         //Track has passed, fill track-level information
         h_effFlux[_wbcBin][0]->Fill(flux);
         h_effNHits[_wbcBin][0]->Fill(nhits_4);
         h_effSpill[_wbcBin][0]->Fill(iSpill);

         float dutFitX=-999., dutFitY=-999.;
         getDUTFitPosition(dutFitX, dutFitY);

         h_effMap[_wbcBin][0]->Fill(dutFitX, dutFitY);
         h_effMapWide[_wbcBin][0]->Fill(dutFitX, dutFitY);

         if(flux > maxFlux) 
            maxFlux = flux;

         //If hit on track, fill hit-level information
         if(goodHits->Contains(entry)) {
            // cout << "\t\t\tgoodHit\n";
            h_effFlux[_wbcBin][1]->Fill(flux);
            h_effNHits[_wbcBin][1]->Fill(nhits_4);
            h_effSpill[_wbcBin][1]->Fill(iSpill);

            h_effMap[_wbcBin][1]->Fill(dutFitX, dutFitY);
            h_effMapWide[_wbcBin][1]->Fill(dutFitX, dutFitY);
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
      h_effMap[iwbc][2]->Divide(h_effMap[iwbc][1],h_effMap[iwbc][0],1,1,"B");
      h_effMapWide[iwbc][2]->Divide(h_effMapWide[iwbc][1],h_effMapWide[iwbc][0],1,1,"B");

      for(int i=0; i<3; i++) {
         h_effFlux[iwbc][i]->Write();
         h_effNHits[iwbc][i]->Write();
         h_effSpill[iwbc][i]->Write();
         h_effMap[iwbc][i]->Write();
         h_effMapWide[iwbc][i]->Write();
         if(i<2)
            for(int iD=0; iD<nD; iD++) {
               h_resSpill[iwbc][iD][i]->Write();
            }
      }
   }
   f->Close();
}

void tbAna::makePlots() {
   utils* util = new utils();

   TCanvas * slide;

   char slideT[256], slideN[256], slideF[256];

   TGraphAsymmErrors *g_effFlux[nWBC];
   TGraphAsymmErrors *g_effNHits[nWBC];
   TGraphAsymmErrors *g_effSpill[nWBC];
   for(int iwbc=0; iwbc<nWBC; iwbc++) {
      g_effFlux[iwbc] = new TGraphAsymmErrors();
      g_effFlux[iwbc]->Divide(h_effFlux[iwbc][1],h_effFlux[iwbc][0],"cl=0.683 b(1,1) mode");
      util->graphSetting(g_effFlux[iwbc],h_effFlux[iwbc][2],
                         WBCcolor[iwbc],WBCstyle[iwbc],
                         TString::Format("Flux"),TString::Format("Efficiency"));

      g_effNHits[iwbc] = new TGraphAsymmErrors();
      g_effNHits[iwbc]->Divide(h_effNHits[iwbc][1],h_effNHits[iwbc][0],"cl=0.683 b(1,1) mode");
      util->graphSetting(g_effNHits[iwbc],h_effNHits[iwbc][2],
                         WBCcolor[iwbc],WBCstyle[iwbc],
                         TString::Format("NHits"),TString::Format("Efficiency"));
      
      g_effSpill[iwbc] = new TGraphAsymmErrors();
      g_effSpill[iwbc]->Divide(h_effSpill[iwbc][1],h_effSpill[iwbc][0],"cl=0.683 b(1,1) mode");
      util->graphSetting(g_effSpill[iwbc],h_effSpill[iwbc][2],
                         WBCcolor[iwbc],WBCstyle[iwbc],
                         TString::Format("Spill"),TString::Format("Efficiency"));
   }

   slide = util->newSlide(TString::Format("eff_vs_flux_%s",_testBoard.c_str()),"");
   TLegend* leg = new TLegend(0.1,0.15,0.25,0.35);
   util->legendSetting(leg);
   h_effFlux[0][2]->GetXaxis()->SetRangeUser(0.,500.);
   h_effFlux[0][2]->Draw();
   for(int iwbc=0; iwbc<nWBC; iwbc++) {
      if(h_effFlux[iwbc][0]->GetEntries() == 0) continue;
      g_effFlux[iwbc]->Draw("pe same");
      leg->AddEntry(g_effFlux[iwbc],Form("WBC=%d",WBCvalue[iwbc]),"p");
   }
   leg->Draw();
   slide->SaveAs(Form("%s/eff_vs_flux_%s.%s",plotDir,_testBoard.c_str(),plotExt));

   slide = util->newSlide(TString::Format("eff_vs_nhits_%s",_testBoard.c_str()),"");
   h_effNHits[0][2]->GetXaxis()->SetRange(1,25);
   h_effNHits[0][2]->Draw();
   for(int iwbc=0; iwbc<nWBC; iwbc++) {
      if(h_effNHits[iwbc][0]->GetEntries() == 0) continue;
      g_effNHits[iwbc]->Draw("pe same");
   }
   leg->Draw();
   slide->SaveAs(Form("%s/eff_vs_nhits_%s.%s",plotDir,_testBoard.c_str(),plotExt));

   int nSpills = 1+_finalSpill-_firstSpill;

   slide = util->newSlide(TString::Format("eff_vs_spill_%s",_testBoard.c_str()),"");
   h_effSpill[0][2]->GetXaxis()->SetNoExponent();
   h_effSpill[0][2]->Draw();
   for(int iwbc=0; iwbc<nWBC; iwbc++) {
      if(h_effSpill[iwbc][0]->GetEntries() == 0) continue;
      g_effSpill[iwbc]->Draw("pe same");
   }
   leg->Draw();
   slide->SaveAs(Form("%s/eff_vs_spill_%s.%s",plotDir,_testBoard.c_str(),plotExt));

   for(int iD=0; iD<nD; iD++) {
      slide = util->newSlide(TString::Format("res%s_vs_spill_%s",D[iD],_testBoard.c_str()),"");
      TH2F *hSpaceRes = new TH2F("hSpaceRes", Form("Mean %s residual vs spill %s",D[iD],_testBoard.c_str()), nSpills, _firstSpill-0.5, _finalSpill+0.5, 100, -0.01, 0.01);
      hSpaceRes->GetXaxis()->SetNoExponent();
      hSpaceRes->SetXTitle("Spill");
      hSpaceRes->SetStats(0);
      hSpaceRes->Draw();
      for(int iwbc=0; iwbc<nWBC; iwbc++) {
         if(h_resSpill[iwbc][iD][0]->GetEntries() == 0) continue;
         h_resSpill[iwbc][iD][0]->Draw("same");
      }
      slide->SaveAs(Form("%s/res%s_vs_spill_%s.%s",plotDir,D[iD],_testBoard.c_str(),plotExt));
      delete hSpaceRes;
   }

   for(int iwbc=0; iwbc<nWBC; iwbc++) {
      if(h_effMap[iwbc][0]->GetEntries() == 0) continue;

      sprintf(slideN,"efficiency_map_%s_WBC%d",_testBoard.c_str(),WBCvalue[iwbc]);
      sprintf(slideT,"Efficiency Map (WBC=%d) %s",WBCvalue[iwbc],_testBoard.c_str());
      sprintf(slideF,"%s/%s.%s",plotDir,slideN,plotExt);
      slide = new TCanvas(slideN, slideT, 0, 0, 1000, 600);
      slide->SetRightMargin(0.10);
      slide->SetBottomMargin(0.135);
      h_effMap[iwbc][2]->SetStats(0);
      h_effMap[iwbc][2]->SetMinimum(0.7);
      h_effMap[iwbc][2]->SetMaximum(1.);
      h_effMap[iwbc][2]->Draw("colz");
      slide->SaveAs(slideF);

      // TProfile* effMapProfX = h_effMap[iwbc][2]->ProfileX();
      // sprintf(slideN,"effMap_profX_%s_WBC%d",_testBoard.c_str(),WBCvalue[iwbc]);
      // sprintf(slideT,"Efficiency Map X Profile (WBC=%d) %s",WBCvalue[iwbc],_testBoard.c_str());
      // sprintf(slideF,"%s/%s.%s",plotDir,slideN,plotExt);
      // slide = new TCanvas(slideN, slideT, 0, 0, 1000, 600);
      // slide->SetRightMargin(0.10);
      // slide->SetBottomMargin(0.135);
      // effMapProfX->SetStats(0);
      // effMapProfX->Draw("e");
      // slide->SaveAs(slideF);

      // TProfile* effMapProfY = h_effMap[iwbc][2]->ProfileY();
      // sprintf(slideN,"effMap_profY_%s_WBC%d",_testBoard.c_str(),WBCvalue[iwbc]);
      // sprintf(slideT,"Efficiency Map Y Profile (WBC=%d) %s",WBCvalue[iwbc],_testBoard.c_str());
      // sprintf(slideF,"%s/%s.%s",plotDir,slideN,plotExt);
      // slide = new TCanvas(slideN, slideT, 0, 0, 1000, 600);
      // slide->SetRightMargin(0.10);
      // slide->SetBottomMargin(0.135);
      // effMapProfY->SetStats(0);
      // effMapProfY->Draw("e");
      // slide->SaveAs(slideF);      

      sprintf(slideN,"effTrack_map_%s_WBC%d",_testBoard.c_str(),WBCvalue[iwbc]);
      sprintf(slideT,"Track Map (WBC=%d) %s",WBCvalue[iwbc],_testBoard.c_str());
      sprintf(slideF,"%s/%s.%s",plotDir,slideN,plotExt);
      slide = new TCanvas(slideN, slideT, 0, 0, 1000, 600);
      slide->SetRightMargin(0.10);
      slide->SetBottomMargin(0.135);
      h_effMap[iwbc][0]->SetStats(0);
      h_effMap[iwbc][0]->Draw("colz");
      slide->SaveAs(slideF);

      sprintf(slideN,"efficiency_mapWide_%s_WBC%d",_testBoard.c_str(),WBCvalue[iwbc]);
      sprintf(slideT,"Efficiency MapWide (WBC=%d) %s",WBCvalue[iwbc],_testBoard.c_str());
      sprintf(slideF,"%s/%s.%s",plotDir,slideN,plotExt);
      slide = new TCanvas(slideN, slideT, 0, 0, 1000, 600);
      slide->SetRightMargin(0.10);
      slide->SetBottomMargin(0.135);
      h_effMapWide[iwbc][2]->SetStats(0);
      h_effMapWide[iwbc][2]->SetMinimum(0.7);
      h_effMapWide[iwbc][2]->SetMaximum(1.);
      h_effMapWide[iwbc][2]->Draw("colz");
      slide->SaveAs(slideF);

      sprintf(slideN,"effTrack_mapWide_%s_WBC%d",_testBoard.c_str(),WBCvalue[iwbc]);
      sprintf(slideT,"Track MapWide (WBC=%d) %s",WBCvalue[iwbc],_testBoard.c_str());
      sprintf(slideF,"%s/%s.%s",plotDir,slideN,plotExt);
      slide = new TCanvas(slideN, slideT, 0, 0, 1000, 600);
      slide->SetRightMargin(0.10);
      slide->SetBottomMargin(0.135);
      h_effMapWide[iwbc][0]->SetStats(0);
      h_effMapWide[iwbc][0]->Draw("colz");
      slide->SaveAs(slideF);

   }

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
   if(!setSlopeCut())
      return false;

   if(!setFiducialCut())
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
   tgaeTP->Divide(hitsTP,tracksTP,"cl=0.683 b(1,1) mode");

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

   int nSpills = 1+_finalSpill-_firstSpill;

   char suffix[16];
   for(int iwbc=0; iwbc<nWBC; iwbc++) {
      for(int i=0; i<3; i++) { //loop over histograms (usually tracks and hits, for efficiency)
         if(0==i)      sprintf(suffix, "tracks");
         else if(1==i) sprintf(suffix, "hits");
         else          sprintf(suffix, "efficiency");

         sprintf(name, "h_effFlux_wbc%d_%s",WBCvalue[iwbc],suffix);
         sprintf(title, "%s vs flux",suffix);
         h_effFlux[iwbc][i] = new TH1F(name, title, 200, 0., 1000.);
         h_effFlux[iwbc][i]->Sumw2();

         sprintf(name, "h_effNHits_wbc%d_%s",WBCvalue[iwbc],suffix);
         sprintf(title, "%s vs pixel hits",suffix);
         h_effNHits[iwbc][i] = new TH1F(name, title, 100, -0.5, 99.5);
         h_effNHits[iwbc][i]->Sumw2();

         sprintf(name, "h_effSpill_wbc%d_%s",WBCvalue[iwbc],suffix);
         sprintf(title, "%s vs spill",suffix);
         h_effSpill[iwbc][i] = new TH1F(name, title, nSpills, _firstSpill-0.5, _finalSpill+0.5);
         h_effSpill[iwbc][i]->Sumw2();

         sprintf(name, "h_effMap_wbc%d_%s",WBCvalue[iwbc],suffix);
         sprintf(title, "%s map (WBC %d)",suffix,WBCvalue[iwbc]);
         h_effMap[iwbc][i] = new TH2F(name, title,60,-4.5,4.5,90,-4.5,4.5);
         h_effMap[iwbc][i]->Sumw2();

         sprintf(name, "h_effMapWide_wbc%d_%s",WBCvalue[iwbc],suffix);
         sprintf(title, "%s map (WBC %d)",suffix,WBCvalue[iwbc]);
         h_effMapWide[iwbc][i] = new TH2F(name, title,30,-4.5,4.5,45,-4.5,4.5);
         h_effMapWide[iwbc][i]->Sumw2();

      }

      for(int i=0; i<2; i++) {//loop over residual histograms
         if(0==i) sprintf(suffix, "mean");
         else     sprintf(suffix, "sigma");

         for(int iD=0; iD<nD; iD++) {
            sprintf(name, "h_res%sSpill_wbc%d_%s",D[iD],WBCvalue[iwbc],suffix);
            sprintf(title, "%s residual %s vs spill (WBC %d)",D[iD],suffix,WBCvalue[iwbc]);
            h_resSpill[iwbc][iD][i] = new TH1F(name, title, nSpills, _firstSpill-0.5, _finalSpill+0.5);
         }
      }

   }//loop over WBC values

}

void tbAna::loadHistogramsFromFile(char* fname) {
   TFile *f = new TFile(fname);
   if(f->IsZombie()) {
      cout << "ERROR: File " << fname << " does not exist\n";
      return;
   }
   else {
      cout << "Loading histograms from " << fname << endl;
   }

   char suffix[16];
   for(int iwbc=0; iwbc<nWBC; iwbc++) {
      for(int i=0; i<3; i++) { //loop over histograms (usually tracks and hits, for efficiency)
         if(0==i)      sprintf(suffix, "tracks");
         else if(1==i) sprintf(suffix, "hits");
         else          sprintf(suffix, "efficiency");

         h_effFlux[iwbc][i] = (TH1F*) f->Get(Form("h_effFlux_wbc%d_%s",WBCvalue[iwbc],suffix));
         h_effNHits[iwbc][i] = (TH1F*) f->Get(Form("h_effNHits_wbc%d_%s",WBCvalue[iwbc],suffix));
         h_effSpill[iwbc][i] = (TH1F*) f->Get(Form("h_effSpill_wbc%d_%s",WBCvalue[iwbc],suffix));

         h_effMap[iwbc][i] = (TH2F*) f->Get(Form("h_effMap_wbc%d_%s",WBCvalue[iwbc],suffix));
         h_effMapWide[iwbc][i] = (TH2F*) f->Get(Form("h_effMapWide_wbc%d_%s",WBCvalue[iwbc],suffix));
      }

      for(int i=0; i<2; i++) {//loop over residual histograms
         if(0==i) sprintf(suffix, "mean");
         else     sprintf(suffix, "sigma");

         for(int iD=0; iD<nD; iD++) {
            h_resSpill[iwbc][iD][i] = (TH1F*) f->Get(Form("h_res%sSpill_wbc%d_%s",D[iD],WBCvalue[iwbc],suffix));
         }
      }

   }//loop over WBC values
   

}
