#include "treeCorrelator.hh"
#include "constants.hh"

#include <iostream>
#include <sstream>
#include <iomanip>
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TH2D.h"
#include <TEntryList.h>
#include <string>


treeCorrelator::treeCorrelator(int spill, string board) 
   : _spill(spill),
     _board(board),
     _isInitialized(false),
     _wbc(0),
     _nBuckets(-1)
{

   _syncMap.clear();
   _tpMap.clear();
   _fluxMap.clear();

   cout << "\tInitializing maps for spill " << _spill << endl;

   ostringstream stream;
   stream << subdir << "/" << _board << "/qie/summary_" << _spill << ".root";
   string QIEfilename = stream.str();
   TFile *f = new TFile(QIEfilename.c_str());
   if(f->IsZombie()) {
      std::cout << "File " << QIEfilename << " does not exist\n";
      return;
   }
   TH1F *hWBC = (TH1F*)f->Get("hsummary_WBC");
   if(hWBC)
      _wbc = hWBC->GetBinContent(1);
   else {
      cout << "Couldn't get WBC. Skipping\n";
      return;
   }
   // _wbc = ((TH1F*)f->Get("hsummary_WBC"))->GetBinContent(1);
   TH1F *hBuckets = (TH1F*)f->Get("hsummary_Buckets_in_intensity");
   if(hBuckets)
      _nBuckets = hBuckets->GetBinContent(1);
   else {
      cout << "Couldn't get number of buckets. Skipping\n";
      return;
   }
   cout << "nBuckets = " << _nBuckets << endl;
   //_nBuckets = ((TH1F*)f->Get("hsummary_Buckets_in_intensity"))->GetBinContent(1);

   stream.str("");
   stream << subdir << "/" << _board << "/timestamps/" << _spill << "_tp.root";
   string TBfilename = stream.str();
   f = new TFile(TBfilename.c_str());
   if(f->IsZombie()) {
      std::cout << "File " << TBfilename << " does not exist\n";
      return;
   }

   tree_summary = new TChain("tree_summary");
   tree_summary->Add(QIEfilename.c_str());
   
   summary_Trigger_count = -1;
   summary_Trigger_RF_onset = 0;
   summary_Trigger_turn_onset = 0;
   summary_Trigger_BeamIntensity_itself = 0; // itself
   summary_Trigger_BeamIntensity_minusWBC = 0;
   summary_Trigger_BeamIntensity = 0;
   summary_Trigger_Nproton_itself = 0.; // itself
   summary_Trigger_Nproton_minusWBC = 0.;
   summary_Trigger_Nproton = 0.;
   tree_summary->SetBranchAddress("summary_Trigger_count",&summary_Trigger_count);
   tree_summary->SetBranchAddress("summary_Trigger_RF_onset",&summary_Trigger_RF_onset);
   tree_summary->SetBranchAddress("summary_Trigger_turn_onset",&summary_Trigger_turn_onset);
   tree_summary->SetBranchAddress("summary_Trigger_BeamIntensity_itself",&summary_Trigger_BeamIntensity_itself);
   tree_summary->SetBranchAddress("summary_Trigger_BeamIntensity_minusWBC",&summary_Trigger_BeamIntensity_minusWBC);
   tree_summary->SetBranchAddress("summary_Trigger_BeamIntensity",&summary_Trigger_BeamIntensity);
   tree_summary->SetBranchAddress("summary_Trigger_Nproton_itself",&summary_Trigger_Nproton_itself);
   tree_summary->SetBranchAddress("summary_Trigger_Nproton_minusWBC",&summary_Trigger_Nproton_minusWBC);
   tree_summary->SetBranchAddress("summary_Trigger_Nproton",&summary_Trigger_Nproton);

   tree = new TChain("tree");
   tree->Add(TBfilename.c_str());

   EventNumber = 0;
   TimeStamp = 0;
   TriggerPhase = 0;
   tree->SetBranchAddress("EventNumber",&EventNumber);
   tree->SetBranchAddress("TimeStamp",&TimeStamp);
   tree->SetBranchAddress("TriggerPhase",&TriggerPhase);

   createMapLong();
   _isInitialized = true;   
}

void treeCorrelator::createMapLong() {

   double timebin = 1000./53.104;  //[us] <-- 18.83ns

   int previousTime = 0;
   int currentTime = 0;
   int nextTime = 0;
   int TimeStampOffset = -1;
   for(int iOffset=0;iOffset<tree->GetEntries();iOffset++) {
      tree->GetEntry(iOffset);
      if(iOffset==0) {
         previousTime = TimeStamp;
         tree->GetEntry(tree->GetEntries()-1);
         if (TimeStamp - previousTime < 5000000){
            TimeStampOffset = TimeStamp;
            break;
         }
      }

      _tpMap[EventNumber] = TriggerPhase;

      if(fabs(TimeStamp - previousTime)>10000000&& TimeStamp>0 && previousTime>0 && iOffset<0.1*tree->GetEntries()) {
         currentTime = TimeStamp;
         if (iOffset+1<tree->GetEntries()) {
            tree->GetEntry(iOffset+1);
            nextTime = TimeStamp;
         }
         tree->GetEntry(iOffset);
         if(nextTime-currentTime< 30) {
            TimeStampOffset = TimeStamp;
            //break;
         }
      }
      previousTime = TimeStamp;
   }

   double TrimTimestamp_ = -1;
   double previousInformation[4] = {0.0, 0, 0.0, 0.0};    // [TBtimestamp, QIEindex, Trigger_turn_onset, Trigger_RF_onset]
   int previousDiff = -999;
   for(int itimestamp=0;itimestamp<tree->GetEntries();itimestamp++) {
      tree->GetEntry(itimestamp);
      if(itimestamp%1000==0) 
         printf("\t\ttreeCorrelator on %1.3f %% (%i th : %i)\r",100.*itimestamp/tree->GetEntries(), itimestamp, TimeStamp);
      if(TimeStampOffset!=-1&&(TimeStamp - TimeStampOffset)>=0) {
          
         if(TimeStamp - TimeStampOffset==0){
            tree_summary->GetEntry(0);
            previousInformation[0] = 0;
            previousInformation[1] = 0;
            previousInformation[2] = summary_Trigger_turn_onset;
            previousInformation[3] = summary_Trigger_RF_onset;
         }
          
         int currentDiff = -999;
         for(int iQIE=(int)previousInformation[1];iQIE<tree_summary->GetEntries();iQIE++) {
            tree_summary->GetEntry(iQIE);
            TrimTimestamp_ = 
               ( (summary_Trigger_turn_onset - previousInformation[2])*588. 
                 + (summary_Trigger_RF_onset - previousInformation[3]) )*timebin/1000. 
               + previousInformation[0];
            if(fabs((TimeStamp - TimeStampOffset) - TrimTimestamp_) < 2) {
               currentDiff = EventNumber - summary_Trigger_count;
               _fluxMap[summary_Trigger_count] = calcFlux();

               previousInformation[0] = TimeStamp - TimeStampOffset;
               previousInformation[1] = iQIE;
               previousInformation[2] = summary_Trigger_turn_onset;
               previousInformation[3] = summary_Trigger_RF_onset;
               break;
            }
         }
          
         if(currentDiff != previousDiff) {
            _syncMap.insert( std::pair<int, int>(EventNumber, currentDiff) );
            previousDiff = currentDiff;
         }          
      }
   }
   printf("\n");

   //_isInitialized = true;
    
   // for(map<int,int>::iterator _syncMapItr = _syncMap.begin(); _syncMapItr != _syncMap.end(); ++_syncMapItr)
   //    cout << _syncMapItr->first << " " << _syncMapItr->second << endl;

}

treeCorrelator::~treeCorrelator() {
   delete tree;
   delete tree_summary;
   _syncMap.clear();
   _tpMap.clear();
   _fluxMap.clear();

}

int treeCorrelator::getQieEvent(int EvtNumr) {

   map<int,int>::iterator _syncMapItr = _syncMap.begin();
   int TrigCount = -1;//EvtNumr - TimestampMapItr_->second.second;
   // cout << " EvtNumr = " << EvtNumr << endl;
   while(_syncMapItr!=_syncMap.end() && EvtNumr>=_syncMapItr->first) {
      // cout << "changeup = " << TimestampMapItr_->second.first;
      TrigCount = EvtNumr - _syncMapItr->second;
      // cout << " TrigCount = " << TrigCount << endl;
      _syncMapItr++;
   }
   
   return TrigCount;
}

int treeCorrelator::getTriggerPhase(int EvtNumr) {
   //cout << " Get Trigger Phase for event " << EvtNumr << " in map of size " << _tpMap.size() << endl;
   map<int, int>::iterator it;
   it = _tpMap.find(EvtNumr);
   if(it != _tpMap.end())
      return it->second;

   return -1;
}

float treeCorrelator::getFlux(int trigCount) {
   map<int, float>::iterator it;
   it = _fluxMap.find(trigCount);
   if(it != _fluxMap.end())
      return it->second;

   return -1.;
}

float treeCorrelator::calcFlux() {

   float qb6 = 0.6755; //qie coverage of beam in x (6mm width)
   float qb5 = 0.7517; //qie coverage of beam in y (5mm width)

   float rbx = 0.495; //fraction of beam seen by roc in x
   float rby = 0.576; //fraction of beam seen by roc in y
   float dx = 0.8; //width of roc in x
   float dy = 0.8; //width of roc in y

   float nproton = summary_Trigger_Nproton 
      + summary_Trigger_Nproton_minusWBC
      - summary_Trigger_Nproton_itself;

   nproton /= (_nBuckets*2.-1.);

   //get flux by normalizing the number of protons by FNAL clock, qie scintillator coverage, and roc coverage
   float flux = nproton * (fnalClock/qb6/qb5) * (rbx*rby/dx/dy);

   return flux;
}
