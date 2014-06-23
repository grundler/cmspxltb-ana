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
#include <vector>

using std::vector;

treeCorrelator::treeCorrelator(int spill, string board, int algo) 
   : _spill(spill),
     _board(board),
     _algo(algo),
     _isInitialized(false),
     _wbc(0),
     _nBuckets(-1)
{

   _syncMap.clear();
   _tpMap.clear();
   _fluxMap.clear();
   _maxMap.clear();

   cout << "\tInitializing maps for spill " << _spill << endl;

   initTpTree();
   if(tree==NULL) {
      cout << "\tCould not get trigger phase information for spill " << _spill << ". Skipping\n";
      return;
   }

   initQieTree();
   if(tree_summary==NULL) {
      cout << "\tCould not get QIE information for spill " << _spill << ". Skipping\n";
      return;
   }

   if(_algo == 0)
      createMapLong();
   else if(_algo == 1)
      createMapQuick();
   else {
      cout << "ERROR: unrecognized option for mapping event numbers\n";
      return;
   }

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
      _tpMap[EventNumber] = TriggerPhase;

      // if(itimestamp%1000==0) 
      //    printf("\t\ttreeCorrelator on %1.3f %% (%i th : %i)\r",100.*itimestamp/tree->GetEntries(), itimestamp, TimeStamp);
      if(TimeStampOffset!=-1&&(TimeStamp - TimeStampOffset)>=0) {
          
         if(TimeStamp - TimeStampOffset==0) {
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
               _fluxMap[summary_Trigger_count] =  summary_Trigger_Nproton 
                                                + summary_Trigger_Nproton_minusWBC 
                                                - summary_Trigger_Nproton_itself;//calcFlux();
               _maxMap[summary_Trigger_count] = summary_Trigger_Nproton_maximum;

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
   // cout << endl;

   // cout << "Synchronization map\n";
   // for(map<int,int>::iterator _syncMapItr = _syncMap.begin(); _syncMapItr != _syncMap.end(); ++_syncMapItr)
   //    cout << _syncMapItr->first << " " << _syncMapItr->second << endl;
   // for(map<int,float>::iterator _fluxMapItr = _fluxMap.begin(); _fluxMapItr != _fluxMap.end(); ++_fluxMapItr)
   //    cout << _fluxMapItr->first << " " << _fluxMapItr->second << endl;

}

void treeCorrelator::createMapQuick() {

   int initEvent = 300; //event to test against when finding first event of new spill

   int spillGap = 100000; //time difference large enough to be considered different spills
   int evtGap = 12; //time difference large enough to note skipped event

   int turnGap = 1; //difference large enough to note skipped turn in QIE

   //Find first event of new spill
   tree->GetEntry(initEvent);
   int time_ref = TimeStamp;
   int evt_t0 = 0;
   // cout << "Reference event: " << EventNumber << ", time: " << time_ref << endl;
   for(int iEvt=0; iEvt<initEvent; iEvt++) {
      tree->GetEntry(iEvt);
      // cout << " iEvt " << iEvt << ", event " << EventNumber << ", time " << TimeStamp << endl;
      if((time_ref-TimeStamp)>spillGap)
         evt_t0 = EventNumber;
   }
   // cout << " Real first event of spill is " << evt_t0 << endl;
   
   //get list of events where gap appears
   vector<int> tb_evtnr;
   for(int iEvt=0; iEvt<tree->GetEntries(); iEvt++) {
      tree->GetEntry(iEvt);
      _tpMap[EventNumber] = TriggerPhase; //make map of trigger phase while we're at it

      if(dtime>evtGap) {
         tb_evtnr.push_back(EventNumber);
         // cout << " dtime > " << evtGap << ", EvtNr " << EventNumber << ", dtime " << dtime << endl;
      }
   }

   //get list of qie events
   int turn_p = 0;
   vector<int> qie_evtnr;
   for(int iQIE=0; iQIE<tree_summary->GetEntries(); iQIE++) {
      tree_summary->GetEntry(iQIE);
      _fluxMap[summary_Trigger_count] =  summary_Trigger_Nproton 
                                       + summary_Trigger_Nproton_minusWBC 
                                       - summary_Trigger_Nproton_itself;//calcFlux();
      _maxMap[summary_Trigger_count] = summary_Trigger_Nproton_maximum;

      int dturn = summary_Trigger_turn_onset-turn_p;
      if(dturn>turnGap || summary_Trigger_turn_onset==1) {
         qie_evtnr.push_back(summary_Trigger_count);
         // cout << " QIE trig " << summary_Trigger_count << ", turn " << summary_Trigger_turn_onset << ", dturn " << dturn << endl;
      }
      turn_p = summary_Trigger_turn_onset;
   }

   //correlate lists
   int notFound = 9999;
   int previousDiff = -notFound;
   for(unsigned i=0; i<tb_evtnr.size(); i++) {
      int mindiff = notFound;
      for(unsigned j=0; j<qie_evtnr.size(); j++) {
         int diff = tb_evtnr[i] - qie_evtnr[j];
         if(abs(diff) < abs(mindiff)) mindiff = diff;
      }
      // cout << "  Event " << tb_evtnr[i] << ", diff " << mindiff << endl;
      if(mindiff!=previousDiff) {
         if(mindiff == notFound) {
            mindiff = previousDiff + 1;
         }
         _syncMap.insert( std::pair<int, int>(tb_evtnr[i], mindiff) );
         previousDiff = mindiff;
      }
   }

   // cout << endl;

   // cout << "Synchronization map\n";
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

   // cout << " Event " << EvtNumr << " corresponds to QIE trigger " << TrigCount << endl;
   
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
      return calcFlux(it->second);

   return -1.;
}

float treeCorrelator::getFluxRatio(int trigCount) {
   float fluxMax = -1.;
   float flux = 1.;

   map<int, float>::iterator it;

   it = _fluxMap.find(trigCount);
   if(it != _fluxMap.end()) 
      flux = it->second;

   it = _maxMap.find(trigCount);
   if(it != _maxMap.end()) 
      fluxMax = it->second;

   return (flux > 0. ? fluxMax/(flux/(_nBuckets*2.-1.)) : -1.);
}

float treeCorrelator::calcFlux(float nproton) {

   float qb6 = 0.6755; //qie coverage of beam in x (6mm width)
   float qb5 = 0.7517; //qie coverage of beam in y (5mm width)

   float rbx = 0.495; //fraction of beam seen by roc in x
   float rby = 0.576; //fraction of beam seen by roc in y
   float dx = 0.8; //width of roc in x
   float dy = 0.8; //width of roc in y

   // float nproton = summary_Trigger_Nproton 
   //    + summary_Trigger_Nproton_minusWBC
   //    - summary_Trigger_Nproton_itself;

   nproton /= (_nBuckets*2.-1.);

   //get flux by normalizing the number of protons by FNAL clock, qie scintillator coverage, and roc coverage
   float flux = nproton * (fnalClock/qb6/qb5) * (rbx*rby/dx/dy);

   return flux;
}

void treeCorrelator::initTpTree() {

   tree = NULL;

   //Open file
   ostringstream stream;
   stream << subdir << "/" << _board << "/timestamps/" << _spill << "_tp.root";
   string TPfilename = stream.str();
   TFile *f = new TFile(TPfilename.c_str());
   if(f->IsZombie()) {
      cout << "File " << TPfilename << " does not exist\n";
      return;
   }

   //Load tree
   tree = new TChain("tree");

   tree->SetBranchStatus("*",kFALSE);

   tree->SetBranchStatus("EventNumber",kTRUE);
   tree->SetBranchStatus("TimeStamp",kTRUE);
   tree->SetBranchStatus("dtime",kTRUE);
   tree->SetBranchStatus("TriggerPhase",kTRUE);

   tree->SetBranchAddress("EventNumber",&EventNumber);
   tree->SetBranchAddress("TimeStamp",&TimeStamp);
   tree->SetBranchAddress("dtime",&dtime);
   tree->SetBranchAddress("TriggerPhase",&TriggerPhase);

   tree->Add(TPfilename.c_str());
}

void treeCorrelator::initQieTree() {

   tree_summary = NULL;

   //Open file
   ostringstream stream;
   stream << subdir << "/" << _board << "/qie/summary_" << _spill << ".root";
   string QIEfilename = stream.str();
   TFile *f = new TFile(QIEfilename.c_str());
   if(f->IsZombie()) {
      cout << "File " << QIEfilename << " does not exist\n";
      return;
   }

   //Get WBC and number of buckets while the file is open
   TH1F *hWBC = (TH1F*)f->Get("hsummary_WBC");
   if(hWBC)
      _wbc = hWBC->GetBinContent(1);
   else {
      cout << "Couldn't get WBC. Skipping\n";
      return;
   }
   TH1F *hBuckets = (TH1F*)f->Get("hsummary_Buckets_in_intensity");
   if(hBuckets)
      _nBuckets = hBuckets->GetBinContent(1);
   else {
      cout << "Couldn't get number of buckets. Skipping\n";
      return;
   }
   cout << "nBuckets = " << _nBuckets << endl;

   //Load tree
   tree_summary = new TChain("tree_summary");

   tree_summary->SetBranchStatus("*",kFALSE);

   tree_summary->SetBranchStatus("summary_Trigger_count",kTRUE);
   tree_summary->SetBranchStatus("summary_Trigger_RF_onset",kTRUE);
   tree_summary->SetBranchStatus("summary_Trigger_turn_onset",kTRUE);
   tree_summary->SetBranchStatus("summary_Trigger_Nproton_itself",kTRUE);
   tree_summary->SetBranchStatus("summary_Trigger_Nproton_minusWBC",kTRUE);
   tree_summary->SetBranchStatus("summary_Trigger_Nproton",kTRUE);
   tree_summary->SetBranchStatus("summary_Trigger_Nproton_maximum",kTRUE);

   summary_Trigger_count = -1;
   tree_summary->SetBranchAddress("summary_Trigger_count",&summary_Trigger_count);
   tree_summary->SetBranchAddress("summary_Trigger_RF_onset",&summary_Trigger_RF_onset);
   tree_summary->SetBranchAddress("summary_Trigger_turn_onset",&summary_Trigger_turn_onset);
   tree_summary->SetBranchAddress("summary_Trigger_Nproton_itself",&summary_Trigger_Nproton_itself);
   tree_summary->SetBranchAddress("summary_Trigger_Nproton_minusWBC",&summary_Trigger_Nproton_minusWBC);
   tree_summary->SetBranchAddress("summary_Trigger_Nproton",&summary_Trigger_Nproton);
   tree_summary->SetBranchAddress("summary_Trigger_Nproton_maximum",&summary_Trigger_Nproton_maximum);

   tree_summary->Add(QIEfilename.c_str());
}
