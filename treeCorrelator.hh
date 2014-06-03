/*
 * ########################################################################## 
 * Maps to correlate different event trees 
 * ########################################################################## 
 *
 * - Author : Ulysses Grundler
 * - Based on LoadSynchronizationMap.h by Jacky Tzeng
 *
 * */

#ifndef treeCorrelator_hh 
#define treeCorrelator_hh

#include <map>

using std::string;
using std::map;

class TChain;

class treeCorrelator{
public:
   treeCorrelator(int spill, string board, int algo=0);
   ~treeCorrelator();

   bool isInitialized() { return _isInitialized; };
   
   int getQieEvent(int EvtNumr);
   int getTriggerPhase(int EvtNumr);
   float getFlux(int trigCount);
   int getWBC() { return _wbc; };

private:
   const int    _spill;
   const string _board;
   const int    _algo;
   bool _isInitialized;

   int _wbc;
   int _nBuckets;

   void createMapLong();
   void createMapQuick();

   float calcFlux();

   TChain *tree;
   
   int EventNumber;
   int TimeStamp;
   int dtime;
   int TriggerPhase;

   TChain *tree_summary;

   int summary_Trigger_count;
   int summary_Trigger_RF_onset;
   int summary_Trigger_turn_onset;
   int summary_Trigger_BeamIntensity_itself;
   int summary_Trigger_BeamIntensity_minusWBC;
   int summary_Trigger_BeamIntensity;
   float summary_Trigger_Nproton_itself;
   float summary_Trigger_Nproton_minusWBC;
   float summary_Trigger_Nproton;


   //map eutelescope event and qie event
   map<int, int> _syncMap;

   //map event and trigger phase
   map<int, int> _tpMap;

   //map qie event and intensity
   map<int, float> _fluxMap;
   
};


#endif
