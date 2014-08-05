#ifndef PLOTTER_HH
#define PLOTTER_HH

#include <TString.h>
#include "constants.hh"

class TCanvas;
class TH1F;
class TH2;
class TH2F;
class TGraph;
class TLegend;

using std::string;
using std::vector;

class plotter {
public:
   plotter(int spill0, int spill1, string outdir);
   ~plotter();

   TH1F* graphSetting(TGraph* g, TString name, TString title="", 
                     int color=-1, int style=-1, 
                     TString xTitle="", TString yTitle="");
   void graphSetting(TGraph* g, TH1F* h,
                     int color=-1, int style=-1, 
                     TString xTitle="", TString yTitle="");

   TCanvas* newSlide(TString name, TString title="");
   void legendSetting(TLegend* l);

   void bookHistos();
   void makePlots(int startWBC=wbc99, int finalWBC=wbc255);
   void compareSpills(int wbc, vector<int> spillList);

   void writeFile(string filename);
   void loadHistogramsFromFile(char* fname);

   //Histograms and graphs

   //tracks, hits per flux to do efficiency

   //vs flux
   TH2F *h_effFlux[nWBC][2]; 
   TH2F *h_effNHits[nWBC][2]; 
   TH2F *h_nHitsFlux[nWBC];

   //vs spill
   TH1F *h_effSpill[nWBC][2];

   TH1F *h_resSpill[nWBC][nD][2]; //mean and sigma of residuals per spill 
   TH2F *h_tpSpill[nWBC][2]; //for efficency vs trigger phase for each spill

   //geometrical
   TH2F *h_effMap[nWBC][3]; //tracks/hits/efficiency maps
   TH2F *h_effMapWide[nWBC][3];

private:
   const int _firstSpill;
   const int _finalSpill;
   const int _nSpills;
   string _outDir;

};

#endif // #ifndef PLOTTER_HH
