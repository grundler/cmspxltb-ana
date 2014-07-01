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

class plotter {
public:
   plotter(int spill0, int spill1, string outdir);
   ~plotter();

   void graphSetting(TGraph* g, TString name, TString title="", 
                     int color=-1, int style=-1, 
                     TString xTitle="", TString yTitle="");
   void graphSetting(TGraph* g, TH1F* h,
                     int color=-1, int style=-1, 
                     TString xTitle="", TString yTitle="");
   void spacerSetting(TH2* h, TGraph* g);
   TCanvas* newSlide(TString name, TString title="");
   void legendSetting(TLegend* l);

   void bookHistos();
   void makePlots();

   void writeFile(string filename);
   void loadHistogramsFromFile(char* fname);

   //Histograms and graphs

   //tracks, hits per flux to do efficiency, third histogram is filler for graph
   //vs flux
   TH1F *h_effFlux[nWBC][3]; 
   TH1F *h_effNHits[nWBC][3]; 
   TH2F *h_nHitsFlux[nWBC];

   //vs spill
   TH1F *h_effSpill[nWBC][3];

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