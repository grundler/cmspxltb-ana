#include <TCut.h>
#include <TChain.h>
#include <TMath.h>

#include "treeCorrelator.hh"

class TH1F;

using std::string;
using std::vector;

static const int nD = 2;
static const char D[nD][2] = {"X","Y"};
static const float resLo[nD] = {-0.5,-0.5};
static const float resHi[nD] = {0.5,0.5};

class tbAna {
public:
   
   tbAna(int dutID, string board, int spill0, int spill1);
   ~tbAna();
   
   void plotEffVsIntensity(TCut myCut="");

   bool initSpill(int spill);

   void bookHistos();

private:

   bool setSlopeCut();
   bool setFiducialCut();
   bool setTriggerPhaseCut();

   void initTrackTree(int spill);

   void loadTrackEntry(int entry);

   //Basic information
   const string _testBoard;
   const int _DUTID;
   const int _firstSpill;
   const int _finalSpill;

   //Info to correlate QIE and telescope data
   treeCorrelator* _tc;

   //Predefined cuts
   TCut chi2_50, chi2_600;
   TCut nearTrack;
   TCut slopeCut;
   TCut isFiducial;
   int  correctTriggerPhase;

   //Histograms and graphs

   //vs flux
   // TH1F *h_tracksIntensity, *h_hitsIntensity;
   TH1F *h_effFlux[2]; //tracks, hits per flux to do efficiency

   //vs spill
   // TH1F *h_tracksRun, *h_hitsRun;
   TH1F *h_effSpill[2]; //tracks, hits per spill to do spill-by-spill efficiency

   TH1F *h_resSpill[nD][2]; //mean and sigma of residuals per spill 

   //Tree information

   //track tree
   TChain *_trackTree;
   Int_t    RunNr, EvtNr;
   Double_t  fitX_0, fitY_0;
   Double_t  fitX_1, fitY_1;
   Double_t  fitX_2, fitY_2;
   Double_t  fitX_3, fitY_3;
   Double_t  fitX_4, fitY_4;
   Double_t  fitX_5, fitY_5;
   Double_t  fitX_6, fitY_6;
   Double_t  fitX_7, fitY_7;
   Double_t  dutX, dutY;
   Int_t Ndf;
   Float_t Chi2;
   Int_t nTrack;
   Double_t nhits_7;

};

Double_t g(Double_t *v, Double_t *par)
{
   Double_t arg = 0;
   arg = (v[0] - par[1]) / par[2];

   Double_t area = par[0]/(par[2]*sqrt(2*3.1415926)) * par[3] ;
   //Double_t norm_area = 1/(par[2]*sqrt(2*3.1415926));                                                                
   Double_t gaus = area*TMath::Exp(-0.5*arg*arg);

   if (gaus<=0) gaus=1e-10;
   return gaus;
}
