#include <TString.h>

class TCanvas;
class TH1F;
class TH2;
class TGraph;
class TLegend;

class utils {
public:
   utils() {};
   ~utils() {};

   void graphSetting(TGraph* g, TString name, TString title="", 
                     int color=-1, int style=-1, 
                     TString xTitle="", TString yTitle="");
   void graphSetting(TGraph* g, TH1F* h,
                     int color=-1, int style=-1, 
                     TString xTitle="", TString yTitle="");
   void spacerSetting(TH2* h, TGraph* g);
   TCanvas* newSlide(TString name, TString title="");
   void legendSetting(TLegend* l);

};
