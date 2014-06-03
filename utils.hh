#include <TString.h>

class TGraph;

class utils {
public:
   utils() {};
   ~utils() {};

   void graphSetting(TGraph* g, TString name, TString title="", 
                     int color=-1, int style=-1, 
                     TString xTitle="", TString yTitle="");
};
