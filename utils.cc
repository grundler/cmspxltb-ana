#include "utils.hh"

#include <TGraph.h>
#include <TAxis.h>

void utils::graphSetting(TGraph* g, TString name, TString title, 
                         int color, int style, 
                         TString xTitle, TString yTitle) {

   Color_t mColor = color>=0 ? color : 1; //black as default
   Style_t mStyle = style>=0 ? style : 8; //scalable dot as default

   g->SetName(name);
   g->SetTitle(title);
   g->SetMarkerColor(mColor);
   g->SetLineColor(mColor);
   g->SetMarkerStyle(mStyle);
   g->SetLineStyle(mStyle);
   g->GetXaxis()->SetTitle(xTitle);
   g->GetYaxis()->SetTitle(yTitle);

}
