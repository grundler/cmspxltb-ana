#include "utils.hh"

#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH2.h>
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

void utils::spacerSetting(TH2* h, TGraph* g) {
   h->SetXTitle(g->GetXaxis()->GetTitle());
   h->SetYTitle(g->GetYaxis()->GetTitle());

   h->GetYaxis()->SetTitleOffset(0.8);
   h->SetStats(0);
}

TCanvas* utils::newSlide(TString name, TString title) {
   TCanvas* slide = new TCanvas(name, title, 0, 0, 1000, 600);

   slide->SetRightMargin(0.10);
   slide->SetBottomMargin(0.135);

   slide->SetGrid();

   return slide;
}

void utils::legendSetting(TLegend* l) {
   l->SetFillColor(10);
   l->SetBorderSize(0);
}
