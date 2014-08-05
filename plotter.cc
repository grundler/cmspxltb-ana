#include "plotter.hh"

#include <iostream>

#include <TFile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TAxis.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>

using std::cout;

plotter::plotter(int spill0, int spill1, string outdir):
   _firstSpill(spill0),
   _finalSpill(spill1),
   _nSpills(1+_finalSpill-_firstSpill),
   _outDir(outdir)
{
   bookHistos();
}

plotter::~plotter() {
   for(int iwbc=0; iwbc<nWBC; iwbc++) {
      for(int i=0; i<3; i++) {
         delete h_effMap[iwbc][i];
         delete h_effMapWide[iwbc][i];
         
         if(i<2) {
            delete h_effFlux[iwbc][i];
            delete h_effNHits[iwbc][i];
            delete h_effSpill[iwbc][i];
            
            delete h_tpSpill[iwbc][i];
            for(int j=0; j<nD; j++)
               delete h_resSpill[iwbc][j][i];
         }
      }
      delete h_nHitsFlux[iwbc];
   }
}

TH1F* plotter::graphSetting(TGraph* g, TString name, TString title, 
                         int color, int style, 
                         TString xTitle, TString yTitle) {

   Color_t mColor = color>=0 ? color : 1; //black as default
   Style_t mStyle = style>=0 ? style : 8; //scalable dot as default

   g->SetMarkerColor(mColor);
   g->SetLineColor(mColor);
   g->SetMarkerStyle(mStyle);
   g->SetLineStyle(mStyle);

   TH1F* h = new TH1F(name,title,100,g->GetXaxis()->GetXmin(),g->GetXaxis()->GetXmax());
   h->GetXaxis()->SetTitle(xTitle);
   h->GetYaxis()->SetTitle(yTitle);
   h->GetYaxis()->SetTitleOffset(0.8);
   h->SetStats(0);
   h->SetMaximum(1.02);
   h->SetMinimum(0.8);
   h->SetNdivisions(110,"Y");
   h->SetNdivisions(110,"X");
   h->SetLineWidth(1);

   g->SetHistogram(h);

   return h;
}

void plotter::graphSetting(TGraph* g, TH1F* h,
                         int color, int style, 
                         TString xTitle, TString yTitle) {

   Color_t mColor = color>=0 ? color : 1; //black as default
   Style_t mStyle = style>=0 ? style : 8; //scalable dot as default

   g->SetMarkerColor(mColor);
   g->SetLineColor(mColor);
   g->SetMarkerStyle(mStyle);
   g->SetLineStyle(mStyle);

   h->GetXaxis()->SetTitle(xTitle);
   h->GetYaxis()->SetTitle(yTitle);
   h->GetYaxis()->SetTitleOffset(0.8);
   h->SetStats(0);
   h->SetMaximum(1.02);
   h->SetMinimum(0.8);
   h->SetNdivisions(110,"Y");
   h->SetNdivisions(110,"X");
   h->SetLineWidth(1);

   g->SetHistogram(h);
}


TCanvas* plotter::newSlide(TString name, TString title) {
   TCanvas* slide = new TCanvas(name, title, 0, 0, 700, 500);

   slide->SetLeftMargin(0.10);
   slide->SetBottomMargin(0.135);

   slide->SetGrid();

   return slide;
}

void plotter::legendSetting(TLegend* l) {
   l->SetFillColor(10);
   l->SetBorderSize(0);

   l->SetTextFont(72);
   l->SetTextSize(0.04);
}


void plotter::bookHistos() {

   char title[128], name[256];

   char suffix[16];
   for(int iwbc=0; iwbc<nWBC; iwbc++) {
      for(int i=0; i<3; i++) { //loop over histograms (usually tracks and hits, for efficiency)
         if(0==i)      sprintf(suffix, "tracks");
         else if(1==i) sprintf(suffix, "hits");
         else          sprintf(suffix, "efficiency");

         sprintf(name, "h_effMap_wbc%d_%s",WBCvalue[iwbc],suffix);
         sprintf(title, "%s map (WBC %d)",suffix,WBCvalue[iwbc]);
         h_effMap[iwbc][i] = new TH2F(name, title,60,-4.5,4.5,90,-4.5,4.5);
         h_effMap[iwbc][i]->Sumw2();

         sprintf(name, "h_effMapWide_wbc%d_%s",WBCvalue[iwbc],suffix);
         sprintf(title, "%s map (WBC %d)",suffix,WBCvalue[iwbc]);
         h_effMapWide[iwbc][i] = new TH2F(name, title,30,-4.5,4.5,45,-4.5,4.5);
         h_effMapWide[iwbc][i]->Sumw2();

         if(i>1) continue;

         sprintf(name, "h_effFlux_wbc%d_%s",WBCvalue[iwbc],suffix);
         sprintf(title, "%s vs flux",suffix);
         h_effFlux[iwbc][i] = new TH2F(name, title, _nSpills, _firstSpill-0.5, _finalSpill+0.5, 200, 0., 1000.);
         h_effFlux[iwbc][i]->Sumw2();

         sprintf(name, "h_effNHits_wbc%d_%s",WBCvalue[iwbc],suffix);
         sprintf(title, "%s vs pixel hits",suffix);
         h_effNHits[iwbc][i] = new TH2F(name, title, _nSpills, _firstSpill-0.5, _finalSpill+0.5, 100, -0.5, 99.5);
         h_effNHits[iwbc][i]->Sumw2();

         sprintf(name, "h_effSpill_wbc%d_%s",WBCvalue[iwbc],suffix);
         sprintf(title, "%s vs spill",suffix);
         h_effSpill[iwbc][i] = new TH1F(name, title, _nSpills, _firstSpill-0.5, _finalSpill+0.5);
         h_effSpill[iwbc][i]->Sumw2();


         sprintf(name, "h_tpSpill_wbc%d_%s",WBCvalue[iwbc],suffix);
         sprintf(title, "%s - spill/trigger phase (WBC %d)",suffix,WBCvalue[iwbc]);
         h_tpSpill[iwbc][i] = new TH2F(name, title,_nSpills,_firstSpill-0.5,_finalSpill+0.5,nPhases,-0.5,nPhases-0.5);
         h_tpSpill[iwbc][i]->Sumw2();

      }

      for(int i=0; i<2; i++) {//loop over residual histograms
         if(0==i) sprintf(suffix, "mean");
         else     sprintf(suffix, "sigma");

         for(int iD=0; iD<nD; iD++) {
            sprintf(name, "h_res%sSpill_wbc%d_%s",D[iD],WBCvalue[iwbc],suffix);
            sprintf(title, "%s residual %s vs spill (WBC %d)",D[iD],suffix,WBCvalue[iwbc]);
            h_resSpill[iwbc][iD][i] = new TH1F(name, title, _nSpills, _firstSpill-0.5, _finalSpill+0.5);
         }
      }

      sprintf(name, "h_nHitsFlux_wbc%d",WBCvalue[iwbc]);
      sprintf(title, "hits vs flux (WBC %d)", WBCvalue[iwbc]);
      h_nHitsFlux[iwbc] = new TH2F(name, title, 200, 0., 1000., 100, -0.5, 99.5);
      h_nHitsFlux[iwbc]->SetXTitle("Flux");
      h_nHitsFlux[iwbc]->SetYTitle("Hits");
      h_nHitsFlux[iwbc]->Sumw2();

   }//loop over WBC values

}

void plotter::loadHistogramsFromFile(char* fname) {
   TFile *f = new TFile(fname);
   if(f->IsZombie()) {
      cout << "ERROR: File " << fname << " does not exist\n";
      return;
   }
   else {
      cout << "Loading histograms from " << fname << endl;
   }

   char suffix[16];
   for(int iwbc=0; iwbc<nWBC; iwbc++) {
      for(int i=0; i<3; i++) { //loop over histograms (usually tracks and hits, for efficiency)
         if(0==i)      sprintf(suffix, "tracks");
         else if(1==i) sprintf(suffix, "hits");
         else          sprintf(suffix, "efficiency");

         h_effMap[iwbc][i] = (TH2F*) f->Get(Form("h_effMap_wbc%d_%s",WBCvalue[iwbc],suffix));
         h_effMapWide[iwbc][i] = (TH2F*) f->Get(Form("h_effMapWide_wbc%d_%s",WBCvalue[iwbc],suffix));

         if(i>1) continue;
         h_effFlux[iwbc][i] = (TH2F*) f->Get(Form("h_effFlux_wbc%d_%s",WBCvalue[iwbc],suffix));
         h_effNHits[iwbc][i] = (TH2F*) f->Get(Form("h_effNHits_wbc%d_%s",WBCvalue[iwbc],suffix));
         h_effSpill[iwbc][i] = (TH1F*) f->Get(Form("h_effSpill_wbc%d_%s",WBCvalue[iwbc],suffix));

         h_tpSpill[iwbc][i] = (TH2F*) f->Get(Form("h_tpSpill_wbc%d_%s",WBCvalue[iwbc],suffix));
      }

      for(int i=0; i<2; i++) {//loop over residual histograms
         if(0==i) sprintf(suffix, "mean");
         else     sprintf(suffix, "sigma");

         for(int iD=0; iD<nD; iD++) {
            h_resSpill[iwbc][iD][i] = (TH1F*) f->Get(Form("h_res%sSpill_wbc%d_%s",D[iD],WBCvalue[iwbc],suffix));
         }
      }

      h_nHitsFlux[iwbc] = (TH2F*) f->Get(Form("h_nHitsFlux_wbc%d",WBCvalue[iwbc]));

   }//loop over WBC values
   

}

void plotter::writeFile(string filename) {
   TFile *f = new TFile(filename.c_str(),"recreate");

   for(int iwbc=0; iwbc<nWBC; iwbc++) {
      h_effMap[iwbc][2]->Divide(h_effMap[iwbc][1],h_effMap[iwbc][0],1,1,"B");
      h_effMapWide[iwbc][2]->Divide(h_effMapWide[iwbc][1],h_effMapWide[iwbc][0],1,1,"B");

      for(int i=0; i<3; i++) {
         h_effMap[iwbc][i]->Write();
         h_effMapWide[iwbc][i]->Write();
         if(i<2) {
            h_effFlux[iwbc][i]->Write();
            h_effNHits[iwbc][i]->Write();
            h_effSpill[iwbc][i]->Write();
            h_tpSpill[iwbc][i]->Write();
            for(int iD=0; iD<nD; iD++) {
               h_resSpill[iwbc][iD][i]->Write();
            }
         }
      }
      h_nHitsFlux[iwbc]->Write();
   }
   f->Close();
}

void plotter::makePlots(int startWBC, int finalWBC) {

   // int startWBC = wbc140;
   // int finalWBC = wbc175;

   TCanvas * slide;

   char slideT[256], slideN[256], slideF[256];

   TGraphAsymmErrors *g_effFlux[nWBC];
   TGraphAsymmErrors *g_effNHits[nWBC];
   TGraphAsymmErrors *g_effSpill[nWBC];

   TH1D *effFlux[nWBC][2];
   TH1F *hFlux[nWBC];

   TH1D *effNHits[nWBC][2];
   TH1F *hNHits[nWBC];

   TH1F * hSpill[nWBC];

   for(int iwbc=0; iwbc<nWBC; iwbc++) {
      for(int i=0; i<2; i++) {
         effFlux[iwbc][i] = h_effFlux[iwbc][i]->ProjectionY();
         effFlux[iwbc][i]->Rebin(5);

         effNHits[iwbc][i] = h_effNHits[iwbc][i]->ProjectionY();
      }

      g_effFlux[iwbc] = new TGraphAsymmErrors();
      g_effFlux[iwbc]->Divide(effFlux[iwbc][1],effFlux[iwbc][0],"cl=0.683 b(1,1) mode");
      hFlux[iwbc] = new TH1F(Form("h_effFlux_wbc%d_efficiency",WBCvalue[iwbc]), 
                             "Efficiency vs Flux", 
                             effFlux[iwbc][0]->GetNbinsX(),
                             effFlux[iwbc][0]->GetXaxis()->GetXmin(),
                             effFlux[iwbc][0]->GetXaxis()->GetXmax());
      graphSetting(g_effFlux[iwbc], hFlux[iwbc],
                   WBCcolor[iwbc],WBCstyle[iwbc],
                   TString::Format("Flux"),TString::Format("Efficiency"));

      g_effNHits[iwbc] = new TGraphAsymmErrors();
      g_effNHits[iwbc]->Divide(effNHits[iwbc][1],effNHits[iwbc][0],"cl=0.683 b(1,1) mode");
      hNHits[iwbc] = new TH1F(Form("h_effNHits_wbc%d_efficiency",WBCvalue[iwbc]), 
                             "Efficiency vs nHits", 
                             effNHits[iwbc][0]->GetNbinsX(),
                             effNHits[iwbc][0]->GetXaxis()->GetXmin(),
                             effNHits[iwbc][0]->GetXaxis()->GetXmax());
      graphSetting(g_effNHits[iwbc], hNHits[iwbc],
                   WBCcolor[iwbc],WBCstyle[iwbc],
                   TString::Format("NHits"),TString::Format("Efficiency"));
      
      g_effSpill[iwbc] = new TGraphAsymmErrors();
      g_effSpill[iwbc]->Divide(h_effSpill[iwbc][1],h_effSpill[iwbc][0],"cl=0.683 b(1,1) mode");
      hSpill[iwbc] = new TH1F(Form("h_effSpill_wbc%d_efficiency",WBCvalue[iwbc]), 
                             "Efficiency vs Spill", 
                             h_effSpill[iwbc][0]->GetNbinsX(),
                             h_effSpill[iwbc][0]->GetXaxis()->GetXmin(),
                             h_effSpill[iwbc][0]->GetXaxis()->GetXmax());
      graphSetting(g_effSpill[iwbc], hSpill[iwbc],
                   WBCcolor[iwbc],WBCstyle[iwbc],
                   TString::Format("Spill"),TString::Format("Efficiency"));
   }

   slide = newSlide("eff_vs_flux","");

   TLegend* leg = new TLegend(0.15,0.20,0.38,0.41);
   legendSetting(leg);

   hFlux[0]->GetXaxis()->SetRangeUser(0.,250.);
   hFlux[0]->Draw();
   for(int iwbc=startWBC; iwbc<=finalWBC; iwbc++) {
      if(effFlux[iwbc][0]->GetEntries() == 0) continue;
      g_effFlux[iwbc]->Draw("pe same");
      leg->AddEntry(g_effFlux[iwbc],Form("WBC=%d",WBCvalue[iwbc]),"p");
   }
   leg->Draw();
    slide->SaveAs(Form("%s/eff_vs_flux.%s",_outDir.c_str(),plotExt));

    slide = newSlide("eff_vs_nhits","");
   hNHits[0]->GetXaxis()->SetRange(1,16);
   hNHits[0]->Draw();
   for(int iwbc=startWBC; iwbc<=finalWBC; iwbc++) {
      if(effNHits[iwbc][0]->GetEntries() == 0) continue;
      g_effNHits[iwbc]->Draw("pe same");
   }
   leg->Draw();
   slide->SaveAs(Form("%s/eff_vs_nhits.%s",_outDir.c_str(),plotExt));

   slide = newSlide("eff_vs_spill","");
   hSpill[0]->GetXaxis()->SetNoExponent();
   hSpill[0]->Draw();
   leg = new TLegend(0.15,0.20,0.38,0.41);
   legendSetting(leg);
   for(int iwbc=0; iwbc<nWBC; iwbc++) {
      if(h_effSpill[iwbc][0]->GetEntries() == 0) continue;
      g_effSpill[iwbc]->Draw("pe same");
      leg->AddEntry(g_effFlux[iwbc],Form("WBC=%d",WBCvalue[iwbc]),"p");
   }
   leg->Draw();
   slide->SaveAs(Form("%s/eff_vs_spill.%s",_outDir.c_str(),plotExt));

   //profile nhits vs flux
   TProfile *p_nHitsFlux[nWBC];
   for(int iwbc=startWBC; iwbc<=finalWBC; iwbc++) {
      if(h_nHitsFlux[iwbc]->GetEntries() == 0) continue;
      h_nHitsFlux[iwbc]->SetNdivisions(110,"Y");
      h_nHitsFlux[iwbc]->SetNdivisions(110,"X");
      p_nHitsFlux[iwbc] = h_nHitsFlux[iwbc]->ProfileX();
      slide = newSlide(TString::Format("nhits_vs_flux_wbc%d",WBCvalue[iwbc]),"");
      slide->SetGrid();
      p_nHitsFlux[iwbc]->SetStats(0);
      p_nHitsFlux[iwbc]->Draw();
      slide->SaveAs(Form("%s/nhits_vs_flux_wbc%d.%s",_outDir.c_str(),WBCvalue[iwbc],plotExt));
   }

   for(int iD=0; iD<nD; iD++) {
      slide = newSlide(TString::Format("res%s_vs_spill",D[iD]),"");
      slide->SetGrid();
      TH1F *hRes = new TH1F("hRes", Form("Mean %s residual vs spill",D[iD]),
                            _nSpills, _firstSpill-0.5, _finalSpill+0.5);
      hRes->GetXaxis()->SetNoExponent();
      hRes->SetXTitle("Spill");
      hRes->SetStats(0);
      hRes->SetMinimum(-0.01);
      hRes->SetMaximum( 0.01);
      hRes->Draw();
      for(int iwbc=0; iwbc<nWBC; iwbc++) {
         if(h_resSpill[iwbc][iD][0]->GetEntries() == 0) continue;
         h_resSpill[iwbc][iD][0]->Draw("same");
      }
      slide->SaveAs(Form("%s/res%s_vs_spill.%s",_outDir.c_str(),D[iD],plotExt));
      delete hRes;
   }

   for(int iwbc=0; iwbc<nWBC; iwbc++) {
      if(h_effMap[iwbc][0]->GetEntries() == 0) continue;

      sprintf(slideN,"efficiency_map_WBC%d",WBCvalue[iwbc]);
      sprintf(slideT,"Efficiency Map (WBC=%d)",WBCvalue[iwbc]);
      sprintf(slideF,"%s/%s.%s",_outDir.c_str(),slideN,plotExt);
      slide = new TCanvas(slideN, slideT, 0, 0, 1000, 600);
      slide->SetRightMargin(0.10);
      slide->SetBottomMargin(0.135);
      h_effMap[iwbc][2]->SetStats(0);
      h_effMap[iwbc][2]->SetMinimum(0.7);
      h_effMap[iwbc][2]->SetMaximum(1.);
      h_effMap[iwbc][2]->Draw("colz");
      slide->SaveAs(slideF);

      // TProfile* effMapProfX = h_effMap[iwbc][2]->ProfileX();
      // sprintf(slideN,"effMap_profX_%s_WBC%d",_testBoard.c_str(),WBCvalue[iwbc]);
      // sprintf(slideT,"Efficiency Map X Profile (WBC=%d) %s",WBCvalue[iwbc],_testBoard.c_str());
      // sprintf(slideF,"%s/%s.%s",_outDir.c_str(),slideN,plotExt);
      // slide = new TCanvas(slideN, slideT, 0, 0, 1000, 600);
      // slide->SetRightMargin(0.10);
      // slide->SetBottomMargin(0.135);
      // effMapProfX->SetStats(0);
      // effMapProfX->Draw("e");
      // slide->SaveAs(slideF);

      // TProfile* effMapProfY = h_effMap[iwbc][2]->ProfileY();
      // sprintf(slideN,"effMap_profY_%s_WBC%d",_testBoard.c_str(),WBCvalue[iwbc]);
      // sprintf(slideT,"Efficiency Map Y Profile (WBC=%d) %s",WBCvalue[iwbc],_testBoard.c_str());
      // sprintf(slideF,"%s/%s.%s",_outDir.c_str(),slideN,plotExt);
      // slide = new TCanvas(slideN, slideT, 0, 0, 1000, 600);
      // slide->SetRightMargin(0.10);
      // slide->SetBottomMargin(0.135);
      // effMapProfY->SetStats(0);
      // effMapProfY->Draw("e");
      // slide->SaveAs(slideF);      

      sprintf(slideN,"effTrack_map_WBC%d",WBCvalue[iwbc]);
      sprintf(slideT,"Track Map (WBC=%d)",WBCvalue[iwbc]);
      sprintf(slideF,"%s/%s.%s",_outDir.c_str(),slideN,plotExt);
      slide = new TCanvas(slideN, slideT, 0, 0, 1000, 600);
      slide->SetRightMargin(0.10);
      slide->SetBottomMargin(0.135);
      h_effMap[iwbc][0]->SetStats(0);
      h_effMap[iwbc][0]->Draw("colz");
      slide->SaveAs(slideF);

      sprintf(slideN,"efficiency_mapWide_WBC%d",WBCvalue[iwbc]);
      sprintf(slideT,"Efficiency MapWide (WBC=%d)",WBCvalue[iwbc]);
      sprintf(slideF,"%s/%s.%s",_outDir.c_str(),slideN,plotExt);
      slide = new TCanvas(slideN, slideT, 0, 0, 1000, 600);
      slide->SetRightMargin(0.10);
      slide->SetBottomMargin(0.135);
      h_effMapWide[iwbc][2]->SetStats(0);
      h_effMapWide[iwbc][2]->SetMinimum(0.7);
      h_effMapWide[iwbc][2]->SetMaximum(1.);
      h_effMapWide[iwbc][2]->Draw("colz");
      slide->SaveAs(slideF);

      sprintf(slideN,"effTrack_mapWide_WBC%d",WBCvalue[iwbc]);
      sprintf(slideT,"Track MapWide (WBC=%d)",WBCvalue[iwbc]);
      sprintf(slideF,"%s/%s.%s",_outDir.c_str(),slideN,plotExt);
      slide = new TCanvas(slideN, slideT, 0, 0, 1000, 600);
      slide->SetRightMargin(0.10);
      slide->SetBottomMargin(0.135);
      h_effMapWide[iwbc][0]->SetStats(0);
      h_effMapWide[iwbc][0]->Draw("colz");
      slide->SaveAs(slideF);

   }

   TGraphAsymmErrors *g_tpSpill[_nSpills];
   TH1F *s_tpSpill[_nSpills];
   for(int ispill=_firstSpill; ispill<=_finalSpill; ispill++) {
      int spillBin = 1 + ispill - _firstSpill;

      g_tpSpill[spillBin-1] = new TGraphAsymmErrors();
      s_tpSpill[spillBin-1] = NULL;
      int wbc = 0;
      for(int iwbc=0; iwbc<nWBC; iwbc++) {
         TH1D *tpHits = h_tpSpill[iwbc][1]->ProjectionY("_py",spillBin,spillBin,"e");
         TH1D *tpTrks = h_tpSpill[iwbc][0]->ProjectionY("_py",spillBin,spillBin,"e");
         if(tpTrks->GetSumOfWeights() > 0) {
            g_tpSpill[spillBin-1]->Divide(tpHits, tpTrks, "cl=0.683 b(1,1) mode");
            wbc = WBCvalue[iwbc];
            s_tpSpill[spillBin-1] = new TH1F(Form("h_tpSpill%d_wbc%d",ispill,wbc),
                                         Form("Efficiency vs trigger phase (%d, WBC %d)",ispill,wbc),
                                         nPhases,-0.5,nPhases-0.5);
            graphSetting(g_tpSpill[spillBin-1],s_tpSpill[spillBin-1],
                               WBCcolor[iwbc],WBCstyle[iwbc],
                               TString::Format("Trigger phase"),TString::Format("Efficiency"));
            break;
         }
      }

      if(s_tpSpill[spillBin-1] != NULL) {
         slide = newSlide(TString::Format("eff_vs_tp_%d_wbc%d",ispill,wbc),"");
         s_tpSpill[spillBin-1]->Draw();
         g_tpSpill[spillBin-1]->Draw("pe same");
         slide->SaveAs(Form("%s/eff_vs_tp_%d_wbc%d.%s",_outDir.c_str(),ispill,wbc,plotExt));
      }
   }

}

void plotter::compareSpills(int wbc, vector<int> spillList) {

   int wbcBin=-1;
   cout << "  WBC=" << wbc << endl;
   for(int iwbc=0; iwbc<nWBC; iwbc++) {
      if(wbc == WBCvalue[iwbc]) {
         wbcBin=iwbc;
         break;
      }
   }
   if(wbcBin < 0) {
      cout << "\tCould not find proper WBC bin\n";
      return;
   }

   TCanvas * slide;

   char slideT[256], slideN[256], slideF[256];

   const int nSpills = (int) spillList.size();
   if (nSpills < 2) return; //what are we comparing?

   TGraphAsymmErrors *g_effFlux[nSpills];
   TGraphAsymmErrors *g_effNHits[nSpills];

   TH1D *effFlux[nSpills][2];
   TH1F *hFlux[nSpills];

   TH1D *effNHits[nSpills][2];
   TH1F *hNHits[nSpills];

   const int nFluxBins = 8;
   double fluxBins[nFluxBins+1] = {0., 20., 30., 40., 50., 70., 100., 150., 250.};//, 1000.}; 

   for(int iSpill=0; iSpill<nSpills; iSpill++) {
      int spillBin = 1 + spillList[iSpill] - _firstSpill;
      for(int i=0; i<2; i++) {
         TH1D *htmp = h_effFlux[wbcBin][i]->ProjectionY("_py",spillBin,spillBin,"e");
         effFlux[iSpill][i] = (TH1D*) htmp->Rebin(nFluxBins,Form("%s_%d",h_effFlux[wbcBin][i]->GetName(),spillList[iSpill]),fluxBins);

         effNHits[iSpill][i] = h_effNHits[wbcBin][i]->ProjectionY("_py",spillBin,spillBin,"e");
      }

      g_effFlux[iSpill] = new TGraphAsymmErrors();
      g_effFlux[iSpill]->Divide(effFlux[iSpill][1],effFlux[iSpill][0],"cl=0.683 b(1,1) mode");
      hFlux[iSpill] = new TH1F(Form("h_effFlux_spill%d_efficiency",spillList[iSpill]), 
                             "Efficiency vs Flux", 
                             effFlux[iSpill][0]->GetNbinsX(),
                             effFlux[iSpill][0]->GetXaxis()->GetXmin(),
                             effFlux[iSpill][0]->GetXaxis()->GetXmax());
      graphSetting(g_effFlux[iSpill], hFlux[iSpill],
                   WBCcolor[iSpill],WBCstyle[iSpill],
                   TString::Format("Flux"),TString::Format("Efficiency"));

      g_effNHits[iSpill] = new TGraphAsymmErrors();
      g_effNHits[iSpill]->Divide(effNHits[iSpill][1],effNHits[iSpill][0],"cl=0.683 b(1,1) mode");
      hNHits[iSpill] = new TH1F(Form("h_effNHits_spill%d_efficiency",spillList[iSpill]), 
                             "Efficiency vs NHits", 
                             effNHits[iSpill][0]->GetNbinsX(),
                             effNHits[iSpill][0]->GetXaxis()->GetXmin(),
                             effNHits[iSpill][0]->GetXaxis()->GetXmax());
      graphSetting(g_effNHits[iSpill], hNHits[iSpill],
                   WBCcolor[iSpill],WBCstyle[iSpill],
                   TString::Format("NHits"),TString::Format("Efficiency"));
   }

   slide = newSlide("eff_vs_flux","");

   TLegend* leg = new TLegend(0.15,0.20,0.38,0.41);
   legendSetting(leg);

   hFlux[0]->GetXaxis()->SetRangeUser(0.,250.);
   hFlux[0]->Draw();
   for(int iSpill=0; iSpill<nSpills; iSpill++) {
      if(effFlux[iSpill][0]->GetEntries() == 0) continue;
      g_effFlux[iSpill]->Draw("pe same");
      leg->AddEntry(g_effFlux[iSpill],Form("%d",spillList[iSpill]),"p");
   }
   leg->Draw();
   slide->SaveAs(Form("%s/eff_vs_flux_bySpill.%s",_outDir.c_str(),plotExt));

   slide = newSlide("eff_vs_nhits","");
   hNHits[0]->GetXaxis()->SetRange(1,16);
   hNHits[0]->Draw();
   for(int iSpill=0; iSpill<nSpills; iSpill++) {
      if(effNHits[iSpill][0]->GetEntries() == 0) continue;
      g_effNHits[iSpill]->Draw("pe same");
   }
   leg->Draw();
   slide->SaveAs(Form("%s/eff_vs_nhits_bySpill.%s",_outDir.c_str(),plotExt));

}
