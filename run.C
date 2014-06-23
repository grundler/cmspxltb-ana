{
   //Settings
   int dutid = 3;
   std::string board("PixelTestBoard1");
   int firstspill = 176620;
   int finalspill = 176997;
   int mapalgo = 0;

   //compile macros
   gSystem->CompileMacro("treeCorrelator.cc","k");
   gSystem->CompileMacro("utils.cc","k");
   gSystem->CompileMacro("tbAna.cc","k");

   //create instance
   tbAna *ta = new tbAna(dutid,board,firstspill,finalspill,mapalgo);

   //disable cuts
   //ta->useCorrectTriggerPhase(false);
   //ta->useFiducial(false);
   //ta->useSlope(false);
   //ta->setMaxFluxRatio(9999999.);

   //run analysis and make plots
   ta->analyze();
   ta->makePlots();
}
