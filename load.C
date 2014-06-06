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

   //load file and make plots
   ostringstream stream;
   stream << "output/" << board << "_DUT" << dutid << "_" << firstspill << "-" << finalspill << "/" << board << "_DUT" << dutid << "_" << firstspill << "-" << finalspill<< ".root";
   string filename = stream.str();
   ta->loadHistogramsFromFile(filename.c_str());
   ta->makePlots();
}
