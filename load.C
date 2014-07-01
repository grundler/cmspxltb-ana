{
   //Settings
   int dutid = 3;
   std::string board("PixelTestBoard1");
   int firstspill = 179410;
   int finalspill = 179620;

   //compile macros
   gSystem->CompileMacro("plotter.cc","k");

   //load file and make plots
   ostringstream stream;
   stream << "output/" << board << "_DUT" << dutid << "_" << firstspill << "-" << finalspill;
   string outdir = stream.str();
   stream << "/" << board << "_DUT" << dutid << "_" << firstspill << "-" << finalspill<< ".root";
   string filename = stream.str();

   plotter *p = new plotter(firstspill,finalspill,outdir);
   p->loadHistogramsFromFile(filename.c_str());
   p->makePlots();
}
