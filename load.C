{
   //Settings
   int dutid = 3;
   std::string board("PixelTestBoard1");
   int firstspill = 179410;
   int finalspill = 179620;
   std::string suffix("");

   //compile macros
   gSystem->CompileMacro("plotter.cc","k");

   //load file and make plots
   ostringstream stream;
   stream << "output/" << board << "_DUT" << dutid << "_" << firstspill << "-" << finalspill;
   if(!suffix.empty())
      stream << "_" << suffix;
   string outdir = stream.str();
   stream << "/" << board << "_DUT" << dutid << "_" << firstspill << "-" << finalspill<< ".root";
   string filename = stream.str();

   plotter *p = new plotter(firstspill,finalspill,outdir);
   p->loadHistogramsFromFile(filename.c_str());
   p->makePlots();
   // std::vector<int> v;
   // v.push_back(179423);
   // v.push_back(179425);
   // v.push_back(179430);
   // p->compareSpills(159,v);
}
