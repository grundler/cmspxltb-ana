{
   gSystem->CompileMacro("treeCorrelator.cc","k");
   gSystem->CompileMacro("utils.cc","k");
   gSystem->CompileMacro("tbAna.cc","k");
   //tbAna *ta = new tbAna(3,"PixelTestBoard1",176800,176820, 0);
   tbAna *ta = new tbAna(3,"PixelTestBoard1",176620,176997, 0);
   ta->analyze();
   ta->makePlots();
}
