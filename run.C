{
   gSystem->CompileMacro("treeCorrelator.cc","k");
   gSystem->CompileMacro("tbAna.cc","k");
   tbAna *ta = new tbAna(3,"PixelTestBoard1",176800,176820);
   ta->analyze();
}
