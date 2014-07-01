cmspxltb-ana
============

Analysis of data from CMS pixel high-rate beam tests

There are three classes:
	  tbAna - main class, runs the analysis
	  		Constructor takes ID of DUT (int), board name (string) used for directory identification, and first and last spill numbers (int), and optional argument to define algorithm for tree correlator to use in matching testboard and QIE data.

	  treeCorrelator - used to correlate information from QIE summary and trigger phase trees to the track tree
	  				 Called by tbAna

	  plotter - used for storing histograms and creating plots			 

constants.hh contains a const static char variable, subdir, indicating the directory where root files are at.
			 track trees assumed to be at "subdir/<TestBoardName>/histograms/<spill>-tracks.root"
			 QIE summary assumed to be at "subdir/<TestBoardName>/qie/summary_<spill>.root"
			 tp trees assumed to be at "subdir/<TestBoardName>/timestamps/<spill>_tp.root"

To simply compile, in root:
   .x compile.C

Edit run.C to change settings (DUT ID, boardname, spill range) and add any cut you wish, e.g.  ta->analyze("Chi2<25.");

To run (includes compilation), in root:
   .x run.C

One can also load histograms from previously run analysis. See load.C