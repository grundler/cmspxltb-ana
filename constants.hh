//Directory where trees are found (in heirarchy <board>/<treeSubDir>)
static const char subdir[128] = "/afs/cern.ch/work/g/grundler/public/testbeam/cmspxltb-submission";

//Place to put output
static const char outDir[64] = "output";

//extension for plots, probably want pdf or png to go in LaTeX
static const char plotExt[8] = "pdf";

static const float fnalClock = 53.104; //MHz
static const float tbClock   = 40.; //MHz

static const int nPhases = 8; //trigger phases

enum WBC{wbc99,
         wbc140,
         wbc159,
         wbc175,
         wbc200,
         wbc225,
         wbc255,
         nWBC
};

static const int WBCvalue[nWBC] = {99, 140, 159, 175, 200, 225, 255};
static const int WBCcolor[nWBC] = {5, 1, 2, 4, 3, 6, 7};
static const int WBCstyle[nWBC] = {24, 20, 21, 22, 23, 25, 26};
