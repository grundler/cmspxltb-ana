static const char subdir[128] = "/afs/cern.ch/work/g/grundler/public/testbeam/cmspxltb-submission";

static const float fnalClock = 53.104; //MHz
static const float tbClock   = 40.; //MHz

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
