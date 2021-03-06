#ifndef THaCodaData_h
#define THaCodaData_h

/////////////////////////////////////////////////////////////////////
//
//   THaCodaData
//   Abstract class of CODA data.
//
//   THaCodaData is an abstract class of CODA data.  Derived
//   classes will be typically either a CODA file (a disk
//   file) or a connection to the ET system.  Public methods
//   allow to open (i.e. set up), read, write, etc.
//
//   author Robert Michaels (rom@jlab.org)
//
/////////////////////////////////////////////////////////////////////


#include <Rtypes.h>
#include "TString.h"
#define CODA_ERROR -1     // Generic error return code
#define CODA_OK  0        // Means return is ok.
#define MAXEVLEN 200000    // Maximum size of events
#define CODA_VERBOSE 1    // Errors explained verbosely (recommended)
#define CODA_DEBUG  1     // Lots of printout (recommend to set = 0)

using namespace std;

class THaCodaData {

public:

   THaCodaData();
   virtual ~THaCodaData();
   virtual int codaOpen(TString filename)=0;
   virtual int codaOpen(TString filename, TString session) {return CODA_OK;};
   virtual int codaOpen(TString filename, TString session, int mode) {return CODA_OK;};
   virtual int codaClose()=0;
   virtual int codaRead()=0; 
   virtual unsigned *getEvBuffer() { return evbuffer; };     
   virtual int getBuffSize() const { return MAXEVLEN; };

private:

   THaCodaData(const THaCodaData &fn);
   THaCodaData& operator=(const THaCodaData &fn);

protected:

   TString filename;
   unsigned *evbuffer;                    // Raw data     

#ifndef STANDALONE
   ClassDef(THaCodaData,0) // Base class of CODA data (file, ET conn, etc)
#endif

};

#endif







