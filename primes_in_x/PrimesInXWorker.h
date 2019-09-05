/* PrimesInXWorker.h -- (C) Mark Rodenkirch, October 2012

   This is the header for the class sets up the call to the PrimesInXWorker GPU function
   and parses the output from the GPU to determine if we have a PrimesInXWorker prime.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _PrimesInXWorker_H
#define _PrimesInXWorker_H

#include "PrimesInXApp.h"
#include "../core/Worker.h"

using namespace std;

class PrimesInXWorker : public Worker
{
public:
   PrimesInXWorker(uint32_t myId, App *theApp);

   ~PrimesInXWorker() {};

   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   void              CleanUp(void);

private:
   void              ExtractFactors(uint64_t p);

   PrimesInXApp     *ip_PrimesInXApp;

   uint32_t         *ii_e3Terms;                        // The string as an array of terms < 1000
   uint32_t         *ii_e6Terms;                        // The string as an array of terms < 1000000 
   uint32_t         *ii_e9Terms;                        // The string as an array of trems < 1000000000
};

#endif

