/* PrimorialWorker.h -- (C) Mark Rodenkirch, July 2018

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _PrimorialWorker_H
#define _PrimorialWorker_H

#include "PrimorialApp.h"
#include "../core/Worker.h"
#include "../core/MpArithVector.h"

// The first prime gap over 300 is at 2e9.  Unlikely anyone will ever search that far
// in the foreseeable future.
#define MAX_GAPS     300

class PrimorialWorker : public Worker
{
public:
   PrimorialWorker(uint32_t myId, App *theApp);

   ~PrimorialWorker() {};

   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   void              CleanUp(void);

protected:
   PrimorialApp     *ip_PrimorialApp;

private:
   void              ExtractFactors(uint64_t p);
   
#ifdef USE_X86
   void              CheckAVXResult(uint64_t *ps, double *dps, uint32_t theN);
   void              VerifyAVXFactor(uint64_t p, uint32_t theN, int32_t theC);
#endif

   uint32_t          ii_MinPrimorial;
   uint32_t          ii_MaxPrimorial;
   
   double           *id_PrimorialPrimes;
   uint32_t         *ip_PrimorialPrimes;
   uint32_t          ii_NumberOfPrimorialPrimes;
   
   uint16_t         *ip_PrimorialPrimeGaps;
   uint16_t          ii_BiggestGap;
   
   MpResVec          ip_ResGaps[MAX_GAPS];
};

#endif

