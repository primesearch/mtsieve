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

using namespace std;

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
   
   void              CheckAVXResult(uint64_t *ps, double *dps, uint32_t theN);
   void              VerifyAVXFactor(uint64_t p, uint32_t theN, int32_t theC);
   
   uint32_t         *ii_Primes;
   double           *id_Primes;
};

#endif

