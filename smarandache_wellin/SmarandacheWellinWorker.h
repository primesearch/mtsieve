/* SmarandacheWellinWorker.h -- (C) Mark Rodenkirch, March 2025

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _SmarandacheWellinWorker_H
#define _SmarandacheWellinWorker_H

#include "SmarandacheWellinApp.h"
#include "../core/Worker.h"

using namespace std;

class SmarandacheWellinWorker : public Worker
{
public:
   SmarandacheWellinWorker(uint32_t myId, App *theApp);

   ~SmarandacheWellinWorker() {};

   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   void              CleanUp(void);

protected:
   void              NotifyPrimeListAllocated(uint32_t primesInList) {}
   
   SmarandacheWellinApp   *ip_SmarandacheWellinApp;
   
   uint32_t          ii_MinN;
   uint32_t          ii_MaxN;
   
   uint32_t         *ip_Primes;
   uint32_t          ii_NumberOfPrimes;
};

#endif

