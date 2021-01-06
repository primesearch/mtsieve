/* GFNDivisorWorker.h -- (C) Mark Rodenkirch, July 2017

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _GFNDivisorWorker_H
#define _GFNDivisorWorker_H

#include "GFNDivisorApp.h"
#include "../core/Worker.h"

using namespace std;

class GFNDivisorWorker : public Worker
{
public:
   GFNDivisorWorker(uint32_t myId, App *theApp);

   ~GFNDivisorWorker(void) {};

   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   void              CleanUp(void);

protected:

private:
   void              RemoveTermsSmallPrime(uint64_t k, uint32_t n, uint64_t prime);
   void              RemoveTermsBigPrime(uint64_t k, uint32_t n, uint64_t prime);
   void              VerifyFactor(uint64_t k, uint32_t n, uint64_t prime);

   GFNDivisorApp    *ip_GFNDivisorApp;
   
   uint64_t          il_KCount;
   uint32_t          ii_NCount;
   uint64_t          il_MinK;
   uint64_t          il_MaxK;
   uint32_t          ii_MinN;
   uint32_t          ii_MaxN;
};

#endif

