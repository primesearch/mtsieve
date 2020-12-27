/* TwinWorker.h -- (C) Mark Rodenkirch, September 2018

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _TwinWorker_H
#define _TwinWorker_H

#include "TwinApp.h"
#include "../core/Worker.h"

using namespace std;

class TwinWorker : public Worker
{
public:
   TwinWorker(uint32_t myId, App *theApp);

   ~TwinWorker(void) {};

   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   void              CleanUp(void);

   bool              VerifyExternalFactor(bool badFactorIsFatal, uint64_t prime, uint64_t k, uint32_t b, uint32_t n, int32_t c);
   
protected:

private:
   void              RemoveTermsSmallPrime(uint64_t prime, uint64_t k, int32_t c);
   void              RemoveTermsBigPrime(uint64_t prime, uint64_t k, int32_t c);
   bool              VerifyFactor(bool badFactorIsFatal, uint64_t prime, uint64_t k, int32_t c, uint64_t bPowNModP);
   
   void              BuildBaseInverses(void);
   uint32_t          EuclidExtendedGCD(uint32_t a, uint32_t base);
   void              DeterminePrimeTermRange(void);

   TwinApp          *ip_TwinApp;
   
   uint32_t         *ii_BaseInverses;
   uint64_t         *il_PrimeList;
   uint32_t         *ii_InverseList;
   
   uint64_t          il_BpowN;
   uint64_t          il_MinK;
   uint64_t          il_MaxK;
   uint32_t          ii_Base;
   uint32_t          ii_N;
};

#endif

