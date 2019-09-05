/* FixedBNCWorker.h -- (C) Mark Rodenkirch, July 2017

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _FixedBNCWorker_H
#define _FixedBNCWorker_H

#include "FixedBNCApp.h"
#include "../core/Worker.h"

using namespace std;

class FixedBNCWorker : public Worker
{
public:
   FixedBNCWorker(uint32_t myId, App *theApp);

   ~FixedBNCWorker(void) {};

   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   void              CleanUp(void);

protected:

private:
   void              RemoveTermsSmallPrime(uint64_t prime, uint64_t k);
   void              RemoveTermsBigPrime(uint64_t prime, uint64_t k);
   void              VerifyFactor(uint64_t prime, uint64_t k, uint64_t bPowNModP);

   void              BuildBaseInverses(void);
   uint32_t          EuclidExtendedGCD(uint32_t a, uint32_t base);
   void              DeterminePrimeTermRange(void);

   FixedBNCApp      *ip_FixedBNCApp;
   
   uint32_t         *ii_BaseInverses;
   uint64_t         *il_PrimeList;
   uint32_t         *ii_InverseList;
   
   uint64_t          il_BpowN;
   uint64_t          il_MinK;
   uint64_t          il_MaxK;
   uint32_t          ii_Base;
   uint32_t          ii_N;
   int32_t           ii_C;
};

#endif

