/* SophieGermainWorker.h -- (C) Mark Rodenkirch, July 2020

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _SophieGermainWorker_H
#define _SophieGermainWorker_H

#include "SophieGermainApp.h"
#include "../core/Worker.h"

using namespace std;

class SophieGermainWorker : public Worker
{
public:
   SophieGermainWorker(uint32_t myId, App *theApp);

   ~SophieGermainWorker(void) {};

   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   void              CleanUp(void);

protected:

private:
   void              RemoveTermsSmallPrime(uint64_t prime, uint64_t k, bool firstOfPair);
   void              RemoveTermsLargePrime(uint64_t prime, uint64_t k, bool firstOfPair);
   void              VerifyFactor(uint64_t prime, uint64_t k, bool firstOfPair);
   
   SophieGermainApp *ip_SophieGermainApp;
      
   uint64_t          il_BpowN;
   uint64_t          il_MinK;
   uint64_t          il_MaxK;
   uint32_t          ii_N;
};

#endif

