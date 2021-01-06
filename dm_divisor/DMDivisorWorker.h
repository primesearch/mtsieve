/* DMDivisorWorker.h -- (C) Mark Rodenkirch, September 2018

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _DMDivisorWorker_H
#define _DMDivisorWorker_H

#include "DMDivisorApp.h"
#include "../core/Worker.h"

using namespace std;

class DMDivisorWorker : public Worker
{
public:
   DMDivisorWorker(uint32_t myId, App *theApp);

   ~DMDivisorWorker(void) {};

   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   void              CleanUp(void);

protected:

private:
   void              RemoveTermsSmallPrime(uint64_t prime, uint64_t k);
   void              RemoveTermsBigPrime(uint64_t prime, uint64_t k);
   
   DMDivisorApp           *ip_DMDivisorApp;
      
   uint64_t          il_MinK;
   uint64_t          il_MaxK;
   uint32_t          ii_N;
};

#endif

