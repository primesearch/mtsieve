/* MultiFactorialWorker.h -- (C) Mark Rodenkirch, October 2012

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _MultiFactorialWorker_H
#define _MultiFactorialWorker_H

#include "MultiFactorialApp.h"
#include "../core/Worker.h"

using namespace std;

class MultiFactorialWorker : public Worker
{
public:
   MultiFactorialWorker(uint32_t myId, App *theApp);

   ~MultiFactorialWorker() {};

   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   void              CleanUp(void);
   
   bool              VerifyFactor(bool badFactorIsFatal, uint64_t p, uint32_t theN, int32_t theC);

protected:
   MultiFactorialApp      *ip_MultiFactorialApp;
   
   uint32_t          ii_MinN;
   uint32_t          ii_MaxN;
   uint32_t          ii_MultiFactorial;

private:
   void              TestFactorial(void);
   void              TestMultiFactorial(void);
};

#endif

