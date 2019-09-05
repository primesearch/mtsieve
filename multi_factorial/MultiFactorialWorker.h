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

protected:
   MultiFactorialApp      *ip_MultiFactorialApp;

private:
   double           *id_Terms;
   
   void              ExtractFactors(uint32_t start, uint64_t p);
   void              CheckAVXResult(uint64_t *ps, double *dps, uint32_t theN);
   void              VerifyAVXFactor(uint64_t p, uint32_t theN, int32_t theC);
};

#endif

