/* CullenWoodallWorker.h -- (C) Mark Rodenkirch, May 2018
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _CullenWoodallWorker_H
#define _CullenWoodallWorker_H

#include "CullenWoodallApp.h"
#include "../core/Worker.h"

using namespace std;

class CullenWoodallWorker : public Worker
{
public:
   CullenWoodallWorker(uint32_t myId, App *theApp);

   ~CullenWoodallWorker() {};

   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   void              CleanUp(void);

private:   
   void              TestSmallPrimesFPU(uint64_t *ps);
   void              TestLargePrimesFPU(uint64_t *ps);
   void              TestPrimesAVX(uint64_t *ps);
   void              CheckAVXResult(uint32_t theN, uint64_t *ps, double *dps);
   
   void              BuildListOfPowers(uint64_t a, uint64_t p, uint32_t count, uint64_t *powers);
   void              VerifyFactor(uint64_t p, uint32_t n, int32_t c);
   uint64_t          ComputeMultiplicativeInverse(uint64_t a, uint64_t p);

   CullenWoodallApp *ip_CullenWoodallApp;
   
   uint32_t          ii_Base;
   uint32_t          ii_MinN;
   uint32_t          ii_MaxN;
   
   uint32_t          ii_MaxTermCount;
   uint32_t         *ii_Terms;
   uint64_t          il_NextTermsBuild;
};

#endif

