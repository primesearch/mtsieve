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
#include "../core/MpArith.h"
#include "../core/MpArithVector.h"

using namespace std;

class CullenWoodallWorker : public Worker
{
public:
   CullenWoodallWorker(uint32_t myId, App *theApp);

   ~CullenWoodallWorker() {};

   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   void              CleanUp(void);
   
protected:
   void              NotifyPrimeListAllocated(uint32_t primesInList) {}
   
private:   
   void              TestSmallPrimes(uint64_t *ps);
   void              TestLargePrimes(uint64_t *ps, MpArithVec mp);

#ifdef USE_X86
   void              TestLargePrimesFPU(uint64_t *ps);
   void              TestPrimesAVX(uint64_t *ps);
   void              CheckAVXResult(uint32_t theN, uint64_t *ps, double *dps);
#endif

   void              BuildListOfPowers(MpRes a, uint64_t p, MpArith mp, uint32_t count, MpRes *powers);

   CullenWoodallApp *ip_CullenWoodallApp;
   
   uint32_t          ii_Base;
   uint32_t          ii_MinN;
   uint32_t          ii_MaxN;
   
   uint32_t          ii_MaxTermCount;
   uint32_t         *ii_Terms;
   uint64_t          il_NextTermsBuild;
};

#endif

