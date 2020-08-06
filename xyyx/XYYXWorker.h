/* XYYXWorker.h -- (C) Mark Rodenkirch, June 2014

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _XYYXWorker_H
#define _XYYXWorker_H

#include "XYYXApp.h"
#include "../core/Worker.h"

using namespace std;

#define MAX_POWERS   50

class XYYXWorker : public Worker
{
public:
   XYYXWorker(uint32_t myId, App *theApp);

   ~XYYXWorker() {};

   void           TestMegaPrimeChunk(void);
   void           TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   void           CleanUp(void);

private:         
   void           FreeTerms(void);
   
   void           TestPrimeChunkFPU(uint64_t &largestPrimeTested, uint64_t &primesTested);   
   void           BuildFpuPowmodTables(base_t *terms, uint32_t minBase, uint32_t maxBase, uint64_t *ps);
   void           CheckForFpuFactors(uint64_t *ps);
   
   void           TestPrimeChunkAVX(uint64_t &largestPrimeTested, uint64_t &primesTested);
   void           BuildAvxXYRemainders(uint64_t *ps, double *dps, double *reciprocals);
   void           CheckAvxXYRemainders(uint64_t *ps, double *dps, double *reciprocals);
   void           BuildAvxListOfPowers(uint32_t base, double *dps, double *reciprocals, uint32_t count);
   void           CheckAvxResult(uint32_t x, uint32_t y, uint64_t *ps, double *dps, power_t *pair);
   
   void           VerifyFactor(uint64_t p, uint32_t x, uint32_t y, int32_t c);

   XYYXApp       *ip_XYYXApp;

   uint32_t       ii_MinX;
   uint32_t       ii_MaxX;
   uint32_t       ii_MinY;
   uint32_t       ii_MaxY;
   uint32_t       ii_XCount;
   uint32_t       ii_YCount;
   uint32_t       ii_MaxPowerDiff;

   bool           ib_IsPlus;
   bool           ib_IsMinus;

   uint64_t       il_NextTermsBuild;

   base_t        *ip_xyTerms;
   base_t        *ip_yxTerms;
   
   uint64_t     **ip_FpuPowers;
   
   double        *ip_AvxPowers[MAX_POWERS+1];
};

#endif

