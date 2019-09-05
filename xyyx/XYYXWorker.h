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

// This is for building a list of even powers for x and y.
#define MAX_POWERS   50

class XYYXWorker : public Worker
{
public:
   XYYXWorker(uint32_t myId, App *theApp);

   ~XYYXWorker() {};

   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   void              CleanUp(void);

private:
   void              TestPrimeChunkFPU(uint64_t &largestPrimeTested, uint64_t &primesTested);   
   void              BuildXYRemainders(uint64_t p);
   void              CheckXYRemainders(uint64_t p);
   void              BuildListOfPowers(uint64_t a, uint64_t p, uint32_t count, uint64_t *powers);
   
   void              TestPrimeChunkAVX(uint64_t &largestPrimeTested, uint64_t &primesTested);
   void              BuildAvxXYRemainders(uint64_t *ps, double *dps, double *reciprocals);
   void              CheckAvxXYRemainders(uint64_t *ps, double *dps, double *reciprocals);
   void              BuildAvxListOfPowers(uint32_t base, double *dps, double *reciprocals, uint32_t count);
   void              CheckAvxResult(uint32_t x, uint32_t y, uint64_t *ps, double *dps);
   
   void              VerifyFactorFPU(uint64_t p, uint32_t x, uint32_t y, int32_t c);
   void              VerifyFactorAVX(uint64_t p, uint32_t x, uint32_t y, int32_t c);

   XYYXApp          *ip_XYYXApp;
   
   uint32_t          ii_MinX;
   uint32_t          ii_MaxX;
   uint32_t          ii_MinY;
   uint32_t          ii_MaxY;
   uint32_t          ii_XCount;
   uint32_t          ii_YCount;
   
   bool              ib_IsPlus;
   bool              ib_IsMinus;
   
   uint32_t          ii_TermsElements;
   
   uint64_t          il_NextTermsBuild;
   
   // The ip_xyTerms list is a flat array.  It starts with the value for x, and is followed by each y
   // for that x that does not have a factor.  The list of y ends with a 0.  Then the next x,
   // and so on.  The list is done when it ends with two zeroes.
   uint32_t         *ip_xyTerms;
   
   // The ip_yxTerms list is a flat array.  It starts with the value for y, and is followed by each x
   // for that y that does not have a factor.  The list of x ends with a 0.  Then the next y,
   // and so on.  The list is done when it ends with two zeroes.
   uint32_t         *ip_yxTerms;
   
   uint64_t       **ip_XYRemainders;
   double        ***ip_AvxXYRemainders;
   
   double          *ip_Powers[MAX_POWERS+1];
};

#endif

