/* XYYXWorker.cpp -- (C) Mark Rodenkirch, September 2012

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <time.h>

#include "XYYXWorker.h"
#include "../x86_asm/fpu-asm-x86.h"
#include "../x86_asm/avx-asm-x86.h"

#define X_INDEX(x)  ((x) - ii_MinX)
#define Y_INDEX(y)  ((y) - ii_MinY)

XYYXWorker::XYYXWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{  
   bases_t  bases;
   uint32_t idx;
   
   ip_XYYXApp = (XYYXApp *) theApp;
   
   ib_IsPlus = ip_XYYXApp->IsPlus();
   ib_IsMinus = ip_XYYXApp->IsMinus();
   
   ii_MinX = ip_XYYXApp->GetMinX();
   ii_MaxX = ip_XYYXApp->GetMaxX();
   ii_MinY = ip_XYYXApp->GetMinY();
   ii_MaxY = ip_XYYXApp->GetMaxY();

   ii_XCount = ip_XYYXApp->GetXCount();
   ii_YCount = ip_XYYXApp->GetYCount();
   
   ip_xyTerms = ip_yxTerms = 0;
   
   ip_XYYXApp->GetTerms(CpuSupportsAvx(), AVX_ARRAY_SIZE, &bases);

   ip_xyTerms = bases.xPowY;
   ip_yxTerms = bases.yPowX;
   
   il_NextTermsBuild = 100;
   ib_Initialized = true;
   ii_MaxPowerDiff = ip_XYYXApp->GetMaxPowerDiff();

   ip_FpuPowers = (uint64_t **) malloc(ii_MaxPowerDiff * sizeof(uint64_t *));
      
   for (idx=0; idx<ii_MaxPowerDiff; idx++)
      ip_FpuPowers[idx] = (uint64_t *) xmalloc(4 * sizeof(uint64_t));
      
   if (ip_XYYXApp->UseAvxIfAvailable() && CpuSupportsAvx())
   {
      for (uint32_t i=0; i<=MAX_POWERS; i++)
         ip_AvxPowers[i] = (double *) xmalloc(AVX_ARRAY_SIZE * sizeof(double));
      
      // The AVX logic does not handle when gcd(x, prime) > 1 or gcd(y, prime) > 1
      if (ii_MaxX < ii_MaxY)
         SetMiniChunkRange(ii_MaxY + 1, PMAX_MAX_52BIT, AVX_ARRAY_SIZE);
      else
         SetMiniChunkRange(ii_MaxX + 1, PMAX_MAX_52BIT, AVX_ARRAY_SIZE);
   }
}

void  XYYXWorker::CleanUp(void)
{
   FreeTerms();

   if (ip_XYYXApp->UseAvxIfAvailable() && CpuSupportsAvx())
   {
      for (uint32_t i=0; i<=MAX_POWERS; i++)
         xfree(ip_AvxPowers[i]);
   }
   
   for (uint32_t idx=0; idx<ii_MaxPowerDiff; idx++)
      xfree(ip_FpuPowers[idx]);
   
   xfree(ip_FpuPowers);
}

void  XYYXWorker::FreeTerms(void)
{
   uint32_t x, y, idx;
   
   if (ip_xyTerms == 0)
      return;
   
   for (x=ii_MinX; x<=ii_MaxX; x++)
   {
      if (ip_xyTerms[x-ii_MinX].powerCount == 0)
         continue;
      
      if (ip_XYYXApp->UseAvxIfAvailable() && CpuSupportsAvx())
         for (idx=0; idx<ip_xyTerms[x-ii_MinX].powerCount; idx++)  
            xfree(ip_xyTerms[x-ii_MinX].powers[idx].avxRemainders);
      
      xfree(ip_xyTerms[x-ii_MinX].neededPowers);
      xfree(ip_xyTerms[x-ii_MinX].powers);
   }
   
   xfree(ip_xyTerms);
   
   for (y=ii_MinY; y<=ii_MaxY; y++)
   {
      if (ip_yxTerms[y-ii_MinY].powerCount == 0)
         continue;
            
      xfree(ip_yxTerms[y-ii_MinY].neededPowers);
      xfree(ip_yxTerms[y-ii_MinY].powers);
   }
   
   xfree(ip_yxTerms);
}

void  XYYXWorker::TestMegaPrimeChunk(void)
{
   uint64_t ps[4];
   uint64_t maxPrime = ip_App->GetMaxPrime();
   
   vector<uint64_t>::iterator it = iv_Primes.begin();
   
   while (it != iv_Primes.end())
   {
      ps[0] = *it;
      it++;
      
      ps[1] = *it;
      it++;
      
      ps[2] = *it;
      it++;
      
      ps[3] = *it;
      it++;
      
      // Every once in a while rebuild the term lists as it will have fewer entries
      // which will speed up testing for the next range of p.
      if (ps[0] > il_NextTermsBuild)
      {
         bases_t bases;
         
         FreeTerms();
         
         ip_XYYXApp->GetTerms(CpuSupportsAvx(), AVX_ARRAY_SIZE, &bases);
               
         ip_xyTerms = bases.xPowY;
         ip_yxTerms = bases.yPowX;
      
         il_NextTermsBuild = (ps[3] << 1);
      }
      
      // Compute x^y for y
      BuildFpuPowmodTables(ip_xyTerms, ii_MinX, ii_MaxX, ps);

      // Compute y^x for x
      BuildFpuPowmodTables(ip_yxTerms, ii_MinY, ii_MaxY, ps);

      CheckForFpuFactors(ps);

      SetLargestPrimeTested(ps[3], 4);
      
      if (ps[3] > maxPrime)
         break;
   }
}

// Pre-compute b^2 thru b^64.  then b^(64*n) for n up to maxPowerDiff.
// Find the other powers that are needed for the base and populate them.  This
// way we only need a mulmod to compute all b^n.
void  XYYXWorker::BuildFpuPowmodTables(base_t *terms, uint32_t minBase, uint32_t maxBase, uint64_t *ps) 
{  
   uint32_t  base;
   base_t   *bPtr;
   uint32_t  idx, n, prevN;
   uint64_t  as[4], bs[4];
   
   for (base=minBase; base<=maxBase; base++)
   {
      bPtr = &terms[base-minBase];
      
      if (bPtr->powerCount == 0)
         continue;

      base = terms->base;

      terms->powers[0].fpuRemainders[0] = base;
      terms->powers[0].fpuRemainders[1] = base;
      terms->powers[0].fpuRemainders[2] = base;
      terms->powers[0].fpuRemainders[3] = base;
      
      fpu_powmod_4b_1n_4p(terms->powers[0].fpuRemainders, terms->powers[0].exponent, ps);

      fpu_push_1divp(ps[3]);
      fpu_push_1divp(ps[2]);
      fpu_push_1divp(ps[1]);
      fpu_push_1divp(ps[0]);
   
      as[0] = as[1] = as[2] = as[3] = base;
      bs[0] = bs[1] = bs[2] = bs[3] = base;

      // Compute b^2*idx thru b^32*idx
      for (idx=1; idx<=32; idx++)
      {
         fpu_mulmod_4a_4b_4p(as, bs, ps);
         
         ip_FpuPowers[idx][0] = as[0];
         ip_FpuPowers[idx][1] = as[1];
         ip_FpuPowers[idx][2] = as[2];
         ip_FpuPowers[idx][3] = as[3];
      }   
      
      bs[0] = ip_FpuPowers[32][0];
      bs[1] = ip_FpuPowers[32][1];
      bs[2] = ip_FpuPowers[32][2];
      bs[3] = ip_FpuPowers[32][3];

      // Compute b^(32*2*idx) thru b^(32*2*idx)
      // This is effectively 2^64, 2^128, ... thru 2^1280
      for (idx=2; idx<=20; idx++)
      {
         fpu_mulmod_4a_4b_4p(as, bs, ps);
         
         ip_FpuPowers[idx*32][0] = as[0];
         ip_FpuPowers[idx*32][1] = as[1];
         ip_FpuPowers[idx*32][2] = as[2];
         ip_FpuPowers[idx*32][3] = as[3];    
      }

      prevN = 0;
      
      // Fill in the table for the remaining powers that we need.
      // If we need b^78, then this is computed as b^64 * b^14, but
      // if we do not need b^79, then we won't compute it.
      for (n=32; n<=terms->maxPowerDiff/2; n++)
      {
         // If we already computed this power, skip it.
         if ((n % 32) == 0)
         {
            bs[0] = ip_FpuPowers[n][0];
            bs[1] = ip_FpuPowers[n][1];
            bs[2] = ip_FpuPowers[n][2];
            bs[3] = ip_FpuPowers[n][3];
            
            prevN = n;
            continue;
         }

         // If we don't need this power, skip it.
         if (terms->neededPowers[n] == 0)
            continue;
         
         ip_FpuPowers[n][0] = ip_FpuPowers[n - prevN][0];
         ip_FpuPowers[n][1] = ip_FpuPowers[n - prevN][1];
         ip_FpuPowers[n][2] = ip_FpuPowers[n - prevN][2];
         ip_FpuPowers[n][3] = ip_FpuPowers[n - prevN][3];

         fpu_mulmod_4a_4b_4p(ip_FpuPowers[n], bs, ps);
         
         prevN = n;
      }
      
      terms->powers[0].fpuRemainders[0] = base;
      terms->powers[0].fpuRemainders[1] = base;
      terms->powers[0].fpuRemainders[2] = base;
      terms->powers[0].fpuRemainders[3] = base;
      
      fpu_powmod_4b_1n_4p(terms->powers[0].fpuRemainders, terms->powers[0].exponent, ps);
      
      for (n=1; n<terms->powerCount; n++)
      {
         terms->powers[n].fpuRemainders[0] = ip_FpuPowers[terms->powers[n].diffPowerIdx][0];
         terms->powers[n].fpuRemainders[1] = ip_FpuPowers[terms->powers[n].diffPowerIdx][1];
         terms->powers[n].fpuRemainders[2] = ip_FpuPowers[terms->powers[n].diffPowerIdx][2];
         terms->powers[n].fpuRemainders[3] = ip_FpuPowers[terms->powers[n].diffPowerIdx][3];
            
         fpu_mulmod_4a_4b_4p(terms->powers[n].fpuRemainders, terms->powers[n-1].fpuRemainders, ps);
      }

      fpu_pop();
      fpu_pop();
      fpu_pop();
      fpu_pop();
      
      terms++;
   }
}

void  XYYXWorker::CheckForFpuFactors(uint64_t *ps)
{
   uint32_t x, y, xyIdx;
   base_t   *xyTerms;
   power_t  *yPowX;
   
   xyTerms = ip_xyTerms;
   
   while (xyTerms->base <= ii_MaxX)
   {
      if (xyTerms->powerCount == 0)
      {
         xyTerms++;
         continue;
      }
      
      x = xyTerms->base;
      
      for (xyIdx=0; xyIdx<=xyTerms->powerCount; xyIdx++)
      {
         y = xyTerms->powers[xyIdx].exponent;
         
         yPowX = (power_t *) xyTerms->powers[xyIdx].pairedPower;
         
         for (uint32_t i=0; i<4; i++)
         {
            if (ib_IsMinus && xyTerms->powers[xyIdx].fpuRemainders[i] == yPowX->fpuRemainders[i])
            {                  
               if (ip_XYYXApp->ReportFactor(ps[i], x, y, -1))
                  VerifyFactor(ps[i], x, y, -1);
            }
            
            if (ib_IsPlus && xyTerms->powers[xyIdx].fpuRemainders[i] == ps[i] - yPowX->fpuRemainders[i])           
            {                  
               if (ip_XYYXApp->ReportFactor(ps[i], x, y, +1))
                  VerifyFactor(ps[i], x, y, +1);
            }
         }
      }
      
      xyTerms++;
   }
}

void  XYYXWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   double __attribute__((aligned(32))) dps[AVX_ARRAY_SIZE];
   double __attribute__((aligned(32))) reciprocals[AVX_ARRAY_SIZE];
   
   // Every once in a while rebuild the term lists as it will have fewer entries
   // which will speed up testing for the next range of p.
   if (miniPrimeChunk[0] > il_NextTermsBuild)
   {
      bases_t bases;
      
      FreeTerms();
      
      ip_XYYXApp->GetTerms(CpuSupportsAvx(), AVX_ARRAY_SIZE, &bases);
            
      ip_xyTerms = bases.xPowY;
      ip_yxTerms = bases.yPowX;
      
      il_NextTermsBuild = (miniPrimeChunk[AVX_ARRAY_SIZE-1] << 1);
   }
   
   // compute the inverse of b (mod p)
   for (int i=0; i<AVX_ARRAY_SIZE; i++)
      dps[i] = (double) miniPrimeChunk[i];
      
   avx_compute_reciprocal(dps, reciprocals);

   BuildAvxXYRemainders(miniPrimeChunk, dps, reciprocals);

   CheckAvxXYRemainders(miniPrimeChunk, dps, reciprocals);
}

// Build a table of x^y mod p for all remaining terms
void  XYYXWorker::BuildAvxXYRemainders(uint64_t *ps, double *dps, double *reciprocals) 
{  
   double    __attribute__((aligned(32))) xPowY[AVX_ARRAY_SIZE];
   uint32_t  x, y, prevY;
   uint32_t  yIndex, powIndex;
   uint32_t  maxPowers;
   base_t   *xyPtr;
   
   // If the range of y is only 50, then we only want to generate
   // up to y^50 instead of y^100  (assuming MAX_POWERS = 50).
   if (ii_YCount < MAX_POWERS * 2)
      maxPowers = (ii_YCount / 2);
   else
      maxPowers = MAX_POWERS;

   for (x=ii_MinX; x<=ii_MaxX; x++)
   {
      xyPtr = &ip_xyTerms[x-ii_MinX];
      
      if (xyPtr->powerCount == 0)
         continue;

      BuildAvxListOfPowers(x, dps, reciprocals, maxPowers);
      
      for (int i=0; i<AVX_ARRAY_SIZE; i++)
         xPowY[i] = (double) x;

      y = xyPtr->powers[0].exponent;

      avx_powmod(xPowY, y, dps, reciprocals);
      
      avx_get_16a(xyPtr->powers[0].avxRemainders);

      prevY = y;
         
      for (yIndex=1; yIndex<xyPtr->powerCount; yIndex++)
      {
         y = xyPtr->powers[yIndex].exponent;
      
         // If x is even then y must be odd and if x is odd
         // then y must be even so y - prevY is always even
         powIndex = (y - prevY) >> 1;
         
         // We have x^prevY (mod p).
         // Now compute x^y (mod p) as (x^prevY * x^(y-prevY) (mod p)
         while (powIndex > maxPowers) 
         {
            avx_set_16b(ip_AvxPowers[maxPowers]);
            avx_mulmod(dps, reciprocals);
            powIndex -= maxPowers;
         };
         
         if (powIndex > 0)
         {
            avx_set_16b(ip_AvxPowers[powIndex]);
            avx_mulmod(dps, reciprocals);
         }
         
         avx_get_16a(xyPtr->powers[yIndex].avxRemainders);
         
         prevY = y;
      }
   }
}

void  XYYXWorker::CheckAvxXYRemainders(uint64_t *ps, double *dps, double *reciprocals) 
{
   double    __attribute__((aligned(32))) yPowX[AVX_ARRAY_SIZE];
   double    __attribute__((aligned(32))) t2[AVX_ARRAY_SIZE];
   uint32_t  x, y, prevX;
   uint32_t  xIndex, powIndex;
   uint32_t  maxPowers;
   base_t   *yxPtr;
   power_t  *pairedPower;

   // If the range of x is only 50, then we only want to generate
   // up to x^50 instead of x^100  (assuming MAX_POWERS = 50).   
   if (ii_XCount < MAX_POWERS * 2)
      maxPowers = (ii_XCount / 2);
   else
      maxPowers = MAX_POWERS;

   for (y=ii_MinY; y<=ii_MaxY; y++)
   {
      yxPtr = &ip_yxTerms[y-ii_MinY];
      
      if (yxPtr->powerCount == 0)
         continue;

      BuildAvxListOfPowers(y, dps, reciprocals, maxPowers);
      
      for (int i=0; i<AVX_ARRAY_SIZE; i++)
         yPowX[i] = (double) y;

      x = yxPtr->powers[0].exponent;

      avx_powmod(yPowX, x, dps, reciprocals);
      
      pairedPower = (power_t *) yxPtr->powers[0].pairedPower;
      
      CheckAvxResult(x, y, ps, dps, pairedPower);
         
      prevX = x;
         
      for (xIndex=1; xIndex<yxPtr->powerCount; xIndex++)
      {
         x = yxPtr->powers[xIndex].exponent;
         
         // If x is even then y must be odd and if x is odd
         // then y must be even so y - prevY is always even
         powIndex = (x - prevX) >> 1;
         
         // We have y^prevX (mod p).
         // Now compute y^x (mod p) as (y^prevX * y^(x-prevX) (mod p)
         while (powIndex > maxPowers) 
         {
            avx_set_16b(ip_AvxPowers[maxPowers]);
            avx_mulmod(dps, reciprocals);
            powIndex -= maxPowers;
         };
         
         if (powIndex > 0)
         {
            avx_set_16b(ip_AvxPowers[powIndex]);
            avx_mulmod(dps, reciprocals);
         }
         
         avx_get_16a(t2);

         pairedPower = (power_t *) yxPtr->powers[xIndex].pairedPower;
      
         CheckAvxResult(x, y, ps, dps, pairedPower);
         
         prevX = x;
      }
   }
}

void  XYYXWorker::BuildAvxListOfPowers(uint32_t base, double *dps, double *reciprocals, uint32_t count)
{
   uint32_t idx;
   double    __attribute__((aligned(32))) a[AVX_ARRAY_SIZE];

   for (int i=0; i<AVX_ARRAY_SIZE; i++)
   {
      a[i] = (double) base;
      ip_AvxPowers[0][i] = 1.0;
   }
      
   avx_set_16a(a);
   avx_set_16b(a);
   
   avx_mulmod(dps, reciprocals);
   
   // ip_AvxPowers[1] = a^2 mod p
   avx_get_16a(ip_AvxPowers[1]);
         
   // Multiply successive terms by a^2 (mod p)
   for (idx=2; idx<=count; idx++)
   {
      avx_set_16b(ip_AvxPowers[1]);
   
      avx_mulmod(dps, reciprocals);
      
      avx_get_16a(ip_AvxPowers[idx]);      
   }         
}

void  XYYXWorker::CheckAvxResult(uint32_t x, uint32_t y, uint64_t *ps, double *dps, power_t *pair)
{
   uint32_t idx;
   double __attribute__((aligned(32))) rems[AVX_ARRAY_SIZE];

   // Only go further if one or more of the 16 primes yielded a factor for this n
   if (ib_IsMinus && avx_pos_compare_16v(pair->avxRemainders) > 0)
   {
      avx_get_16a(rems);
            
      for (idx=0; idx<AVX_ARRAY_SIZE; idx++)
         if (rems[idx] == pair->avxRemainders[idx])
         {
            if (ip_XYYXApp->ReportFactor(ps[idx], x, y, -1))
               VerifyFactor(ps[idx], x, y, -1);
         }
   }
      
   // Only go further if one or more of the 16 primes yielded a factor for this n
   if (ib_IsPlus && avx_neg_compare_16v(pair->avxRemainders, dps))
   {
      avx_get_16a(rems);
      
      for (idx=0; idx<AVX_ARRAY_SIZE; idx++)
         if (rems[idx] == dps[idx] - pair->avxRemainders[idx])
         {
            if (ip_XYYXApp->ReportFactor(ps[idx], x, y, +1))
               VerifyFactor(ps[idx], x, y, +1);
         }
   }
}

void  XYYXWorker::VerifyFactor(uint64_t p, uint32_t x, uint32_t y, int32_t c)
{
   uint64_t xPowY, yPowX;
   
   fpu_push_1divp(p);
   
   xPowY = fpu_powmod(x, y, p);
   yPowX = fpu_powmod(y, x, p);
      
   if (c == -1 && xPowY != yPowX)
      FatalError("%" PRIu64" does not divide %u^%u-%u^%u", p, x, y, y, x);
   
   if (c == +1 && xPowY + yPowX != p)
      FatalError("%" PRIu64" does not divide %u^%u+%u^%u", p, x, y, y, x);
   
   fpu_pop();
}

