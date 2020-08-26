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
   il_NextTermsBuild = 0;
   
   ib_Initialized = true;

   for (uint32_t i=0; i<=MAX_POWERS; i++)
      ip_FpuPowers[i] = (uint64_t *) malloc(4 * sizeof(uint64_t));
            
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
   
   for (uint32_t i=0; i<=MAX_POWERS; i++)
      xfree(ip_FpuPowers[i]);
}

void  XYYXWorker::FreeTerms(void)
{
   uint32_t x, y;
   
   if (ip_xyTerms == 0)
      return;
   
   for (x=ii_MinX; x<=ii_MaxX; x++)
   {
      if (ip_xyTerms[x-ii_MinX].powerCount == 0)
         continue;
      
      if (ib_HaveFpuRemainders)
         xfree(ip_xyTerms[x-ii_MinX].fpuRemaindersPtr);

      if (ib_HaveAvxRemainders)
         xfree(ip_xyTerms[x-ii_MinX].avxRemaindersPtr);
      
      xfree(ip_xyTerms[x-ii_MinX].powersOfX);
   }
   
   for (y=ii_MinY; y<=ii_MaxY; y++)
   {
      if (ip_yxTerms[y-ii_MinY].powerCount == 0)
         continue;

      xfree(ip_yxTerms[y-ii_MinY].powersOfY);
   }
      
   xfree(ip_xyTerms);
   xfree(ip_yxTerms);
   
   ip_xyTerms = 0;
   ip_yxTerms = 0;
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
         
         ib_HaveFpuRemainders = true;
         ib_HaveAvxRemainders = false;
         ip_XYYXApp->GetTerms(4, 0, &bases);
               
         ip_xyTerms = bases.xPowY;
         ip_yxTerms = bases.yPowX;
      
         il_NextTermsBuild = (ps[3] << 2);
      }
      
      // Compute x^y for y
      BuildFpuXYRemainders(ps);

      CheckFpuXYRemainders(ps);

      SetLargestPrimeTested(ps[3], 4);
      
      if (ps[3] > maxPrime)
         break;
   }
}

// Build a table of x^y mod p for all remaining terms
void  XYYXWorker::BuildFpuXYRemainders(uint64_t *ps) 
{  
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

      BuildFpuListOfPowers(x, ps, maxPowers);
      
      xyPtr->powersOfX[0].fpuRemainders[0] = x;
      xyPtr->powersOfX[0].fpuRemainders[1] = x;
      xyPtr->powersOfX[0].fpuRemainders[2] = x;
      xyPtr->powersOfX[0].fpuRemainders[3] = x;

      y = xyPtr->powersOfX[0].y;
      
      fpu_powmod_4b_1n_4p(xyPtr->powersOfX[0].fpuRemainders, y, ps);

      fpu_push_1divp(ps[3]);
      fpu_push_1divp(ps[2]);
      fpu_push_1divp(ps[1]);
      fpu_push_1divp(ps[0]);
      
      prevY = y;
         
      for (yIndex=1; yIndex<xyPtr->powerCount; yIndex++)
      {
         y = xyPtr->powersOfX[yIndex].y;
      
         // If x is even then y must be odd and if x is odd
         // then y must be even so y - prevY is always even
         powIndex = (y - prevY) >> 1;
         
         xyPtr->powersOfX[yIndex].fpuRemainders[0] = xyPtr->powersOfX[yIndex-1].fpuRemainders[0];
         xyPtr->powersOfX[yIndex].fpuRemainders[1] = xyPtr->powersOfX[yIndex-1].fpuRemainders[1];
         xyPtr->powersOfX[yIndex].fpuRemainders[2] = xyPtr->powersOfX[yIndex-1].fpuRemainders[2];
         xyPtr->powersOfX[yIndex].fpuRemainders[3] = xyPtr->powersOfX[yIndex-1].fpuRemainders[3];
         
         // We have x^prevY (mod p).
         // Now compute x^y (mod p) as (x^prevY * x^(y-prevY) (mod p)
         while (powIndex > maxPowers) 
         {
            fpu_mulmod_4a_4b_4p(xyPtr->powersOfX[yIndex].fpuRemainders, ip_FpuPowers[maxPowers], ps);
            powIndex -= maxPowers;
         };
         
         if (powIndex > 0)
            fpu_mulmod_4a_4b_4p(xyPtr->powersOfX[yIndex].fpuRemainders, ip_FpuPowers[powIndex], ps);
         
         prevY = y;
      }
      
      fpu_pop();
      fpu_pop();
      fpu_pop();
      fpu_pop();
   }
}

void  XYYXWorker::CheckFpuXYRemainders(uint64_t *ps) 
{
   uint32_t    x, y, prevX;
   uint32_t    xIndex, powIndex;
   uint32_t    maxPowers;
   uint64_t    yPowXRemainders[4];
   base_t     *yxPtr;
   powerofx_t *powerOfX;

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

      BuildFpuListOfPowers(y, ps, maxPowers);
      
      yPowXRemainders[0] = y;
      yPowXRemainders[1] = y;
      yPowXRemainders[2] = y;
      yPowXRemainders[3] = y;

      x = yxPtr->powersOfY[0].x;
      
      fpu_powmod_4b_1n_4p(yPowXRemainders, x, ps);
      
      powerOfX = yxPtr->powersOfY[0].powerOfX;
      
      CheckFpuResult(x, y, ps, powerOfX->fpuRemainders, yPowXRemainders);
         
      fpu_push_1divp(ps[3]);
      fpu_push_1divp(ps[2]);
      fpu_push_1divp(ps[1]);
      fpu_push_1divp(ps[0]);

      prevX = x;
      
      for (xIndex=1; xIndex<yxPtr->powerCount; xIndex++)
      {
         x = yxPtr->powersOfY[xIndex].x;
                  
         // If x is even then y must be odd and if x is odd
         // then y must be even so y - prevY is always even
         powIndex = (x - prevX) >> 1;
         
         // We have x^prevY (mod p).
         // Now compute x^y (mod p) as (x^prevY * x^(y-prevY) (mod p)
         while (powIndex > maxPowers) 
         {
            fpu_mulmod_4a_4b_4p(yPowXRemainders, ip_FpuPowers[maxPowers], ps);
            powIndex -= maxPowers;
         };
         
         if (powIndex > 0)
            fpu_mulmod_4a_4b_4p(yPowXRemainders, ip_FpuPowers[powIndex], ps);
                  
         powerOfX = yxPtr->powersOfY[xIndex].powerOfX;
      
         CheckFpuResult(x, y, ps, powerOfX->fpuRemainders, yPowXRemainders);
         
         prevX = x;
      }
      
      fpu_pop();
      fpu_pop();
      fpu_pop();
      fpu_pop();
   }
}

void  XYYXWorker::BuildFpuListOfPowers(uint32_t base, uint64_t *ps, uint32_t count)
{
   uint32_t idx;

   fpu_push_1divp(ps[3]);
   fpu_push_1divp(ps[2]);
   fpu_push_1divp(ps[1]);
   fpu_push_1divp(ps[0]);
   
   ip_FpuPowers[0][0] = base;
   ip_FpuPowers[0][1] = base;
   ip_FpuPowers[0][2] = base;
   ip_FpuPowers[0][3] = base;

   ip_FpuPowers[1][0] = base;
   ip_FpuPowers[1][1] = base;
   ip_FpuPowers[1][2] = base;
   ip_FpuPowers[1][3] = base;
   
   fpu_mulmod_4a_4b_4p(ip_FpuPowers[1], ip_FpuPowers[0], ps);

   // Multiply successive terms by a^2 (mod p)
   for (idx=2; idx<=count; idx++)
   {
      ip_FpuPowers[idx][0] = ip_FpuPowers[idx-1][0];
      ip_FpuPowers[idx][1] = ip_FpuPowers[idx-1][1];
      ip_FpuPowers[idx][2] = ip_FpuPowers[idx-1][2];
      ip_FpuPowers[idx][3] = ip_FpuPowers[idx-1][3];
      
      fpu_mulmod_4a_4b_4p(ip_FpuPowers[idx], ip_FpuPowers[1], ps);
   }
   
   fpu_pop();
   fpu_pop();
   fpu_pop();
   fpu_pop();
}

void  XYYXWorker::CheckFpuResult(uint32_t x, uint32_t y, uint64_t *ps, uint64_t *powerOfX, uint64_t *powerOfY)
{
   uint32_t idx;

   if (ib_IsMinus)
   {            
      for (idx=0; idx<4; idx++)
         if (powerOfX[idx] == powerOfY[idx])
         {
            if (ip_XYYXApp->ReportFactor(ps[idx], x, y, -1))
               VerifyFactor(ps[idx], x, y, -1);
         }
   }

   if (ib_IsPlus)
   {
      for (idx=0; idx<4; idx++)
         if (powerOfX[idx] == ps[idx] - powerOfY[idx])
         {
            if (ip_XYYXApp->ReportFactor(ps[idx], x, y, +1))
               VerifyFactor(ps[idx], x, y, +1);
         }
   }
}

void  XYYXWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   double __attribute__((aligned(32))) dps[AVX_ARRAY_SIZE];
   double __attribute__((aligned(32))) reciprocals[AVX_ARRAY_SIZE];
   
   // Every once in a while rebuild the term lists as it will have fewer entries
   // which will speed up testing for the next range of p.
   if (!ib_HaveAvxRemainders || miniPrimeChunk[0] > il_NextTermsBuild)
   {
      bases_t bases;
      
      FreeTerms();
      
      ib_HaveFpuRemainders = false;
      ib_HaveAvxRemainders = true;
      
      ip_XYYXApp->GetTerms(0, AVX_ARRAY_SIZE, &bases);
            
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

      y = xyPtr->powersOfX[0].y;

      avx_powmod(xPowY, y, dps, reciprocals);
      
      avx_get_16a(xyPtr->powersOfX[0].avxRemainders);

      prevY = y;
         
      for (yIndex=1; yIndex<xyPtr->powerCount; yIndex++)
      {
         y = xyPtr->powersOfX[yIndex].y;
      
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
         
         avx_get_16a(xyPtr->powersOfX[yIndex].avxRemainders);
         
         prevY = y;
      }
   }
}

void  XYYXWorker::CheckAvxXYRemainders(uint64_t *ps, double *dps, double *reciprocals) 
{
   double      __attribute__((aligned(32))) yPowX[AVX_ARRAY_SIZE];
   uint32_t    x, y, prevX;
   uint32_t    xIndex, powIndex;
   uint32_t    maxPowers;
   base_t     *yxPtr;
   powerofx_t *powerOfX;

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

      x = yxPtr->powersOfY[0].x;

      avx_powmod(yPowX, x, dps, reciprocals);
      
      powerOfX = yxPtr->powersOfY[0].powerOfX;
      
      CheckAvxResult(x, y, ps, dps, powerOfX->avxRemainders);
         
      prevX = x;
         
      for (xIndex=1; xIndex<yxPtr->powerCount; xIndex++)
      {
         x = yxPtr->powersOfY[xIndex].x;
         
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
         
         powerOfX = yxPtr->powersOfY[xIndex].powerOfX;
      
         CheckAvxResult(x, y, ps, dps, powerOfX->avxRemainders);
         
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

void  XYYXWorker::CheckAvxResult(uint32_t x, uint32_t y, uint64_t *ps, double *dps, double *powerOfX)
{
   uint32_t idx;
   double __attribute__((aligned(32))) rems[AVX_ARRAY_SIZE];

   // Only go further if one or more of the 16 primes yielded a factor for this n
   if (ib_IsMinus && avx_pos_compare_16v(powerOfX) > 0)
   {
      avx_get_16a(rems);
            
      for (idx=0; idx<AVX_ARRAY_SIZE; idx++)
         if (rems[idx] == powerOfX[idx])
         {
            if (ip_XYYXApp->ReportFactor(ps[idx], x, y, -1))
               VerifyFactor(ps[idx], x, y, -1);
         }
   }
      
   // Only go further if one or more of the 16 primes yielded a factor for this n
   if (ib_IsPlus && avx_neg_compare_16v(powerOfX, dps))
   {
      avx_get_16a(rems);
      
      for (idx=0; idx<AVX_ARRAY_SIZE; idx++)
         if (rems[idx] == dps[idx] - powerOfX[idx])
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

