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

   ii_TermsElements = ip_XYYXApp->GetTermCount() + 2*ii_XCount + 2*ii_YCount;
   
   ip_xyTerms = (uint32_t *) xmalloc(ii_TermsElements * sizeof(uint32_t));
   ip_yxTerms = (uint32_t *) xmalloc(ii_TermsElements * sizeof(uint32_t));

   ip_AvxXYRemainders = 0;
   ip_XYRemainders = (uint64_t **) xmalloc(ii_XCount * sizeof(uint64_t *));

   for (uint32_t x=ii_MinX; x<=ii_MaxX; x++)
      ip_XYRemainders[X_INDEX(x)] = (uint64_t *) xmalloc(ii_YCount * sizeof(uint64_t));

   if (ip_XYYXApp->UseAvxIfAvailable() && CpuSupportsAvx())
   {
      ip_AvxXYRemainders = (double ***) xmalloc(ii_XCount * sizeof(double **));

      for (uint32_t x=ii_MinX; x<=ii_MaxX; x++)
      {
         ip_AvxXYRemainders[X_INDEX(x)] = (double **) xmalloc(ii_YCount * sizeof(double *));
         
         for (uint32_t y=ii_MinY; y<=ii_MaxY; y++)
            ip_AvxXYRemainders[X_INDEX(x)][Y_INDEX(y)] = (double *) xmalloc(AVX_ARRAY_SIZE * sizeof(double));
      }
            
      for (uint32_t i=0; i<=MAX_POWERS; i++)
         ip_Powers[i] = (double *) xmalloc(AVX_ARRAY_SIZE * sizeof(double));
   }
   
   il_NextTermsBuild = 0;
   ib_Initialized = true;
   
   if (ip_XYYXApp->UseAvxIfAvailable() && CpuSupportsAvx())
   {
      // The AVX logic does not handle when gcd(x, prime) > 1 or gcd(y, prime) > 1
      if (ii_MaxX < ii_MaxY)
         SetMiniChunkRange(ii_MaxY + 1, PMAX_MAX_52BIT, AVX_ARRAY_SIZE);
      else
         SetMiniChunkRange(ii_MaxX + 1, PMAX_MAX_52BIT, AVX_ARRAY_SIZE);
   }
}

void  XYYXWorker::CleanUp(void)
{
   xfree(ip_xyTerms);
   xfree(ip_yxTerms);

   if (CpuSupportsAvx())
   {
      for (uint32_t x=ii_MinX; x<=ii_MaxX; x++)
      {
         for (uint32_t y=ii_MinY; y<=ii_MaxY; y++)
            xfree(ip_AvxXYRemainders[X_INDEX(x)][Y_INDEX(y)]);
         
         xfree(ip_AvxXYRemainders[X_INDEX(x)]);
      }
      
      xfree(ip_AvxXYRemainders);
      
      for (uint32_t i=0; i<=MAX_POWERS; i++)
         xfree(ip_Powers[i]);
   }
   
   for (uint32_t x=ii_MinX; x<=ii_MaxX; x++)
      xfree(ip_XYRemainders[X_INDEX(x)]);
   
   xfree(ip_XYRemainders);
}

void  XYYXWorker::TestMegaPrimeChunk(void)
{
   uint64_t currentPrime = 0;
   uint64_t maxPrime = ip_App->GetMaxPrime();
   
   vector<uint64_t>::iterator it = iv_Primes.begin();

   while (it != iv_Primes.end())
   {
      currentPrime = *it;
      it++;
      
      // Every once in a while rebuild the term lists as it will have fewer entries
      // which will speed up testing for the next range of p.
      if (currentPrime > il_NextTermsBuild)
      {
         memset(ip_xyTerms, 0, ii_TermsElements * sizeof(uint32_t));
         memset(ip_yxTerms, 0, ii_TermsElements * sizeof(uint32_t));
         
         ip_XYYXApp->GetTerms(ip_xyTerms, ip_yxTerms);
         
         il_NextTermsBuild = (currentPrime << 1);
      }
      
      fpu_push_1divp(currentPrime);
      
      BuildXYRemainders(currentPrime);

      CheckXYRemainders(currentPrime);

      SetLargestPrimeTested(currentPrime, 1);
      
      fpu_pop();
      
      if (currentPrime >= maxPrime)
         break;
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
      memset(ip_xyTerms, 0, ii_TermsElements * sizeof(uint32_t));
      memset(ip_yxTerms, 0, ii_TermsElements * sizeof(uint32_t));
      
      ip_XYYXApp->GetTerms(ip_xyTerms, ip_yxTerms);
      
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
void  XYYXWorker::BuildXYRemainders(uint64_t p) 
{  
   uint32_t  x, y, prevY;
   uint32_t  pIndex;
   uint32_t *termList;
   uint64_t  xPowY, xPowN;
   uint64_t  powers[MAX_POWERS+1];
   uint32_t  maxPowers;
   
   // If the range of y is only 50, then we only want to generate
   // up to y^50 instead of y^100  (assuming MAX_POWERS = 50).
   if (ii_YCount < MAX_POWERS * 2)
      maxPowers = (ii_YCount / 2);
   else
      maxPowers = MAX_POWERS;

   termList = ip_xyTerms;

   while (*termList)
   {
      x = *termList;
      termList++;

      if (*termList)
      {
         y = *termList;
         termList++;

         BuildListOfPowers(x, p, maxPowers, powers);

         xPowY = fpu_powmod(x, y, p);

         ip_XYRemainders[X_INDEX(x)][Y_INDEX(y)] = xPowY; 
         
         prevY = y;
         
         while (*termList)
         {
            y = *termList;
            termList++;
      
            // If x is even then y must be odd and if x is odd
            // then y must be even so y - prevY is always even
            pIndex = (y - prevY) >> 1;
            
            if (pIndex >= maxPowers)
               xPowN = fpu_powmod(x, y - prevY, p);
            else
               xPowN = powers[pIndex];

            // We have x^prevY (mod p).
            // Now compute x^y (mod p) as (x^prevY * x^(y-prevY) (mod p)
            xPowY = fpu_mulmod(xPowY, xPowN, p);

            ip_XYRemainders[X_INDEX(x)][Y_INDEX(y)] = xPowY;
            
            prevY = y;
         }
      }

      // Skip the 0 term
      termList++;
   }
}

void  XYYXWorker::CheckXYRemainders(uint64_t p) 
{
   uint32_t  x, y, prevX;
   uint32_t  pIndex;
   uint64_t  yPowX, xPowY, yPowN;
   uint32_t *termList;
   uint64_t  powers[MAX_POWERS+1];
   uint32_t  maxPowers;

   // If the range of x is only 50, then we only want to generate
   // up to x^50 instead of x^100  (assuming MAX_POWERS = 50).   
   if (ii_XCount < MAX_POWERS * 2)
      maxPowers = (ii_XCount / 2);
   else
      maxPowers = MAX_POWERS;

   termList = ip_yxTerms;

   while (*termList)
   {
      y = *termList;
      termList++;
      
      if (*termList)
      {
         x = *termList;
         termList++;
         
         BuildListOfPowers(y, p, maxPowers, powers);

         yPowX = fpu_powmod(y, x, p);

         xPowY = ip_XYRemainders[X_INDEX(x)][Y_INDEX(y)];
                     
         // Note the terms where (xPowY + yPowX) (mod p) = 0 have already been
         // removed as gcd(x, y) > 1 for those x and y
         if (ib_IsPlus && xPowY + yPowX == p)
            if (ip_XYYXApp->ReportFactor(p, x, y, +1))
               VerifyFactorFPU(p, x, y, +1);
            
         if (ib_IsMinus && xPowY == yPowX)
            if (ip_XYYXApp->ReportFactor(p, x, y, -1))
               VerifyFactorFPU(p, x, y, -1);
         
         prevX = x;
         
         while (*termList)
         {
            x = *termList;
            termList++;
            
            // If x is even then y must be odd and if x is odd
            // then y must be even so y - prevY is always even
            pIndex = (x - prevX) >> 1;
            
            if (pIndex >= maxPowers)
               yPowN = fpu_powmod(y, x - prevX, p);
            else
               yPowN = powers[pIndex];

            // We have y^prevX (mod p).
            // Now compute y^x (mod p) as (y^prevX * y^(x-prevX) (mod p)
            yPowX = fpu_mulmod(yPowX, yPowN, p);

            xPowY = ip_XYRemainders[X_INDEX(x)][Y_INDEX(y)];
            
            // Note the terms where (xPowY + yPowX) (mod p) = 0 have already been
            // removed as gcd(x, y) > 1 for those x and y
            if (ib_IsPlus && xPowY + yPowX == p)
               if (ip_XYYXApp->ReportFactor(p, x, y, +1))
                  VerifyFactorFPU(p, x, y, +1);
               
            if (ib_IsMinus && xPowY == yPowX)
               if (ip_XYYXApp->ReportFactor(p, x, y, -1))
                  VerifyFactorFPU(p, x, y, -1);
            
            prevX = x;
         }
      }
      
      // Skip the 0 term
      termList++;
   }
}

// Build a list of powers for a from a^0 thru a^n for all even n up to (count * 2).
void  XYYXWorker::BuildListOfPowers(uint64_t a, uint64_t p, uint32_t count, uint64_t *powers)
{
   uint32_t index;
   uint64_t squaredMod;
   
   memset(powers, 0, MAX_POWERS * sizeof(uint64_t));

   squaredMod = fpu_mulmod(a, a, p);
   
   fpu_push_adivb(squaredMod, p);
   
   powers[0] = 1;
   
   // powers[1] = a^2 mod p
   powers[1] = squaredMod;
   
   // Multiply successive terms by a^2 (mod p)
   for (index=2; index<=count; index++)
       powers[index] = fpu_mulmod_iter(powers[index-1], squaredMod, p);

   fpu_pop();
}

// Build a table of x^y mod p for all remaining terms
void  XYYXWorker::BuildAvxXYRemainders(uint64_t *ps, double *dps, double *reciprocals) 
{  
   double    __attribute__((aligned(32))) xPowY[AVX_ARRAY_SIZE];
   uint32_t  x, y, prevY;
   uint32_t  pIndex;
   uint32_t *termList;
   uint32_t  maxPowers;
   
   // If the range of y is only 50, then we only want to generate
   // up to y^50 instead of y^100  (assuming MAX_POWERS = 50).
   if (ii_YCount < MAX_POWERS * 2)
      maxPowers = (ii_YCount / 2);
   else
      maxPowers = MAX_POWERS;

   termList = ip_xyTerms;

   while (*termList)
   {
      x = *termList;
      termList++;

      if (*termList)
      {
         y = *termList;
         termList++;

         BuildAvxListOfPowers(x, dps, reciprocals, maxPowers);

         for (int i=0; i<AVX_ARRAY_SIZE; i++)
            xPowY[i] = (double) x;

         avx_powmod(xPowY, y, dps, reciprocals);
         
         avx_get_16a(ip_AvxXYRemainders[X_INDEX(x)][Y_INDEX(y)]);

         prevY = y;
         
         while (*termList)
         {
            y = *termList;
            termList++;
      
            // If x is even then y must be odd and if x is odd
            // then y must be even so y - prevY is always even
            pIndex = (y - prevY) >> 1;
            
            // We have x^prevY (mod p).
            // Now compute x^y (mod p) as (x^prevY * x^(y-prevY) (mod p)
            while (pIndex > maxPowers) 
            {
               avx_set_16b(ip_Powers[maxPowers]);
               avx_mulmod(dps, reciprocals);
               pIndex -= maxPowers;
            };
            
            if (pIndex > 0)
            {
               avx_set_16b(ip_Powers[pIndex]);
               avx_mulmod(dps, reciprocals);
            }
            
            avx_get_16a(ip_AvxXYRemainders[X_INDEX(x)][Y_INDEX(y)]);

            prevY = y;
         }
      }

      // Skip the 0 term
      termList++;
   }
}

void  XYYXWorker::CheckAvxXYRemainders(uint64_t *ps, double *dps, double *reciprocals) 
{
   double    __attribute__((aligned(32))) yPowX[AVX_ARRAY_SIZE];
   double    __attribute__((aligned(32))) t2[AVX_ARRAY_SIZE];
   uint32_t  x, y, prevX;
   uint32_t  pIndex;
   uint32_t *termList;
   uint32_t  maxPowers;

   // If the range of x is only 50, then we only want to generate
   // up to x^50 instead of x^100  (assuming MAX_POWERS = 50).   
   if (ii_XCount < MAX_POWERS * 2)
      maxPowers = (ii_XCount / 2);
   else
      maxPowers = MAX_POWERS;

   termList = ip_yxTerms;

   while (*termList)
   {
      y = *termList;
      termList++;
      
      if (*termList)
      {
         x = *termList;
         termList++;
    
         BuildAvxListOfPowers(y, dps, reciprocals, maxPowers);
         
         //if (y == 2067)
         //   for (int i=0; i<AVX_ARRAY_SIZE; i++)
         //      for (int j=0; j<maxPowers; j++)
         //         printf("%d^%d mod %llu %.0lf \n", y, j*2, ps[i], ip_Powers[i][j]);
        
         for (int i=0; i<AVX_ARRAY_SIZE; i++)
            yPowX[i] = (double) y;
         
         avx_powmod(yPowX, x, dps, reciprocals);
         
         CheckAvxResult(x, y, ps, dps);
         
         prevX = x;
         
         while (*termList)
         {
            x = *termList;
            termList++;
            
            // If x is even then y must be odd and if x is odd
            // then y must be even so y - prevY is always even
            pIndex = (x - prevX) >> 1;
            
            // We have y^prevX (mod p).
            // Now compute y^x (mod p) as (y^prevX * y^(x-prevX) (mod p)
            while (pIndex > maxPowers) 
            {
               avx_set_16b(ip_Powers[maxPowers]);
               avx_mulmod(dps, reciprocals);
               pIndex -= maxPowers;
            };
            
            if (pIndex > 0)
            {
               avx_set_16b(ip_Powers[pIndex]);
               avx_mulmod(dps, reciprocals);
            }
            
            avx_get_16a(t2);

            CheckAvxResult(x, y, ps, dps);
            
            prevX = x;
         }
      }
      
      // Skip the 0 term
      termList++;
   }
}

void  XYYXWorker::BuildAvxListOfPowers(uint32_t base, double *dps, double *reciprocals, uint32_t count)
{
   uint32_t idx;
   double    __attribute__((aligned(32))) a[AVX_ARRAY_SIZE];

   for (int i=0; i<AVX_ARRAY_SIZE; i++)
   {
      a[i] = (double) base;
      ip_Powers[0][i] = 1.0;
   }
      
   avx_set_16a(a);
   avx_set_16b(a);
   
   avx_mulmod(dps, reciprocals);
   
   // ip_Powers[1] = a^2 mod p
   avx_get_16a(ip_Powers[1]);
         
   // Multiply successive terms by a^2 (mod p)
   for (idx=2; idx<=count; idx++)
   {
      avx_set_16b(ip_Powers[1]);
   
      avx_mulmod(dps, reciprocals);
      
      avx_get_16a(ip_Powers[idx]);      
   }         
}

void  XYYXWorker::CheckAvxResult(uint32_t x, uint32_t y, uint64_t *ps, double *dps)
{
   uint32_t idx;
   double __attribute__((aligned(32))) rems[AVX_ARRAY_SIZE];

   // Only go further if one or more of the 16 primes yielded a factor for this n
   if (ib_IsMinus && avx_pos_compare_16v(ip_AvxXYRemainders[X_INDEX(x)][Y_INDEX(y)]) > 0)
   {
      avx_get_16a(rems);
            
      for (idx=0; idx<AVX_ARRAY_SIZE; idx++)
         if (rems[idx] == ip_AvxXYRemainders[X_INDEX(x)][Y_INDEX(y)][idx])
         {
            if (ip_XYYXApp->ReportFactor(ps[idx], x, y, -1))
               VerifyFactorAVX(ps[idx], x, y, -1);
         }
   }
      
   // Only go further if one or more of the 16 primes yielded a factor for this n
   if (ib_IsPlus && avx_neg_compare_16v(ip_AvxXYRemainders[X_INDEX(x)][Y_INDEX(y)], dps))
   {
      avx_get_16a(rems);
      
      for (idx=0; idx<AVX_ARRAY_SIZE; idx++)
         if (rems[idx] == dps[idx] - ip_AvxXYRemainders[X_INDEX(x)][Y_INDEX(y)][idx])
         {
            if (ip_XYYXApp->ReportFactor(ps[idx], x, y, +1))
               VerifyFactorAVX(ps[idx], x, y, +1);
         }
   }
}            
            
void  XYYXWorker::VerifyFactorFPU(uint64_t p, uint32_t x, uint32_t y, int32_t c)
{
   uint64_t xPowY, yPowX;
   
   xPowY = fpu_powmod(x, y, p);
   yPowX = fpu_powmod(y, x, p);
      
   if (c == -1 && xPowY != yPowX)
      FatalError("%" PRIu64" does not divide %u^%u-%u^%u", p, x, y, y, x);
   
   if (c == +1 && xPowY + yPowX != p)
      FatalError("%" PRIu64" does not divide %u^%u+%u^%u", p, x, y, y, x);
}
            
void  XYYXWorker::VerifyFactorAVX(uint64_t p, uint32_t x, uint32_t y, int32_t c)
{
   fpu_push_1divp(p);
   
   VerifyFactorFPU(p, x, y, c);
   
   fpu_pop();
}

