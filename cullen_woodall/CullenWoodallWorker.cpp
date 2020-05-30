/* CullenWoodallWorker.cpp -- (C) Mark Rodenkirch, May 2018

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <time.h>

#include "CullenWoodallWorker.h"
#include "../x86_asm/fpu-asm-x86.h"
#include "../x86_asm/avx-asm-x86.h"

// This is for building a list of even powers for b
#define MAX_POWERS   50

#define X_INDEX(x)  ((x) - ii_MinX)
#define Y_INDEX(y)  ((y) - ii_MinY)

CullenWoodallWorker::CullenWoodallWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{ 
   ip_CullenWoodallApp = (CullenWoodallApp *) theApp;
   
   ii_Base = ip_CullenWoodallApp->GetBase();
   ii_MinN = ip_CullenWoodallApp->GetMinN();
   ii_MaxN = ip_CullenWoodallApp->GetMaxN();
   
   ii_MaxTermCount = ip_CullenWoodallApp->GetMaxN() - ip_CullenWoodallApp->GetMinN() + 10;

   // Allocate enough memory to hold all of the terms.
   ii_Terms = (uint32_t *) xmalloc(ii_MaxTermCount*sizeof(int32_t));

   il_NextTermsBuild = 0;
   ib_Initialized = true;
   
   if (CpuSupportsAvx())
   {
      if (ii_Base < ii_MaxN)
         SetMiniChunkRange(ii_MaxN + 1, PMAX_MAX_52BIT, AVX_ARRAY_SIZE);
      else
         SetMiniChunkRange(ii_Base + 1, PMAX_MAX_52BIT, AVX_ARRAY_SIZE);
   }
}

void  CullenWoodallWorker::CleanUp(void)
{   
   xfree(ii_Terms);
}

void  CullenWoodallWorker::TestMegaPrimeChunk(void)
{
   uint64_t maxPrime = ip_CullenWoodallApp->GetMaxPrime();
   uint64_t ps[4];
   uint32_t maxPForSmallPrimeLogic = ((ii_Base < ii_MaxN) ? (ii_MaxN + 1) : (ii_Base + 1));
   
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
      // Unlike the GPU, we will put all terms into one group.
      if (ps[0] > il_NextTermsBuild)
      {
         ip_CullenWoodallApp->GetTerms(ii_Terms, ii_MaxTermCount, ii_MaxTermCount);
         
         il_NextTermsBuild = (ps[3] << 1);
      }
      
      if (ps[0] < maxPForSmallPrimeLogic)
         TestSmallPrimesFPU(ps);
      else
         TestLargePrimesFPU(ps);
      
      SetLargestPrimeTested(ps[3], 4);
      
      if (ps[3] > maxPrime)
         break;
   }
}

void  CullenWoodallWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   // Every once in a while rebuild the term lists as it will have fewer entries
   // which will speed up testing for the next range of p.
   // Unlike the GPU, we will put all terms into one group.
   if (miniPrimeChunk[0] > il_NextTermsBuild)
   {
      ip_CullenWoodallApp->GetTerms(ii_Terms, ii_MaxTermCount, ii_MaxTermCount);
      
      il_NextTermsBuild = (miniPrimeChunk[AVX_ARRAY_SIZE-1] << 1);
   }

   TestPrimesAVX(miniPrimeChunk);
}

// For small p, we need to iterate from minn to maxn as the large prime algorithm
// does not work correctly when p < maxN.
void  CullenWoodallWorker::TestSmallPrimesFPU(uint64_t *ps)
{
   uint32_t theN, prevN;
   int32_t  termIndex, maxIndex;
   uint32_t power;
   uint64_t thePrime;
   uint64_t powrem, cwrem, powers[MAX_POWERS+1];

   maxIndex = 0;
   while (ii_Terms[maxIndex] > 0)
      maxIndex++;
      
   for (int i=0; i<4; i++)
   {
      thePrime = ps[i];
      
      if (ii_Base % thePrime == 0)
         continue;
      
      fpu_push_1divp(thePrime);
      
      termIndex = maxIndex - 1;

      prevN = theN = ii_Terms[termIndex];
         
      BuildListOfPowers(ii_Base, thePrime, MAX_POWERS+1, powers);
   
      // Compute our starting term
      powrem = fpu_powmod(ii_Base, theN, thePrime);

      termIndex--;
      
      // Note that the terms start at max N and decrease
      termIndex--;
      while (termIndex >= 0)
      {
         theN = ii_Terms[termIndex];
         
         power = theN - prevN;

         // Our table isn't infinite in size, so we'll break down the power into smaller
         // chunks.  At worst this might cause an extra mulmod or two every once in a while.
         while (power > MAX_POWERS)
         {
            powrem = fpu_mulmod(powrem, powers[MAX_POWERS], thePrime);
            power -= MAX_POWERS;
         }

         powrem = fpu_mulmod(powrem, powers[power], thePrime);

         cwrem = fpu_mulmod(powrem, theN, thePrime);
         
         if (cwrem == +1)
            if (ip_CullenWoodallApp->ReportFactor(thePrime, theN, -1))
               VerifyFactor(thePrime, theN, -1);
            
         if (cwrem == thePrime - 1)
            if (ip_CullenWoodallApp->ReportFactor(thePrime, theN, +1))
               VerifyFactor(thePrime, theN, +1);

         prevN = theN;
         termIndex--;
      };
      
      fpu_pop();
   }
}

// We're trying to determine if n*b^n (mod p) = +1/-1.  This requires an two
// operations, an exponentiation (b^n) and a multiplcation (n).  We can make a
// change to eliminate one of those operations.  Doing this requires us to
// compute the multiplicative inverse of b.  The multiplicated inverse of
// b (which I will call B) is the number such at b*B = 1 mod p.  Here is how
// the algorithm works.
//
// We start with our number:
//       n*b^n (mod p) = +1/-1
//    -> n*b^n*B^n (mod p) = +1/-1 * B^n (mod p)
//    -> n (mod p) = +1/-1 B^n (mod p)
//    -> +n/-n (mod p) = B^n (mod p)
// And voila, if B^n (mod p) = -n, then p is a factor of n*b^n+1.
// and if B^n (mod p) = +n, then p is a factor of n*b^n-1.
//
// To get then next term (which I will call N) where n > N, we need to
// multiply B^n by B^(N-n):
//    -> +n/-n (mod p) = B^n * B(N-n) (mod p)
//
// Since N-n is negative we could compute (p-1)+N-n, but that could be
// a very large exponent.  Instead since B is the multiplicative inverse
// of b, we can use that, i.e. b^(n-N), which will have a much smaller
// exponent and that allows us to put b^(n-N) into a table so that we
// can use mulmods instead of expmods as we iterate through n.
//    -> +n/-n (mod p) = B^n * b^n (mod p)
void  CullenWoodallWorker::TestLargePrimesFPU(uint64_t *ps)
{
   uint32_t theN, prevN;
   uint32_t termIndex;
   uint64_t powinvs[4];
   uint64_t thePrime;
   uint64_t powers[MAX_POWERS+1][4];
   uint64_t rems[4];
   uint32_t power;
   
   // compute the inverse of b (mod p)
   powinvs[0] = ComputeMultiplicativeInverse(ii_Base, ps[0]);
   powinvs[1] = ComputeMultiplicativeInverse(ii_Base, ps[1]);
   powinvs[2] = ComputeMultiplicativeInverse(ii_Base, ps[2]);
   powinvs[3] = ComputeMultiplicativeInverse(ii_Base, ps[3]);

   fpu_powmod_4b_1n_4p(powinvs, ii_Terms[0], ps);
   
   for (int i=0; i<4; i++)
   {
      thePrime = ps[i];
      
      theN = ii_Terms[0];
      
      if (powinvs[i] == theN)
         if (ip_CullenWoodallApp->ReportFactor(thePrime, theN, -1))
            VerifyFactor(thePrime, theN, -1);
         
      if (powinvs[i] == thePrime - theN)
         if (ip_CullenWoodallApp->ReportFactor(thePrime, theN, +1))
            VerifyFactor(thePrime, theN, +1);
   }
   
   fpu_push_1divp(ps[3]);
   fpu_push_1divp(ps[2]);
   fpu_push_1divp(ps[1]);
   fpu_push_1divp(ps[0]);
   
   powers[1][0] = ii_Base;
   powers[1][1] = ii_Base;
   powers[1][2] = ii_Base;
   powers[1][3] = ii_Base;
   
   // Multiply successive terms by a (mod p)
   for (uint32_t idx=2; idx<MAX_POWERS+1; idx++)
   {
      // If the bsae is odd, then all n must be even and thus the difference
      // between any two remaining n for the base must also be even.
      if (idx > 2 && ii_Base & 1)
      {
         if (idx & 1)
            continue;
         
         rems[0] = powers[idx][0] = powers[idx-2][0];
         rems[1] = powers[idx][1] = powers[idx-2][1];
         rems[2] = powers[idx][2] = powers[idx-2][2];
         rems[3] = powers[idx][3] = powers[idx-2][3];
         
         fpu_mulmod_4a_4b_4p(powers[idx], powers[2], ps);
      }
      else 
      {
         rems[0] = powers[idx][0] = powers[idx-1][0];
         rems[1] = powers[idx][1] = powers[idx-1][1];
         rems[2] = powers[idx][2] = powers[idx-1][2];
         rems[3] = powers[idx][3] = powers[idx-1][3];
      
         fpu_mulmod_4a_4b_4p(powers[idx], powers[1], ps);
      }
   }
   
   rems[0] = powinvs[0];
   rems[1] = powinvs[1];
   rems[2] = powinvs[2];
   rems[3] = powinvs[3];
   
   prevN = ii_Terms[0];
   
   // Note that the terms start at max N and decrease
   termIndex = 1;
   while (ii_Terms[termIndex] > 0)
   {
      theN = ii_Terms[termIndex];

      power = prevN - theN;

      // Our table isn't infinite in size, so we'll break down the power into smaller
      // chunks.  At worst this might cause an extra mulmod or two every once in a while.
      while (power > MAX_POWERS)
      {
         fpu_mulmod_4a_4b_4p(rems, powers[MAX_POWERS], ps);
         power -= MAX_POWERS;
      }

      fpu_mulmod_4a_4b_4p(rems, powers[power], ps);

      // At this point we have computed (1/b)^n (mod p).
      // If (1/b)^n (mod p) == n then we have a Woodall factor.
      // If (1/b)^n (mod p) == thePrime - n then we have a Cullen factor.

      if (rems[0] == theN)
         if (ip_CullenWoodallApp->ReportFactor(ps[0], theN, -1))
            VerifyFactor(ps[0], theN, -1);
         
      if (rems[0] == ps[0] - theN)
         if (ip_CullenWoodallApp->ReportFactor(ps[0], theN, +1))
            VerifyFactor(ps[0], theN, +1);
         
      if (rems[1] == theN)
         if (ip_CullenWoodallApp->ReportFactor(ps[1], theN, -1))
            VerifyFactor(ps[1], theN, -1);
         
      if (rems[1] == ps[1] - theN)
         if (ip_CullenWoodallApp->ReportFactor(ps[1], theN, +1))
            VerifyFactor(ps[1], theN, +1);
         
      if (rems[2] == theN)
         if (ip_CullenWoodallApp->ReportFactor(ps[2], theN, -1))
            VerifyFactor(ps[2], theN, -1);
         
      if (rems[2] == ps[2] - theN)
         if (ip_CullenWoodallApp->ReportFactor(ps[2], theN, +1))
            VerifyFactor(ps[2], theN, +1);
         
      if (rems[3] == theN)
         if (ip_CullenWoodallApp->ReportFactor(ps[3], theN, -1))
            VerifyFactor(ps[3], theN, -1);
         
      if (rems[3] == ps[3] - theN)
         if (ip_CullenWoodallApp->ReportFactor(ps[3], theN, +1))
            VerifyFactor(ps[3], theN, +1);
         
      prevN = theN;
      termIndex++;
   };
   
   fpu_pop();
   fpu_pop();
   fpu_pop();
   fpu_pop();
}

// Same as TestLargePrimesFPU, but using AVX
void  CullenWoodallWorker::TestPrimesAVX(uint64_t *ps)
{
   uint32_t theN, prevN;
   uint32_t termIndex;
   uint32_t power;
   double __attribute__((aligned(32))) powers[MAX_POWERS+1][AVX_ARRAY_SIZE];
   double __attribute__((aligned(32))) dps[AVX_ARRAY_SIZE];
   double __attribute__((aligned(32))) reciprocals[AVX_ARRAY_SIZE];
   double __attribute__((aligned(32))) multinvs[AVX_ARRAY_SIZE];
   
   // compute the inverse of b (mod p)
   for (int i=0; i<AVX_ARRAY_SIZE; i++)
   {      
      dps[i] = (double) ps[i];
      multinvs[i] = (double) ComputeMultiplicativeInverse(ii_Base, ps[i]);
      powers[1][i] = (double) ii_Base;
   }
   
   avx_compute_reciprocal(dps, reciprocals);
   
   avx_powmod(multinvs, ii_Terms[0], dps, reciprocals);

   CheckAVXResult(ii_Terms[0], ps, dps);
   
   // Multiply successive terms by a (mod p)
   for (uint32_t idx=2; idx<MAX_POWERS+1; idx++)
   {
      // If the bsae is odd, then all n must be even and thus the difference
      // between any two remaining n for the base must also be even.
      if (idx > 2 && ii_Base & 1)
      {
         if (idx & 1)
            continue;
         
         avx_set_16a(powers[idx-2]);
         avx_set_16b(powers[2]);
         avx_mulmod(dps, reciprocals);
         avx_get_16a(powers[idx]);
      }
      else 
      {
         avx_set_16a(powers[idx-1]);
         avx_set_16b(powers[1]);
         avx_mulmod(dps, reciprocals);
         avx_get_16a(powers[idx]);
      }
   }
   
   avx_set_16a(multinvs);
 
   prevN = ii_Terms[0];
   
   // Note that the terms start at max N and decrease
   termIndex = 1;
   while (ii_Terms[termIndex] > 0)
   {
      theN = ii_Terms[termIndex];

      power = prevN - theN;

      // Our table isn't infinite in size, so we'll break down the power into smaller
      // chunks.  At worst this might cause an extra mulmod or two every once in a while.
      while (power > MAX_POWERS)
      {
         avx_set_16b(powers[MAX_POWERS]);
         avx_mulmod(dps, reciprocals);
         power -= MAX_POWERS;
      }

      avx_set_16b(powers[power]);
      avx_mulmod(dps, reciprocals);
      
      // At this point we have computed (1/b)^n (mod p).
      // If (1/b)^n (mod p) == n then we have a Woodall factor.
      // If (1/b)^n (mod p) == thePrime - n then we have a Cullen factor.
      CheckAVXResult(theN, ps, dps);
          
      prevN = theN;
      termIndex++;
   };
}

void  CullenWoodallWorker::CheckAVXResult(uint32_t theN, uint64_t *ps, double *dps)
{
   uint32_t idx;
   double __attribute__((aligned(32))) comparator[1]; 
   double __attribute__((aligned(32))) rems[AVX_ARRAY_SIZE];
   
   comparator[0] = (double) theN;
   
   // Only go further if one or more of the 16 primes yielded a factor for this n
   if (avx_pos_compare_1v(comparator) > 0)
   {
      avx_get_16a(rems);
      
      for (idx=0; idx<AVX_ARRAY_SIZE; idx++)
         if (rems[idx] == comparator[0])
         {
            if (ip_CullenWoodallApp->ReportFactor(ps[idx], theN, -1))
               VerifyFactor(ps[idx], theN, -1);
         }
   }
      
   // Only go further if one or more of the 16 primes yielded a factor for this n
   if (avx_neg_compare_1v(comparator, dps))
   {
      avx_get_16a(rems);
      
      for (idx=0; idx<AVX_ARRAY_SIZE; idx++)
         if (rems[idx] == dps[idx] - comparator[0])
         {
            if (ip_CullenWoodallApp->ReportFactor(ps[idx], theN, +1))
               VerifyFactor(ps[idx], theN, +1);
         }
   }
}
   
uint64_t CullenWoodallWorker::ComputeMultiplicativeInverse(uint64_t a, uint64_t p)
{
   uint64_t m = p;
   int64_t  y = 0, x = 1;
   uint64_t q, t;

   while (a > 1)
   {
      // q is quotient
      q = a / m;
      t = m;

      // m is remainder now, process same as
      // Euclid's algo
      m = a % m;
      a = t;
      t = y;

      // Update y and x
      y = x - q * y;
      x = t;
   }

   // Make x positive
   if (x < 0)
      x += p;

   return x;
}

// Build a list of powers for a from a^0 thru a^n for all even n up to count.
void  CullenWoodallWorker::BuildListOfPowers(uint64_t a, uint64_t p, uint32_t count, uint64_t *powers)
{
   uint32_t index;
   
   fpu_push_adivb(a, p);
   
   powers[0] = 1;
   powers[1] = a;
   
   // Multiply successive terms by a (mod p)
   for (index=2; index<count; index++)
       powers[index] = fpu_mulmod_iter(powers[index-1], a, p);

   fpu_pop();
}

void  CullenWoodallWorker::VerifyFactor(uint64_t p, uint32_t n, int32_t c)
{
   uint64_t rem;
      
   fpu_push_1divp(p);
   
   rem = fpu_powmod(ii_Base, n, p);
   rem = fpu_mulmod(rem, n, p);
   
   fpu_pop();
   
   if (c == -1 && rem != +1)
      FatalError("%" PRIu64" does not divide %u*%u^%u-1", p, n, ii_Base, n);
   
   if (c == +1 && rem != p-1)
      FatalError("%" PRIu64" does not divide %u*%u^%u+1", p, n, ii_Base, n);
}
