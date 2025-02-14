/* CullenWoodallWorker.cpp -- (C) Mark Rodenkirch, May 2018

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <time.h>

#include "CullenWoodallWorker.h"

#ifdef USE_X86
#include "../x86_asm/fpu-asm-x86.h"
#include "../x86_asm/avx-asm-x86.h"
#endif

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
   ii_Terms = (uint32_t *) xmalloc(ii_MaxTermCount, sizeof(int32_t), "terms");

#ifdef USE_X86
   if (ip_CullenWoodallApp->UseAvxIfAvailable() && CpuSupportsAvx())
   {
      if (ii_Base < ii_MaxN)
         SetMiniChunkRange(ii_MaxN + 1, PMAX_MAX_52BIT, AVX_ARRAY_SIZE);
      else
         SetMiniChunkRange(ii_Base + 1, PMAX_MAX_52BIT, AVX_ARRAY_SIZE);
   }
#endif

   il_NextTermsBuild = 0;
   ib_Initialized = true;
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
   
   for (uint32_t pIdx=0; pIdx<ii_PrimesInList; pIdx+=4)
   {
      ps[0] = il_PrimeList[pIdx+0];
      ps[1] = il_PrimeList[pIdx+1];
      ps[2] = il_PrimeList[pIdx+2];
      ps[3] = il_PrimeList[pIdx+3];
      
      // Every once in a while rebuild the term lists as it will have fewer entries
      // which will speed up testing for the next range of p.
      // Unlike the GPU, we will put all terms into one group.
      if (ps[0] > il_NextTermsBuild)
      {
         ip_CullenWoodallApp->GetTerms(ii_Terms, ii_MaxTermCount, ii_MaxTermCount);
         
         il_NextTermsBuild = (ps[3] << 1);
      }
      
      if (ps[0] < maxPForSmallPrimeLogic)
         TestSmallPrimes(ps);
      else
      {
#ifdef USE_X86
         TestLargePrimesFPU(ps);
#else
         MpArithVec mp(ps);
         TestLargePrimes(ps, mp);
#endif
      }
      
      SetLargestPrimeTested(ps[3], 4);
      
      if (ps[3] > maxPrime)
         break;
   }
}

void  CullenWoodallWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
#ifdef USE_X86
   // Every once in a while rebuild the term lists as it will have fewer entries
   // which will speed up testing for the next range of p.
   // Unlike the GPU, we will put all terms into one group.
   if (miniPrimeChunk[0] > il_NextTermsBuild)
   {
      ip_CullenWoodallApp->GetTerms(ii_Terms, ii_MaxTermCount, ii_MaxTermCount);
      
      il_NextTermsBuild = (miniPrimeChunk[AVX_ARRAY_SIZE-1] << 1);
   }

   TestPrimesAVX(miniPrimeChunk);
#else
   FatalError("CullenWoodallWorker::TestMiniPrimeChunk not implemented");
#endif
}

// For small p, we need to iterate from minn to maxn as the large prime algorithm
// does not work correctly when p < maxN.
void  CullenWoodallWorker::TestSmallPrimes(uint64_t *ps)
{
   uint32_t theN, prevN;
   int32_t  termIndex, maxIndex;
   uint32_t power;
   uint64_t thePrime;
   
   MpRes    powers[MAX_POWERS+1];
   MpRes    res, resCW;
   
   if (ip_CullenWoodallApp->GetTermCount() == 0)
      return;
   
   maxIndex = 0;
   while (ii_Terms[maxIndex] > 0)
      maxIndex++;
      
   for (int i=0; i<4; i++)
   {
      thePrime = ps[i];
      
      if (ii_Base % thePrime == 0)
         continue;
      
      MpArith mp(thePrime);
      const MpRes resBase = mp.nToRes(ii_Base);
      const MpRes pOne = mp.one();
      const MpRes mOne = mp.sub(mp.zero(), pOne);
      
      termIndex = maxIndex - 1;

      prevN = theN = ii_Terms[termIndex];
         
      BuildListOfPowers(resBase, thePrime, mp, MAX_POWERS+1, powers);
      
      // Compute our starting term
      res = mp.pow(resBase, theN);
      
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
            res = mp.mul(res, powers[MAX_POWERS]);
            power -= MAX_POWERS;
         }

         res = mp.mul(res, powers[power]);

         resCW = mp.mul(res, mp.nToRes(theN));
         
         if (resCW == pOne)
            ip_CullenWoodallApp->ReportFactor(thePrime, theN, -1);

         if (resCW == mOne)
            ip_CullenWoodallApp->ReportFactor(thePrime, theN, +1);

         prevN = theN;
         termIndex--;
      };
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
void  CullenWoodallWorker::TestLargePrimes(uint64_t *ps, MpArithVec mp)
{
   uint32_t theN, prevN;
   uint32_t termIndex;
   uint64_t powInv[4], temp[4];
   uint32_t power;
   
   MpResVec powers[MAX_POWERS+1];
   MpResVec resC[MAX_POWERS+1];
   MpResVec resW[MAX_POWERS+1];
   MpResVec res, cullen, woodall;
   
   if (ip_CullenWoodallApp->GetTermCount() == 0)
      return;
   
   // compute the inverse of b (mod p)
   powInv[0] = InvMod32(ii_Base, ps[0]);
   powInv[1] = InvMod32(ii_Base, ps[1]);
   powInv[2] = InvMod32(ii_Base, ps[2]);
   powInv[3] = InvMod32(ii_Base, ps[3]);
   
   res = mp.pow(mp.nToRes(powInv), ii_Terms[0]);
   
   theN = ii_Terms[0];

   // res = b^minN % p

   temp[0] = ps[0] - theN;
   temp[1] = ps[1] - theN;
   temp[2] = ps[2] - theN;
   temp[3] = ps[3] - theN;
   
   cullen = mp.nToRes(temp);

   temp[0] = theN;
   temp[1] = theN;
   temp[2] = theN;
   temp[3] = theN;
   
   woodall = mp.nToRes(temp);
   
   for (size_t k = 0; k < VECTOR_SIZE; ++k)
   {
      // If res = p - miN, then we have a factor
      if (res[k] == cullen[k])
         ip_CullenWoodallApp->ReportFactor(ps[k], theN, +1);

      // If res = miN, then we have a factor
      if (res[k] == woodall[k])
         ip_CullenWoodallApp->ReportFactor(ps[k], theN, -1);
   }
   
   powers[1] = mp.nToRes(ii_Base);   
   powers[2] = mp.mul(powers[1], powers[1]);
   
   resC[1] = mp.sub(mp.nToRes(ps), mp.one());
   resW[1] = mp.one();
   
   resC[2] = mp.add(resC[1], resC[1]);
   resW[2] = mp.add(resW[1], resW[1]);
   
   if (ii_Base & 1)
   {
      // If the base is odd, then all n must be even and thus the difference
      // between any two remaining n for the base must also be even.
      for (size_t idx=4; idx<MAX_POWERS+1; idx+=2)
      {
         powers[idx] = mp.mul(powers[idx-2], powers[2]);
         
         resC[idx] = mp.add(resC[idx-2], resC[2]);
         resW[idx] = mp.add(resW[idx-2], resW[2]);
      }
   }
   else
   {
      for (size_t idx=3; idx<MAX_POWERS+1; idx+=2)
      {
         powers[idx] = mp.mul(powers[idx-1], powers[1]);
         
         resC[idx] = mp.add(resC[idx-1], resC[1]);
         resW[idx] = mp.add(resW[idx-1], resW[1]);
      }
   }
   
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
         res = mp.mul(res, powers[MAX_POWERS]);
         
         cullen = mp.sub(cullen, resC[MAX_POWERS]);
         woodall = mp.sub(woodall, resW[MAX_POWERS]);
         
         power -= MAX_POWERS;
      }

      res = mp.mul(res, powers[power]);

      cullen = mp.sub(cullen, resC[power]);
      woodall = mp.sub(woodall, resW[power]);

      // At this point we have computed (1/b)^n (mod p).
      // If (1/b)^n (mod p) == n then we have a Woodall factor.
      // If (1/b)^n (mod p) == thePrime - n then we have a Cullen factor.

      for (size_t k = 0; k < VECTOR_SIZE; ++k)
      {
         if (res[k] == cullen[k])
            ip_CullenWoodallApp->ReportFactor(ps[k], theN, +1);
            
         if (res[k] == woodall[k]) 
            ip_CullenWoodallApp->ReportFactor(ps[k], theN, -1);
      }
                  
      prevN = theN;
      termIndex++;
   };
}

#ifdef USE_X86
void  CullenWoodallWorker::TestLargePrimesFPU(uint64_t *ps)
{
   uint32_t theN, prevN;
   uint32_t termIndex;
   uint64_t powinvs[4];
   uint64_t thePrime;
   uint64_t powers[MAX_POWERS+1][4];
   uint64_t rems[4];
   uint32_t power;
   
   if (ip_CullenWoodallApp->GetTermCount() == 0)
      return;
   
   // compute the inverse of b (mod p)
   powinvs[0] = InvMod32(ii_Base, ps[0]);
   powinvs[1] = InvMod32(ii_Base, ps[1]);
   powinvs[2] = InvMod32(ii_Base, ps[2]);
   powinvs[3] = InvMod32(ii_Base, ps[3]);

   fpu_powmod_4b_1n_4p(powinvs, ii_Terms[0], ps);
   
   for (int i=0; i<4; i++)
   {
      thePrime = ps[i];
      
      theN = ii_Terms[0];
      
      if (powinvs[i] == theN)
         ip_CullenWoodallApp->ReportFactor(thePrime, theN, -1);
         
      if (powinvs[i] == thePrime - theN)
         ip_CullenWoodallApp->ReportFactor(thePrime, theN, +1);
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
         ip_CullenWoodallApp->ReportFactor(ps[0], theN, -1);
         
      if (rems[0] == ps[0] - theN)
         ip_CullenWoodallApp->ReportFactor(ps[0], theN, +1);
         
      if (rems[1] == theN)
         ip_CullenWoodallApp->ReportFactor(ps[1], theN, -1);
         
      if (rems[1] == ps[1] - theN)
         ip_CullenWoodallApp->ReportFactor(ps[1], theN, +1);
         
      if (rems[2] == theN)
         ip_CullenWoodallApp->ReportFactor(ps[2], theN, -1);
         
      if (rems[2] == ps[2] - theN)
         ip_CullenWoodallApp->ReportFactor(ps[2], theN, +1);
         
      if (rems[3] == theN)
         ip_CullenWoodallApp->ReportFactor(ps[3], theN, -1);
         
      if (rems[3] == ps[3] - theN)
         ip_CullenWoodallApp->ReportFactor(ps[3], theN, +1);
         
      prevN = theN;
      termIndex++;
   };
   
   fpu_pop();
   fpu_pop();
   fpu_pop();
   fpu_pop();
}

// Same as TestLargePrimes, but using AVX
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
      multinvs[i] = (double) InvMod32(ii_Base, ps[i]);
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
            ip_CullenWoodallApp->ReportFactor(ps[idx], theN, -1);
   }
      
   // Only go further if one or more of the 16 primes yielded a factor for this n
   if (avx_neg_compare_1v(comparator, dps))
   {
      avx_get_16a(rems);
      
      for (idx=0; idx<AVX_ARRAY_SIZE; idx++)
         if (rems[idx] == dps[idx] - comparator[0])
            ip_CullenWoodallApp->ReportFactor(ps[idx], theN, +1);
   }
}
#endif

// Build a list of powers for a from a^0 thru a^n for all even n up to count.
void  CullenWoodallWorker::BuildListOfPowers(MpRes a, uint64_t p, MpArith mp, uint32_t count, MpRes *powers)
{
   uint32_t index;

   powers[0] = mp.one();
   powers[1] = a;
   
   // Multiply successive terms by a (mod p)
   for (index=2; index<count; index++)
       powers[index] = mp.mul(powers[index-1], a);
}
