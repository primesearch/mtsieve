/* AlgebraicFactorHelper.cpp -- (C) Mark Rodenkirch, January 2018

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <inttypes.h>
#include <memory.h>
#include <time.h>
#include <stdarg.h>
#include "AlgebraicFactorHelper.h"

#include "../core/inline.h"
#include "../sieve/primesieve.hpp"

#define SMALL_PRIME_COUNT  100000

AlgebraicFactorHelper::AlgebraicFactorHelper(App *theApp, uint32_t base, uint32_t minN, uint32_t maxN)
{
   // Generate primes for this worker
   primesieve::generate_n_primes(SMALL_PRIME_COUNT, 1, &iv_SmallPrimes);

   ip_App = theApp;
   
   ip_AlgebraicFactorFile = 0;
   
   ii_Base = base;
   ii_MinN = minN;
   ii_MaxN = maxN;
   
   GetRoot(ii_Base, &ii_BRoot, &ii_BPower);
}

AlgebraicFactorHelper::~AlgebraicFactorHelper(void)
{
}

uint64_t AlgebraicFactorHelper::RemoveTermsWithAlgebraicFactors(seq_t *seqPtr)
{
   uint64_t  startingCount, removedCount = 0;
   
   if (seqPtr->c > 1 || seqPtr->c < -1)
      return 0;
   
   GetRoot(seqPtr->k, &ii_KRoot, &ii_KPower);
   
   startingCount = 1 + ii_MaxN - ii_MinN;

   CheckForSpecialForm(seqPtr);

   removedCount = RemoveSimpleTerm(seqPtr);
   
   if (startingCount - removedCount > 0)
      removedCount += RemoveTermsWithKPowers(seqPtr);
   
   if (startingCount - removedCount > 0)
      removedCount += RemoveTermsWithKAndBPowers(seqPtr);
   
   if (startingCount - removedCount > 0)
      removedCount += RemoveComplexRoot(seqPtr);
         
   if (startingCount - removedCount > 0)
      removedCount += CheckPower4(seqPtr);

   if (startingCount - removedCount > 0)
      removedCount += CheckBase2(seqPtr);
   
   if (ip_AlgebraicFactorFile != 0)
   {
      fclose(ip_AlgebraicFactorFile);
      ip_AlgebraicFactorFile = 0;
   }
   
   return removedCount;
}

// If c = -1 and k=2^f and b=2^g for any f and g, then this is a Mersenne number.
// If c = 1 and k=x^f and b=x^g then we have a Generalized Fermat Number.
// 
// Note that this will only give a warning, but will not remove terms since it
// is not looking for algebraic factors.
void    AlgebraicFactorHelper::CheckForSpecialForm(seq_t *seqPtr)
{
   // Needs to be k*b^n-1 or k*b^n+1
   if (seqPtr->c != -1 && seqPtr->c != 1)
      return;
   
   // If c = -1, then b must be 2^g for some g
   if (seqPtr->c == -1 && ii_BRoot != 2)
      return;

   if (seqPtr->k > 1)
   {      
      // If c = -1, then k must be 2^f for some f
      if (seqPtr->c == -1 && ii_KRoot != 2)
         return;
      
      // If c = +1, then k and b must have the same root
      if (seqPtr->c == +1 && ii_KRoot != ii_BRoot)
         return;
   }
   else
      ii_KPower = 0;

   if (seqPtr->c == -1)
      ip_App->WriteToConsole(COT_OTHER, "(sf) Sequence %" PRIu64"*%u^n-1 is the form of a Mersenne number", seqPtr->k, ii_Base);
   else
   {
      if (ii_KPower == 0)
      {
         if (ii_BPower == 1)
            ip_App->WriteToConsole(COT_OTHER, "(sf) Sequence %" PRIu64"*%u^n+1 as it is a GFN", seqPtr->k, ii_Base);
         else
            ip_App->WriteToConsole(COT_OTHER, "(sf) Sequence %" PRIu64"*%u^n+1 as it is a GFN --> %u^(%u*n)+1", seqPtr->k, ii_Base, ii_BRoot, ii_BPower);
      }
      else
      {
         if (ii_BPower == 1)
            ip_App->WriteToConsole(COT_OTHER, "(sf) Sequence %" PRIu64"*%u^n+1 as it is a GFN --> %u^(n+%u)+1", seqPtr->k, ii_Base, ii_BRoot, ii_KPower);
         else
            ip_App->WriteToConsole(COT_OTHER, "(sf) Sequence %" PRIu64"*%u^n+1 as it is a GFN --> %u^(%u*n+%u)+1", seqPtr->k, ii_Base, ii_BRoot, ii_BPower, ii_KPower);
      }
   }
}

// If k=x^f and b=x^g then we have a special form
//
// If c=-1, then x-1 will factor all terms
// If c=+1, then x+1 will factor terms where (f+n*g) is odd
uint64_t  AlgebraicFactorHelper::RemoveSimpleTerm(seq_t *seqPtr)
{
   uint32_t  removedCount = 0;
   uint32_t  n;
      
   // For ii_BRoot = 2 and seqPtr->c = -1, we have 1 as a divisor, ignore it.
   if (ii_BRoot + seqPtr->c == 1)
      return 0;

   if (seqPtr->k == 1)
   {
      ii_KPower = 0;
      
      ip_App->WriteToConsole(COT_OTHER, "(st) Sequence has algebraic factorization: %" PRIu64"*%u^n%+" PRId64" -> b = %u^%u", seqPtr->k, ii_Base, seqPtr->c,
               ii_BRoot, ii_BPower);
   }
   else 
   {
      if (ii_KRoot != ii_BRoot)
         return 0;
      
      ip_App->WriteToConsole(COT_OTHER, "(st) Sequence has algebraic factorization: %" PRIu64"*%u^n%+" PRId64" -> k = %u^%u and b = %u^%u", seqPtr->k, ii_Base, seqPtr->c,
               ii_KRoot, ii_KPower, ii_BRoot, ii_BPower);
   }   

   for (n=ii_MinN; n<=ii_MaxN; n++)
   {
      // For c = +1, the simple factor will only divide the term if ii_KPower + n*ii_BPower is odd
      if (seqPtr->c == +1)
         if ((ii_KPower + n*ii_BPower) % 2 == 0)
            continue;
      
      removedCount += CheckAndLogAlgebraicFactor(seqPtr, n, "%" PRId64"", ii_BRoot + seqPtr->c);
   }
   
   ip_App->WriteToConsole(COT_OTHER, "(st) Sequence %" PRIu64"*%u^n%+" PRId64" has %d terms removed has they have the factor %" PRId64"", 
         seqPtr->k, ii_Base, seqPtr->c, removedCount, ii_BRoot + seqPtr->c);

   return removedCount;
}

// If k=x^f then look for algebraic factors.
// Note that all n will be covered by these factors.
//
// If c=-1, z divides f and n*g, any n*g, then x^(f/z)*y^((n*g)/z)-1 will be a factor.
// If c=+1, z divides f and n*g, odd n*g, then x^(f/z)*y^((n*g)/z)+1 will be a factor.
uint64_t  AlgebraicFactorHelper::RemoveTermsWithKPowers(seq_t *seqPtr)
{
   uint32_t  removedCount = 0, loopRemovedCount = 0;
   uint32_t  n, idx;
   uint32_t  curroot, curpower;
   
   // If k = 1, then RemoveAlgebraicSimpleTerm() will have handled it
   if (seqPtr->k == 1)
      return 0;

   // If ii_KRoot == ii_BRoot, then RemoveAlgebraicSimpleTerm() will have handled it
   if (ii_KRoot == ii_BRoot)
      return 0;

   // seqPtr->k must be x^f for some x and f must be greater than 1
   if (ii_KPower == 1)
      return 0;


   
   for (idx=1; idx<ii_KPower; idx++)
   {
      // Given k=ii_KRoot^ii_KPower, find all idx where ii_KPower%idx = 0.
      // If ii_KPower == 6, then we look for algebraic factors with the
      // forms ii_KRoot^(6/2), ii_KRoot^(6/2), and ii_KRoot^(6/6).
      if (ii_KPower % idx != 0)
         continue;
      
      curpower = ii_KPower / idx;
      curroot = ii_KRoot;
      for (n=1; n<idx; n++)
         curroot *= ii_KRoot;
      
      if (seqPtr->c == +1)
      {
         // x^1 is not a divisor of x^n+1, so skip this idx
         if (idx == ii_KPower)
            continue;
         
         // x^y for even y is not a divisor of x^(y*n)+1, so skip this idx.
         if (curpower%2 == 0)
            continue;
      }

      // Since k = kroot^kpower we can remove all terms n where n%kpower = 0
      // For example 27*10^n+1 --> (3^3)*10^n+1 has factors in the form 3*10^n/3+1 
      ip_App->WriteToConsole(COT_OTHER, "(kp) Sequence has algebraic factorization: %" PRIu64"*%u^n%+" PRId64" -> (%u^%u)*%u^n%+" PRId64"", seqPtr->k, ii_Base, seqPtr->c,
               curroot, curpower, ii_Base, seqPtr->c);
               
      loopRemovedCount = 0;

      for (n=ii_MinN; n<=ii_MaxN; n++)
      {
         if (n % curpower != 0)
            continue;
         
         loopRemovedCount += CheckAndLogAlgebraicFactor(seqPtr, n, "(%u*%u^%u%+" PRId64")", curroot, ii_Base, n/curpower, seqPtr->c);
      }

      ip_App->WriteToConsole(COT_OTHER, "(kp) Sequence %" PRIu64"*%u^n%+" PRId64" has %d terms removed due to algebraic factors of the form %u*%u^(n/%u)%+" PRId64"", 
               seqPtr->k, ii_Base, seqPtr->c, loopRemovedCount, curroot, ii_Base, curpower, seqPtr->c);

      removedCount += loopRemovedCount;
   }

   return removedCount;
}

// If k=x^f and b=y^g then look for algebraic factors.
// Note that all n will be covered by these factors.
//
// If c=-1, z divides f and n*g, any n*g, then x^(f/z)*y^((n*g)/z)-1 will be a factor.
// If c=+1, z divides f and n*g, odd n*g, then x^(f/z)*y^((n*g)/z)+1 will be a factor.
uint64_t  AlgebraicFactorHelper::RemoveTermsWithKAndBPowers(seq_t *seqPtr)
{
   uint32_t  removedCount = 0, loopRemovedCount = 0;
   uint32_t  n, idx;
   uint32_t  curroot, curpower;
   
   // If k = 1, then RemoveAlgebraicSimpleTerm() will have handled it
   if (seqPtr->k == 1)
      return 0;

   // If ii_KRoot == ii_BRoot, then RemoveAlgebraicSimpleTerm() will have handled it
   if (ii_KRoot == ii_BRoot)
      return 0;

   // seqPtr->k must be x^f for some x and f must be greater than 1
   if (ii_KPower == 1)
      return 0;

   if (ii_BPower == 1)
      return removedCount;   
   
   for (idx=2; idx<=ii_KPower; idx++)
   {
      // Given k=ii_KRoot^ii_KPower, find all idx where ii_KPower%idx = 0.
      // If ii_KPower == 6, then we look for algebraic factors with the
      // forms ii_KRoot^(6/2), ii_KRoot^(6/2), and ii_KRoot^(6/6).
      if (ii_KPower % idx != 0)
         continue;
      
      curpower = ii_KPower / idx;
      curroot = ii_KRoot;
      for (n=1; n<idx; n++)
         curroot *= ii_KRoot;
      
      if (seqPtr->c == +1)
      {
         // x^1 is not a divisor of x^n+1, so skip this idx
         if (idx == ii_KPower)
            continue;
         
         // x^y for even y is not a divisor of x^(y*n)+1, so skip this idx.
         if (curpower%2 == 0)
            continue;
      }

      ip_App->WriteToConsole(COT_OTHER, "(kbp) Sequence has algebraic factorization: %" PRIu64"*%u^n%+" PRId64" -> (%u^%u)*%u^(%u*n)%+" PRId64"", 
                  seqPtr->k, ii_Base, seqPtr->c, curroot, curpower, ii_BRoot, ii_BPower, seqPtr->c);
               
      loopRemovedCount = 0;

      for (n=ii_MinN; n<=ii_MaxN; n++)
      {
         if ((n*ii_BPower) % curpower != 0)
            continue;
         
         // For c = +1, x^y does not divide x^(y*n)+1 when y*n is even
         if (seqPtr->c == +1)
            if (((n*ii_BPower) / curpower) % 2 == 0)
               continue;
         
         loopRemovedCount += CheckAndLogAlgebraicFactor(seqPtr, n, "(%u*%u^%u%+" PRId64")", curroot, ii_BRoot, (ii_BPower*n)/curpower, seqPtr->c);
      }

      if (idx == ii_KPower)
        ip_App->WriteToConsole(COT_OTHER, "(kbp) Sequence %" PRIu64"*%u^n%+" PRId64" has %u terms removed due to algebraic factors of the form %u*%u^((%u*n)/%u)%+" PRId64"", 
                  seqPtr->k, ii_Base, seqPtr->c, loopRemovedCount, curroot, ii_BRoot, ii_BPower, curpower, seqPtr->c);
      else
        ip_App->WriteToConsole(COT_OTHER, "(kbp) Sequence %" PRIu64"*%u^n%+" PRId64" has %u terms removed due to algebraic factors of the form (%u^%u)*%u^((%u*n)/%u)%+" PRId64"", 
                  seqPtr->k, ii_Base, seqPtr->c, loopRemovedCount, curroot, curpower, ii_BRoot, ii_BPower, curpower, seqPtr->c);

      removedCount += loopRemovedCount;
   }

   return removedCount;
}
   
// If seqPtr->k*b^ii = root^m for root > 1 and m > 1, then we can remove algebraic factors
// because we can find a simple root for some of the terms.
//
//   1)  base^n = b^(bpow*n)
//
//   2)  seqPtr->k -> k*b^y
// 
//   3)  seqPtr->k*base^n -> (k*b^y)*(b^(bpow*n))
//                    -> k*b^(bpow*n+y)
//
//   4)  For each ii where k*b^ii = root^m for root > 1 and m > 1:
//                    -> (k*b^ii)*b^(bpow*n+y-ii)
//                    -> root^m*b^(bpow*n+y-ii)
//
//   5)  Find prime factors of m -> [f1,f2,...,fn]
//
//   6)  For each prime factor f of m, if (bpow*n+y-ii)%f = 0, then:
//          if c = +1 and f is odd: root^f*b^(bpow*n+y-ii)%f+1 has a factor of root*b^(bpow*n+y-ii)/f+1
//          if c = -1: root^f*b^(bpow*n+y-ii)%f-1 has a factor of root*b^(bpow*n+y-ii)/f-1
uint64_t  AlgebraicFactorHelper::RemoveComplexRoot(seq_t *seqPtr)
{
   uint32_t  removedCount = 0, loopRemovedCount;
   uint32_t  ii, m, n, z, z1, z2, idx;
   uint32_t  kf_count, bf_count, kbf_count;
   uint32_t  bf_factor[50], bf_power[50];
   uint32_t  kf_factor[50], kf_power[50];
   uint32_t  kbf_factor[50], kbf_power[50];
   char      root[50], part[50], addParenthesis;

   // We want seqPtr->c = +1 or -1
   if (seqPtr->c != 1 && seqPtr->c != -1)
      return 0;

   // base and seqPtr->k must share a divisor
   if (gcd64(ii_Base, seqPtr->k) == 1)
      return 0;
   
   // If k is r^x and b is r^y for some r, x, and y, then we have simple roots that
   // are handled elsewhere.
   if (ii_BRoot == ii_KRoot)
      return 0;
   
   // Step 1:  Factorize b
   bf_count = GetFactorList(ii_Base, bf_factor, bf_power);
   
   // Step 2:  Factorize k
   kf_count = GetFactorList(seqPtr->k, kf_factor, kf_power);
      
   for (z1=0; z1<kf_count; z1++)
   {
      kbf_factor[z1] = kf_factor[z1];
      kbf_power[z1] = kf_power[z1];
   }
   
   kbf_count = kf_count;
      
   for (ii=1; ii<10; ii++)
   {
      // We could multiply out seqPtr->k*base^ii and factor that value, but
      // seqPtr->k*base^ii will likely exceed what can be stored in a 64-bit integer.
      // We'll just mulitply k*b^(ii-1) by b
      for (z2=0; z2<bf_count; z2++)
      {
         idx = 99999;

         for (z1=0; z1<kbf_count; z1++)
            if (kbf_factor[z1] == bf_factor[z2])
            {
               kbf_power[z1] += bf_power[z2];
               idx = z1;
            }
         
         if (idx == 99999)
         {
            kbf_factor[kbf_count] = bf_factor[z2];
            kbf_power[kbf_count] = bf_power[z2];
            kbf_count++;            
         }
      }

      // We now have seqPtr->k*base^ii = f1^p1 * f2^p2 ... fn^pn
      // Step 4:  Determine m such that m = gcd(p1, p2, ..., pn)
      m = kbf_power[0];

      for (idx=1; idx<kbf_count; idx++)
         m = gcd32(m, kbf_power[idx]);
      
      // If m then seqPtr->k*base^ii is not a perfect power
      if (m == 1)
         continue;

      for (idx=2; idx<m; idx++)
      {         
         // Given k*b^ii=root^m, find all idx where m%idx = 0.
         // If m == 6, then we look for algebraic factors with the
         // forms root^2 and root^3.
         if (m % idx != 0)
            continue;
         
         z = m / idx;
                  
         root[0] = 0;

         for (z1=0; z1<kbf_count; z1++)
         {
            if (kbf_power[z1] == 0)
               continue;
            
            if (z1 == 0)
            {
               if (kbf_power[z1] == z)
                  sprintf(root, "%u", kbf_factor[z1]);
               else
               {
                  addParenthesis = 1;
                  sprintf(root, "%u^%u", kbf_factor[z1], kbf_power[z1] / z);
               }
            }
            else
            {
               addParenthesis = 1;
               sprintf(part, "%s", root);
               
               if (kbf_power[z1] == z)
                  sprintf(root, "%s*%u", part, kbf_factor[z1]);
               else
                  sprintf(root, "%s*%u^%u", part, kbf_factor[z1], kbf_power[z1] / z);
            }
         }
         
         if (addParenthesis)
         {
            if (ii == 1)
              ip_App->WriteToConsole(COT_OTHER, "(cr) Sequence has algebraic factorization: %" PRIu64"*%u^n%+" PRId64" -> %" PRIu64"*%u = (%s)^%u",
                                     seqPtr->k, ii_Base, seqPtr->c, seqPtr->k, ii_Base, root, z);
            else
              ip_App->WriteToConsole(COT_OTHER, "(cr) Sequence has algebraic factorization: %" PRIu64"*%u^n%+" PRId64" -> %" PRIu64"*%u^%u = (%s)^%u",
                                     seqPtr->k, ii_Base, seqPtr->c, seqPtr->k, ii_Base, ii, root, z);
         }
         else
         {
            if (ii == 1)
              ip_App->WriteToConsole(COT_OTHER, "(cr) Sequence has algebraic factorization: %" PRIu64"*%u^n%+" PRId64" -> %" PRIu64"*%u = %s^%u",
                                     seqPtr->k, ii_Base, seqPtr->c, seqPtr->k, ii_Base, root, z);
            else
              ip_App->WriteToConsole(COT_OTHER, "(cr) Sequence has algebraic factorization: %" PRIu64"*%u^n%+" PRId64" -> %" PRIu64"*%u^%u = %s^%u",
                                     seqPtr->k, ii_Base, seqPtr->c, seqPtr->k, ii_Base, ii, root, z);
         }

         // Don't allow n to go negative
         if (ii_MinN < ii)
            n = ii_MinN;
         else
            n = ii_MinN - ii;

         // n-ii must be greater than 0
         if (n < ii+1) n = ii + 1;
         
         // To avoid 2^1-1 divisors
         if (!strcmp(root, "2") && seqPtr->c == -1 && n-ii==1) n++;
         
         // For c = +1, if z is odd, then double it to make it even as we can only
         // provide factors for odd n.
         if (seqPtr->c == +1 && z%2 == 1)
            z *= 2;

         // Find smallest n >= ii_MinN where (n-ii)%z = 0
         while ((n - ii) % z != 0) n++;
         
         // For c = +1, n must be odd in order for root^(n/f)+1 to be a factor
         if (seqPtr->c == +1 && (n - ii)%2 == 0)
         {
            // If z and n are even, then n%z will never be 0 so we can skip this z.
            if (z%2 == 0)
               continue;
            
            n += z;
         }
                     
         loopRemovedCount = 0;

         for (; n<=ii_MaxN; n+=z)
            if (addParenthesis)
               loopRemovedCount += CheckAndLogAlgebraicFactor(seqPtr, n, "((%s)*%u^%u%+" PRId64")", root, ii_Base, (n-ii)/z, seqPtr->c);
            else
               loopRemovedCount += CheckAndLogAlgebraicFactor(seqPtr, n, "(%s*%u^%u%+" PRId64")", root, ii_Base, (n-ii)/z, seqPtr->c);
         
         if (addParenthesis)
           ip_App->WriteToConsole(COT_OTHER, "(cr) Sequence %" PRIu64"*%u^n%+" PRId64" has %u terms removed due to algebraic factors of the form (%s)*%u^((n-%d)/%d)%+" PRId64"", 
                     seqPtr->k, ii_Base, seqPtr->c, loopRemovedCount, root, ii_Base, ii, z, seqPtr->c);
         else
           ip_App->WriteToConsole(COT_OTHER, "(cr) Sequence %" PRIu64"*%u^n%+" PRId64" has %u terms removed due to algebraic factors of the form %s*%u^((n-%d)/%d)%+" PRId64"", 
                     seqPtr->k, ii_Base, seqPtr->c, loopRemovedCount, root, ii_Base, ii, z, seqPtr->c);
         
         removedCount += loopRemovedCount;
      }     
   }

   return removedCount;
}

uint64_t  AlgebraicFactorHelper::CheckBase2(seq_t *seqPtr)
{
   uint64_t  k;
   uint32_t  removedCount = 0;
   uint32_t  m, n, root, y, idx;
   uint32_t  kf_count;
   uint32_t  b, bpow;
   uint32_t  kf_factor[50], kf_power[50];

   // We want c = 1
   if (seqPtr->c != 1)
      return 0;
   
   bpow = 0;
   b = ii_Base;
   while (b > 1)
   {
      // if b is odd and greater than 1, then we can't use it
      if (b & 1)
         return 0;
      
      b >>= 1;
      bpow++;
   }

   y = 0;
   k = seqPtr->k;
   while (k % 2 == 0)
   {
      k >>= 1;
      y++;
   }

   // To get rid of compiler warnings
   kf_power[0] = 0;

   // Now that the base is removed, refactorize k
   kf_count = GetFactorList(k, kf_factor, kf_power);
   
   // We now have k = f1^p1 * f2^p2 ... fn^pn
   // Now determine m such that m = gcd(p1, p2, ..., pn)
   m = kf_power[0];
      
   for (idx=1; idx<kf_count; idx++)
      m = gcd32(m, kf_power[idx]);

   // We want k where k=root^(4*m)
   if (m % 4 != 0)
      return 0;

   root = 1;
   for (idx=0; idx<kf_count; idx++)
      root *= (uint32_t) pow((double) kf_factor[idx], (double) (kf_power[idx] / 4));

   n = ii_MinN;
   
   // Exit the function right away if n*bpow+y is not divisible
   // by 4 and the gcd of bpow and 4 is not 1 because we'll get
   // an infinite loop in the while statement right after this.
   if ((((n*bpow+y) % 4) != 2) && gcd32(bpow, 4) > 1)
      return 0;
         
   while ((n*bpow+y)% 4 != 2) n++;
   while (n > ii_MinN && n > 4) n -= 4;
   while (n < ii_MinN) n += 4;

   for ( ; n<=ii_MaxN; n+=4)
      removedCount += CheckAndLogAlgebraicFactor(seqPtr, n, "(%u^2*2^%u-%u*2^%u+1)",
            root, (n*bpow+y)/2, root, (n*bpow+y+2)/4);

   if (removedCount > 0)
   {
      ip_App->WriteToConsole(COT_OTHER, "(b2) Removed %u algebraic factors for %" PRIu64"*%u^n+1 of the form (%u^2)*2^(n/2)-%u*2^((n+2)/4))+1 when n%%4=2",
               removedCount, seqPtr->k, ii_Base, root, root);
   }

   return removedCount;
}

// If k = 4*x^4 and n%4 = 0 then k*b^n+1 has 2*x^2*b^n/2-2*x*b^n/4+1 and 2*x^2*b^n/2+2*x*b^n/4+1  as factors
uint64_t  AlgebraicFactorHelper::CheckPower4(seq_t *seqPtr)
{
   uint32_t  removedCount = 0;
   uint32_t  n, ninc;
   uint32_t  b1, bexp;

   // We want c = +1
   if (seqPtr->c != 1)
      return 0;

   // We want k to be divisible by 4
   if (seqPtr->k % 4 != 0)
      return 0;
   
   if (ii_BPower % 4 == 0)
      ninc = 1;
   
   else if (ii_BPower % 4 == 2)
   {
      ninc = 2;
      
      if (ii_BPower > 2)
      {
         // If base = x^(y*2) then compute ii_BRoot as x^y
         // In other words set ii_BRoot = sqrt(base)
         bexp = ii_BPower / 2;
         
         b1 = 1;
         while (bexp > 0)
         {
            b1 *= ii_BRoot;
            bexp--;
         }
      
         ii_BPower = 2;
         ii_BRoot = b1;
      }
   }
   else
   {
      ii_BRoot = ii_Base;
      ninc = 4;
   }
   
   if (seqPtr->k == 4)
      ii_KRoot = 1;
   else
   {
      // Now we have seqPtr->k = 4*ii_KRoot^ii_KPower
      GetRoot(seqPtr->k/4, &ii_KRoot, &ii_KPower);
      
      // We want ii_KPower to be a multiple of 4
      if (ii_KPower % 4 != 0)
         return 0;
   }

   n = ii_MinN;
   while (n % ninc > 0) n++;
   
   for (; n<=ii_MaxN; n+=ninc)
   {
      if (ii_KRoot == 1)
      {
         if (ninc == 1)
            removedCount += CheckAndLogAlgebraicFactor(seqPtr, n, "(2*(%u^%u)^%u+2*(%u^%u)^%u+1)",
                  ii_BRoot, ii_BPower/2, n, ii_BRoot, ii_BPower/4, n);
         else if (ninc == 2)
            removedCount += CheckAndLogAlgebraicFactor(seqPtr, n, "(2*%u^%u+2*%u^%u+1)",
                  ii_BRoot, n, ii_BRoot, n/2);
         else
            removedCount += CheckAndLogAlgebraicFactor(seqPtr, n, "(2*%u^%u+2*%u^%u+1)",
                  ii_Base, n/2, ii_BRoot, n/4);
      }
      else
      {
         if (ninc == 1)
            removedCount += CheckAndLogAlgebraicFactor(seqPtr, n, "(2*%u^%u*(%u^%u)^%u+2*%u^%u*(%u^%u)^%u+1)",
                  ii_KRoot, ii_KPower/2, ii_BRoot, ii_BPower/2, n, ii_KRoot, ii_KPower/4, ii_BRoot, ii_BPower/4, n);
         else if (ninc == 2)
            removedCount += CheckAndLogAlgebraicFactor(seqPtr, n, "(2*%u^%u*%u^%u+2*%u^%u*%u^%u+1)",
                  ii_KRoot, ii_KPower/2, ii_BRoot, n, ii_KRoot, ii_KPower/4, ii_BRoot, n/2);
         else
            removedCount += CheckAndLogAlgebraicFactor(seqPtr, n, "(2*%u^%u*%u^%u+2*%u^%u*%u^%u+1)",
                  ii_KRoot, ii_KPower/2, ii_Base, n/2, ii_KRoot, ii_KPower/4, ii_Base, n/4);
      }
   }

   if (removedCount)
   {
      if (ii_KRoot == 1)
      {
         if (ninc == 1)
           ip_App->WriteToConsole(COT_OTHER, "(p4) Removed %u algebraic factors for %" PRIu64"*%u^n%+" PRId64" of the form (2*(%u^%u)^n+2*(%u^%u)^n+1)",
                  removedCount, seqPtr->k, ii_Base, seqPtr->c, ii_BRoot, ii_BPower/2, ii_BRoot, ii_BPower/4);
         else if (ninc == 2)
           ip_App->WriteToConsole(COT_OTHER, "(p4) Removed %d algebraic factors for %" PRIu64"*%u^n%+" PRId64" of the form (2*%u^%u*(%u^%u)^n+2*%u^%u*(%u^%u)^n+1)",
                  removedCount, seqPtr->k, ii_Base, seqPtr->c, ii_KRoot, ii_KPower/2, ii_BRoot, n, ii_KRoot, ii_KPower/4, ii_BRoot, n/2);

         else
           ip_App->WriteToConsole(COT_OTHER, "(p4) Removed %u algebraic factors for %" PRIu64"*%u^n%+" PRId64" of the form (2*%u^(n/2)+2*%u^(n/4)+1)",
                  removedCount, seqPtr->k, ii_Base, seqPtr->c, ii_BRoot, ii_BRoot);
      }
      else
      {
         if (ninc == 1)
           ip_App->WriteToConsole(COT_OTHER, "(p4) Removed %u algebraic factors for %" PRIu64"*%u^n%+" PRId64" of the form (2*%u^%u*(%u^%u)^n+2*%u^%u*(%u^%u)^n+1)",
                  removedCount, seqPtr->k, ii_Base, seqPtr->c, ii_KRoot, ii_KPower/2, ii_BRoot, ii_BPower/2, ii_KRoot, ii_KPower/4, ii_BRoot, ii_BPower/4);
         else if (ninc == 2)
           ip_App->WriteToConsole(COT_OTHER, "(p4) Removed %u algebraic factors for %" PRIu64"*%u^n%+" PRId64" of the form (2*%u^%u*(%u^%u)^n+2*%u^%u*(%u^%u)^n+1)",
                  removedCount, seqPtr->k, ii_Base, seqPtr->c, ii_KRoot, ii_KPower/2, ii_BRoot, n, ii_KRoot, ii_KPower/4, ii_BRoot, n/2);
         else
           ip_App->WriteToConsole(COT_OTHER, "(p4) Removed %u algebraic factors for %" PRIu64"*%u^n%+" PRId64" of the form (2*%u^%u*%u^(n/2)+2*%u^%u*%u^(n/4)+1)",
                  removedCount, seqPtr->k, ii_Base, seqPtr->c, ii_KRoot, ii_KPower/2, ii_BRoot, ii_KRoot, ii_KPower/4, ii_BRoot);
      }
   }

   return removedCount;
}

uint32_t   AlgebraicFactorHelper::CheckAndLogAlgebraicFactor(seq_t *seqPtr, uint32_t n, const char *fmt, ...)
{   
   va_list args;

   if (n < ii_MinN || n > ii_MaxN)
      return 0;

   if (!seqPtr->nTerms[n-ii_MinN])
      return 0;
   
   seqPtr->nTerms[n-ii_MinN] = false;

   if (!ip_AlgebraicFactorFile)
   {   
      char  fName[50];
      
      sprintf(fName, "alg_%" PRIu64"_%u_%+" PRId64".log", seqPtr->k, ii_Base, seqPtr->c);
      
      ip_AlgebraicFactorFile = fopen(fName, "a");
   }
   
   va_start(args,fmt);
   vfprintf(ip_AlgebraicFactorFile, fmt, args);
   va_end(args);
   
   fprintf(ip_AlgebraicFactorFile, " | %" PRIu64"*%u^%d%+" PRId64"\n", seqPtr->k, ii_Base, n, seqPtr->c);
   
   return 1;
}

uint32_t AlgebraicFactorHelper::GetFactorList(uint64_t theNumber, uint32_t *factorList, uint32_t *powerList)
{
   uint32_t  distinctFactors = 0;
   uint32_t  power;
   uint64_t  thePrime;

   vector<uint64_t>::iterator it = iv_SmallPrimes.begin(); 
         
   // Get the prime factorization of the input number.  Note that since the vector only has
   // 100000 primes, some numbers might not get fully factored.
   while (it != iv_SmallPrimes.end())
   {
      thePrime = *it;
      it++;
      
      power = 0;
      while (theNumber % thePrime == 0)
      {
         theNumber /= thePrime;
         power++;
      }

      if (power > 0)
      {
         factorList[distinctFactors] = thePrime;
         powerList[distinctFactors] = power;
         distinctFactors++;
         
         if (theNumber < thePrime * thePrime)
            break;
      }
   }

   // If then number wasn't fully factored, that's fine.
   if (theNumber > 1)
   {
      factorList[distinctFactors] = theNumber;
      powerList[distinctFactors] = 1;
      distinctFactors++;
   }

   return distinctFactors;
}

// If c = -1 and k=2^f and b=2^g for any f and g, then this is a Mersenne number.
// If c = 1 and k=x^f and b=x^g then we have a Generalized Fermat Number.
bool     AlgebraicFactorHelper::IsGfnOrMersenneForm(seq_t *seqPtr)
{   
   // Needs to be k*b^n-1 or k*b^n+1
   if (seqPtr->c != -1 && seqPtr->c != 1)
      return false;
   
   // If c = -1, then b must be 2^g for some g
   if (seqPtr->c == -1 && ii_BRoot != 2)
      return false;

   if (seqPtr->k > 1)
   {
      // If c = -1, then k must be 2^f for some f
      if (seqPtr->c == -1 && ii_KRoot != 2)
         return false;
      
      // If c = +1, then k and b must have the same root
      if (seqPtr->c == +1 && ii_BRoot != ii_KRoot)
         return false;
   }

   return true;
}

// Find root and power such that root^power = number
void  AlgebraicFactorHelper::GetRoot(uint64_t number, uint32_t *root, uint32_t *power)
{
   uint32_t  idx, r, rpow, rf_count;
   uint32_t  rf_factor[50], rf_power[50];

   // To get rid of compiler warnings
   rf_power[0] = 0;
   
   rf_count = GetFactorList(number, rf_factor, rf_power);
   
   rpow = rf_power[0];
   for (idx=1; idx<rf_count; idx++)
      rpow = gcd32(rpow, rf_power[idx]);
      
   r = 1;
   for (idx=0; idx<rf_count; idx++)
      r *= (uint32_t) pow((double) rf_factor[idx], (double) (rf_power[idx] / rpow));
   
   *root = r;
   *power = rpow;
}