/* AlgebraicFactorApp.cpp -- (C) Mark Rodenkirch, January 2018

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
#include "inline.h"
#include "Clock.h"
#include "AlgebraicFactorApp.h"

#include "../sieve/primesieve.hpp"

#define SMALL_PRIME_COUNT  100000

AlgebraicFactorApp::AlgebraicFactorApp(void)
{
   // Generate primes for this worker
   primesieve::generate_n_primes(SMALL_PRIME_COUNT, 1, &iv_SmallPrimes);
}

AlgebraicFactorApp::~AlgebraicFactorApp(void)
{

}

// Find root and power such that root^power = number
void  AlgebraicFactorApp::GetRoot(uint64_t theNumber, uint64_t *root, uint32_t *power)
{
   uint32_t  idx, rpow, rf_count;
   uint64_t  r, rf_factor[50];
   uint32_t  rf_power[50];

   for (idx=0; idx<50; idx++)
      rf_power[idx] = 0;

   rf_count = GetFactorList(theNumber, rf_factor, rf_power);
   
   rpow = rf_power[0];
   for (idx=1; idx<rf_count; idx++)
      rpow = gcd32(rpow, rf_power[idx]);
      
   r = 1;
   for (idx=0; idx<rf_count; idx++)
      r *= (uint32_t) pow((double) rf_factor[idx], (double) (rf_power[idx] / rpow));
   
   *root = r;
   *power = rpow;
}

uint32_t AlgebraicFactorApp::GetFactorList(uint64_t theNumber, uint64_t *factorList, uint32_t *powerList)
{
   uint32_t  distinctFactors = 0;
   uint32_t  power;
   uint64_t  thePrime;

   std::vector<uint64_t>::iterator it = iv_SmallPrimes.begin(); 
         
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
bool     AlgebraicFactorApp::IsGfnOrMersenneForm(uint64_t k, uint32_t base, int32_t c)
{
   uint64_t  broot, kroot;
   uint32_t  bpower, kpower;
   
   GetRoot(base, &broot, &bpower);
      
   // Needs to be k*b^n-1 or k*b^n+1
   if (c != -1 && c != 1)
      return false;
   
   // If c = -1, then b must be 2^g for some g
   if (c == -1 && broot != 2)
      return false;

   if (k > 1)
   {
      GetRoot(k, &kroot, &kpower);
   
      // If c = -1, then k must be 2^f for some f
      if (c == -1 && kroot != 2)
         return false;
      
      // If c = +1, then k and b must have the same root
      if (c == +1 && broot != kroot)
         return false;
   }

   return true;
}

