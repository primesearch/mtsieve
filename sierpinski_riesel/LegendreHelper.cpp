/* LegendreHelper.cpp -- (C) Mark Rodenkirch, July 2019

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <stdint.h>

#include "../core/inline.h"
#include "../sieve/primesieve.hpp"
#include "../x86_asm/fpu-asm-x86.h"
#include "../x86_asm/sse-asm-x86.h"
#include "LegendreHelper.h"

LegendreHelper::LegendreHelper(App *theApp, uint32_t base, seq_t *sequences, uint32_t sequenceCount)
{
   uint64_t max_k = 0;
   uint32_t seqIdx;
   
   ip_App = theApp;
   
   ii_Base = base;
   ip_Sequences = sequences;
   ii_SequenceCount = sequenceCount;
   
   iv_Primes.clear();
   
   for (seqIdx=0; seqIdx<ii_SequenceCount; seqIdx++)
      max_k = MAX(max_k, ip_Sequences[seqIdx].k);
   
   // We're going to use this later.
   primesieve::generate_n_primes(sqrt(max_k) + 1, &iv_Primes);
}

// Set fields for sequence k*b^n+c so that bit (p/2)%mod is set if and only if:
//   (-ck/p)=1 for sequences with all n even,
//   (-bck/p)=1 for sequences with all n odd,
//   (-ck/p)=1 or (-bck/p)=1 for sequences with both odd and even n.
//
// In the worst case the table for k*b^n+c could be 4*b*k bits long.
bool  LegendreHelper::BuildLegendreTables(void)
{
   uint32_t b, i, y;
   uint64_t k, m;
   uint32_t seqIdx;
   int64_t r;
   int64_t c;

   b = (uint32_t) SquareFreeFactorization(ii_Base);
      
   // For each sequence build a map to hold the results of the jacobi check.   
   // We need to be careful about the size of the map for two reasons:
   //    1) It would require a large amount of memory for the sequence
   //    2) It would take a long time to build the map
   //
   // For now, if the map is too large, we won't build it for any sequences.  This won't
   // prevent the usage of this program.  It only means that it won't perform as well.
   //
   // In the future the building of the map could be off-loaded to a GPU to reduce the time,
   // but there would still be concerns about how much memory is required to hold the map.
   for (seqIdx=0; seqIdx<ii_SequenceCount; seqIdx++)
   {
      k = ip_Sequences[seqIdx].k;
      c = ip_Sequences[seqIdx].c;
      m = SquareFreeFactorization(k);

      // We may need the signed product c*b*m, and the unsigned product m*b*2.
      r = (c < 0) ? m : -m;

      switch (ip_Sequences[seqIdx].parity)
      {
         // odd n, test for (-bck/p) == 1
         case 1:
            m *= b;
            r *= b;
            // fall through
            
         // even n, test for (-ck/p) == 1
         case 0:
            if ((r < 0 && (-r) % 4 != 3) || (r > 0 && r % 4 != 1))
               m *= 2;

            if ((2*m+1) >= INT32_MAX)
            {
               theApp->WriteToConsole(COT_OTHER, "Square-free part of sequence %" PRIu64"*%u^n%+d is too large so the Legendre logic has been disabled", k, ii_Base, c);
               return false;
            }

            theApp->WriteToConsole(COT_SIEVE, "Buiding Legendre symbol lookup table for %" PRIu64"*%u^n%+d", k, ii_Base, c);
            
            ip_Sequences[seqIdx].legendreMap.resize(m);
            
            for (i=0; i<m; i++)
               if (jacobi(r, 2*i+1) == 1)
                  ip_Sequences[seqIdx].legendreMap(i) = true;
               
            break;

         // even and odd n, test for (-ck/p)==1 or (-bck/p)==1
         case 2: 
            m = m*2*b;

            if ((2*m+1) >= INT32_MAX || r*b >= INT32_MAX)
            {
               theApp->WriteToConsole(COT_OTHER, "Square-free part of sequence %" PRIu64"*%u^n%+d is too large so the Legendre logic has been disabled", k, ii_Base, c);
               return false;
            }

            theApp->WriteToConsole(COT_SIEVE, "Buiding Legendre symbol lookup table for %" PRIu64"*%u^n%+d", k, ii_Base, c);
            
            ip_Sequences[seqIdx].legendreMap.resize(m);
            
            for (i=0; i<m; i++)
               if (jacobi(r, 2*i+1) == 1 || jacobi(r*b, 2*i+1) == 1)
                  ip_Sequences[seqIdx].legendreMap(i) = true;
               
            break;
      }
   }
}

// This will return the product of the prime factors of the square free part of n.
//
// Examples:
//    n =  18 --> 3^2 * 2       --> return 2
//    n =  27 --> 3^2 * 3       --> return 3
//    n =  28 --> 2^2 * 7       --> return 7
//    n =  36 --> 2^2 * 3^2 * 1 --> return 1
//    n =  91 --> 7 * 13        --> return 91
//    n = 180 --> 2^2 * 3^2 * 5 --> return 5
uint64_t    LegendreHelper::SquareFreeFactorization(uint64_t n)
{
   uint64_t c = 1, d, q, r;
   
   vector<uint64_t>::iterator it = iv_Primes.begin();
   
   r = sqrt(n);
   
   while (it != iv_Primes.end())
   {
      q = *it;
      it++;
      
      if (n % q != 0)
         continue;
   
      while (n % q == 0)
      {
         n /= q;
         
         if (n % q != 0)
         {
            c *= q;
            break;
         }
         
         n /= q;
      };
      
      r = sqrt(n);
      
      if (r*r == n)
         return c;
   }
   
   return c * n;
}
