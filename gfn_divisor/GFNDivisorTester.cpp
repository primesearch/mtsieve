/* GFNDivisorTester.cpp -- (C) Mark Rodenkirch, February 2021

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <time.h>
#include "GFNDivisorTester.h"
#include "GFNDivisorApp.h"
#include "../core/Clock.h"
#include "../core/inline.h"
#include "../x86_asm_ext/asm-ext-x86.h"

#define BIT(k)          (((k) - il_MinK) >> 1)

// Set PRE_SQUARE=N to compute 2^2^n as (2^2^N)^2^(n-N), which saves N
// sqrmods at a cost of more time in mpn_tdiv_qr().  N must satisfy 0 <= N <= 5.
#define PRE_SQUARE 5

// Handle the possibility that mp_limb_t and uint64_t are different typedefs.
// If they are they the compiler will not be able to cast.
#ifdef WIN32
static const unsigned long long _MONTGOMERY_DATA[14] = {1,0,0,0,0,0,0,0,0,0,0,0,0,UINT64_C(1)<<(1<<PRE_SQUARE)};
#else
static const unsigned long _MONTGOMERY_DATA[14] = {1,0,0,0,0,0,0,0,0,0,0,0,0,UINT64_C(1)<<(1<<PRE_SQUARE)};
#endif

#define ONE    (_MONTGOMERY_DATA+0)
#define TWO128 (_MONTGOMERY_DATA+11)
#define TWO192 (_MONTGOMERY_DATA+10)
#define TWO256 (_MONTGOMERY_DATA+9)
#define TWO320 (_MONTGOMERY_DATA+8)
#define TWO384 (_MONTGOMERY_DATA+7)
#define TWO448 (_MONTGOMERY_DATA+6)
#define TWO512 (_MONTGOMERY_DATA+5)
#define TWO576 (_MONTGOMERY_DATA+4)
#define TWO640 (_MONTGOMERY_DATA+3)
#define TWO704 (_MONTGOMERY_DATA+2)
#define TWO768 (_MONTGOMERY_DATA+1)

GFNDivisorTester::GFNDivisorTester(App *theApp)
{
   ip_App = theApp;
   
   GFNDivisorApp *gfnApp = (GFNDivisorApp *) theApp;
   
   il_MinK = gfnApp->GetMinK();
   il_MaxK = gfnApp->GetMaxK();
   
   ii_MinN = gfnApp->GetMinN();
   ii_MaxN = gfnApp->GetMaxN();
   
   iv_Terms = gfnApp->GetTerms();
}

void  GFNDivisorTester::StartedSieving(void)
{
   il_StartSievingUS = Clock::GetCurrentMicrosecond();
}

void  GFNDivisorTester::TestRemainingTerms(uint64_t totalTerms, uint64_t termsInChunk, uint64_t termCount)
{
   time_t   lastCheckPointTime = time(NULL);
   uint64_t sievingUS, startTestingUS, currentUS;
   uint64_t calculationUS;
   uint64_t termsTested = 0;
   uint64_t termsEvaluatedInChunk = 0;
   uint64_t termsTestedPerSecond;
   double   percentSievingTimeSlice;
   double   percentTermsRequiringTest;
   double   percentChunkTested, percentRangeTested;
   FILE    *fPtr;
   mpz_t    rem, fermat, nTemp, kTemp, factor, minus1;

   startTestingUS = Clock::GetCurrentMicrosecond();
   sievingUS = startTestingUS - il_StartSievingUS;
            
   mpz_init(rem);
   mpz_init(fermat);
   mpz_init(nTemp);
   mpz_init(kTemp);
   mpz_init(factor);
   mpz_init(minus1);

   mpz_set_ui(fermat, 2);
   
   for (uint32_t n=ii_MinN; n<=ii_MaxN; n++)
   {
      uint64_t bit = BIT(il_MinK);
         
      mpz_set_ui(nTemp, 2);
      mpz_pow_ui(nTemp, nTemp, n);
         
      for (uint64_t k=il_MinK; k<=il_MaxK; k+=2, bit++)
      {
         il_TotalTermsEvaluated++;
         termsEvaluatedInChunk++;
         
         if (!iv_Terms[n-ii_MinN][bit])
            continue;

         termsTested++;

         if (time(NULL) > lastCheckPointTime + 60)
         {            
            percentSievingTimeSlice = ((double) termsEvaluatedInChunk) / (double) termsInChunk;
            percentChunkTested = 100.0 * percentSievingTimeSlice;
            
            currentUS = Clock::GetCurrentMicrosecond();
            
            // Add the time to test this range to the time to sieve this range.  So if it took 180 seconds
            // to sieve this chunk and we have tested 20 percent of this chunk then "assign" 36 seconds
            // (as 36 is 20% of 180) of sieving time to this chunk.
            calculationUS = (sievingUS * percentSievingTimeSlice) + (currentUS - startTestingUS);
                        
            percentTermsRequiringTest = (100.0 * (double) termCount) / (double) termsInChunk;
            
            termsTestedPerSecond = termsEvaluatedInChunk / (calculationUS / 1000000);
            
            // We really didn't evaluate even k, but we count against the rate anyways.
            termsTestedPerSecond *= 2;

            if (termsInChunk == totalTerms)
               ip_App->WriteToConsole(COT_SIEVE, "Tested %5.2f pct of range at %" PRIu64" terms per second (%5.2f pct terms passed sieving)",
                                      percentChunkTested, termsTestedPerSecond, percentTermsRequiringTest);
            else
            {
               percentRangeTested = (100.0 * (double) il_TotalTermsEvaluated) / (double) totalTerms;
            
               ip_App->WriteToConsole(COT_SIEVE, "Tested %5.2f pct of chunk at %" PRIu64" terms per second (%5.2f pct terms passed sieving) (%5.2f pct of range)",
                                      percentChunkTested, termsTestedPerSecond, percentTermsRequiringTest, percentRangeTested);
            }
      
            lastCheckPointTime = time(NULL);
         }
         
         if (IsFermatDivisor(k, n))
         {
#ifdef WIN32
            // Even though build with 64-bit limbs, mpz_set_ui doesn't
            // populate kTemp correctly when k > 32 bits.
            mpz_set_ui(kTemp, k >> 32);
            mpz_mul_2exp(kTemp, kTemp, 32);
            mpz_add_ui(kTemp, kTemp, k & (0xffffffff));
#else
            mpz_set_ui(kTemp, k);
#endif

            mpz_mul_2exp(minus1, kTemp, n);
            mpz_add_ui(factor, minus1, 1);

            mpz_powm(rem, fermat, nTemp, factor);;
         
            if (mpz_cmp(rem, minus1) == 0)
            {
               ip_App->WriteToConsole(COT_OTHER, "Found factor %" PRIu64"*2^%u+1 of 2^(2^%u)+1", k, n, n-2);
               
               fPtr = fopen("gfn_factors.txt", "a+");
               fprintf(fPtr, "Found factor %" PRIu64"*2^%u+1 of 2^(2^%u)+1\n", k, n, n-2);
               fclose(fPtr);
            }
            else if (mpz_cmp_ui(rem, 1) == 0)
            {
               mpz_set_ui(rem, 2);
               
               for (uint32_t m=1; m<=n-2; m++)
               {
                  mpz_powm_ui(rem, rem, 2, factor);
                  
                  if (mpz_cmp(rem, minus1) == 0)
                  {
                     ip_App->WriteToConsole(COT_OTHER, "Found factor %" PRIu64"*2^%u+1 of 2^(2^%u)+1", k, n, m);
                     
                     fPtr = fopen("gfn_factors.txt", "a+");
                     fprintf(fPtr, "Found factor %" PRIu64"*2^%u+1 of 2^(2^%u)+1\n", k, n, m);
                     fclose(fPtr);
                  }
               }
            }
         }
      }
   }
   
   mpz_clear(rem);
   mpz_clear(fermat);
   mpz_clear(nTemp);
   mpz_clear(kTemp);
   mpz_clear(factor);
   mpz_clear(minus1);

   currentUS = Clock::GetCurrentMicrosecond();
   
   // If we took less than 5 seconds to test the range, then the terms per second
   // calculation is rather meaningless so we won't show it.
   if (currentUS - il_StartSievingUS > 50000000)
   {
      percentChunkTested = (100.0 * (double) termsEvaluatedInChunk) / (double) termsInChunk;
      percentTermsRequiringTest = (100.0 * (double) termCount) / (double) termsInChunk;
      
      termsTestedPerSecond = termsEvaluatedInChunk / ((currentUS - il_StartSievingUS) / 1000000);

      // We really didn't evaluate even k, but we count against the rate anyways.
      termsTestedPerSecond *= 2;
               
      if (termsInChunk == totalTerms)
         ip_App->WriteToConsole(COT_SIEVE, "Tested %5.2f pct of range at %" PRIu64" terms per second (%5.2f pct terms passed sieving)",
                                percentChunkTested, termsTestedPerSecond, percentTermsRequiringTest);
      else
      {
         percentRangeTested = (100.0 * (double) il_TotalTermsEvaluated) / (double) totalTerms;
      
         ip_App->WriteToConsole(COT_SIEVE, "Tested %5.2f pct of chunk at %" PRIu64" terms per second (%5.2f pct terms passed sieving) (%5.2f pct of range)",
                                percentChunkTested, termsTestedPerSecond, percentTermsRequiringTest, percentRangeTested);
      }
   }
}

// Return 1 iff k*2^n+1 is a Fermat divisor.
// This code is from fermat_redc.c of GMP-Fermat
bool  GFNDivisorTester::IsFermatDivisor(uint64_t k, uint32_t n)
{
// Handle the possibility that mp_limb_t and uint64_t are different typedefs.
// We don't use mp_limb_t here, but the variables are passed to our assembler
// ext functions and to GMP functions, so theu must be compatible.
#ifdef WIN32
   unsigned long long inv, A[14], N[12];
#else
   unsigned long inv, A[14], N[12];
#endif

   uint32_t i;


   // If 2^2^(n-1) = 1 (mod k*2^n+1) then k*2^n+1 divides 2^2^(n-1)-1 and
   // hence divides at least one of 2^2^1+1, 2^2^2+1, ..., 2^2^(n-2)+1.
   switch (n/64)
   {
      case 0: /* n < 64, so k*2^n+1 < 2^127 */
         /* N <-- k*2^n+1 */
         N[0] = (k << n) | 1;
         N[1] = (k >> (64-n));

         /* Encode 2 into Montgomery form, A <-- 2*2^128 mod N */
         A[1] = 0;
         
         mpn_tdiv_qr(A+2,A,0L,TWO128,3,N,(N[1]>0)?2:1);

         /* A <-- 2^2^(n-1) mod N, in Montgomery form */
         inv = invproth(k,n);
         for (i = n-1-PRE_SQUARE; i > 0; i--)
            sqrmod128_proth0(A,N,inv);

         /* Decode result from Montgomery form */
         mulmod128(A,ONE,A,N,inv);
         CheckRedc(A,2,2,n-1,k,n);

         /* Return 1 iff k*2^n+1 is a Fermat divisor */
         return (A[0] == 1 && A[1] == 0);

      case 1: /* n < 128 */
         if (n == 64 || (k >> (128-n)) == 0) /* 2^64 < k*2^n+1 < 2^128 */
         {
            N[0] = 1;
            N[1] = (k << (n-64));
            
            mpn_tdiv_qr(A+2,A,0L,TWO128,3,N,2);
            
            for (i = n-1-PRE_SQUARE; i > 0; i--)
               sqrmod128_proth1(A,N);
         
            mulmod128_proth1(A,ONE,A,N);
            CheckRedc(A,2,2,n-1,k,n);
            
            return (A[0] == 1 && A[1] == 0);
         }
         else /* 2^128 < k*2^n+1 < 2^191 */
         {
            N[0] = 1;
            N[1] = (k << (n-64));
            N[2] = (k >> (128-n));
            
            mpn_tdiv_qr(A+3,A,0L,TWO192,4,N,3);
            
            for (i = n-1-PRE_SQUARE; i > 0; i--)
               sqrmod192_proth0(A,N);
         
            mulmod192_proth0(A,ONE,A,N);
            
            CheckRedc(A,3,2,n-1,k,n);
            
            return (A[0] == 1 && A[1] == 0 && A[2] == 0);
         }

      case 2: /* n < 192 */
         if (n == 128 || (k >> (192-n)) == 0) /* 2^128 < k*2^n+1 < 2^192 */
         {
            N[0] = 1; N[1] = 0;
            N[2] = (k << (n-128));
            
            mpn_tdiv_qr(A+3,A,0L,TWO192,4,N,3);
            
            for (i = n-1-PRE_SQUARE; i > 0; i--)
               sqrmod192_proth1(A,N);
         
            mulmod192_proth1(A,ONE,A,N);
            
            CheckRedc(A,3,2,n-1,k,n);
            
            return (A[0] == 1 && A[1] == 0 && A[2] == 0);
         }
         else /* 2^192 < k*2^n+1 < 2^255 */
         {
            N[0] = 1; N[1] = 0;
            N[2] = (k << (n-128));
            N[3] = (k >> (192-n));
            
            mpn_tdiv_qr(A+4,A,0L,TWO256,5,N,4);
            
            for (i = n-1-PRE_SQUARE; i > 0; i--)
               sqrmod256_proth0(A,N);
            
            mulmod256_proth0(A,ONE,A,N);
            
            CheckRedc(A,4,2,n-1,k,n);
            
            return (A[0] == 1 && A[1] == 0 && A[2] == 0 && A[3] == 0);
         }

      case 3: /* n < 256 */
         if (n == 192 || (k >> (256-n)) == 0) /* 2^192 < k*2^n+1 < 2^256 */
         {
            N[0] = 1; N[1] = 0; N[2] = 0;
            N[3] = (k << (n-192));
            
            mpn_tdiv_qr(A+4,A,0L,TWO256,5,N,4);
            
            for (i = n-1-PRE_SQUARE; i > 0; i--)
               sqrmod256_proth1(A,N);
            
            mulmod256_proth1(A,ONE,A,N);
            
            CheckRedc(A,4,2,n-1,k,n);
            
            return (A[0] == 1 && A[1] == 0 && A[2] == 0 && A[3] == 0);
         }
         else /* 2^256 < k*2^n+1 < 2^319 */
         {
            N[0] = 1; N[1] = 0; N[2] = 0;
            N[3] = (k << (n-192));
            N[4] = (k >> (256-n));
            
            mpn_tdiv_qr(A+5,A,0L,TWO320,6,N,5);
            
            for (i = n-1-PRE_SQUARE; i > 0; i--)
               mulmod320_proth0(A,A,A,N);
         
            mulmod320_proth0(A,ONE,A,N);
            
            CheckRedc(A,5,2,n-1,k,n);
            
            return (A[0] == 1 && A[1] == 0 && A[2] == 0 && A[3] == 0 && A[4] == 0);
         }

      case 4: /* n < 320 */
         if (n == 256 || (k >> (320-n)) == 0) /* 2^256 < k*2^n+1 < 2^320 */
         {
            N[0] = 1; N[1] = 0; N[2] = 0; N[3] = 0;
            N[4] = (k << (n-256));
            
            mpn_tdiv_qr(A+5,A,0L,TWO320,6,N,5);
            
            for (i = n-1-PRE_SQUARE; i > 0; i--)
               mulmod320_proth1(A,A,A,N);
         
            mulmod320_proth1(A,ONE,A,N);
            
            CheckRedc(A,5,2,n-1,k,n);
            
            return (A[0] == 1 && A[1] == 0 && A[2] == 0 && A[3] == 0 && A[4] == 0);
         }
         else /* 2^320 < k*2^n+1 < 2^383 */
         {
            N[0] = 1; N[1] = 0; N[2] = 0; N[3] = 0;
            N[4] = (k << (n-256));
            N[5] = (k >> (320-n));
            
            mpn_tdiv_qr(A+6,A,0L,TWO384,7,N,6);
            
            for (i = n-1-PRE_SQUARE; i > 0; i--)
               mulmod384_proth0(A,A,A,N);
            
            mulmod384_proth0(A,ONE,A,N);
            
            CheckRedc(A,6,2,n-1,k,n);
            
            return (A[0] == 1 && A[1] == 0 && A[2] == 0 && A[3] == 0 && A[4] == 0 && A[5] == 0);
         }

      case 5: /* n < 384 */
         if (n == 320 || (k >> (384-n)) == 0) /* 2^320 < k*2^n+1 < 2^384 */
         {
            N[0] = 1; N[1] = 0; N[2] = 0; N[3] = 0; N[4] = 0;
            N[5] = (k << (n-320));
            
            mpn_tdiv_qr(A+6,A,0L,TWO384,7,N,6);
            
            for (i = n-1-PRE_SQUARE; i > 0; i--)
               mulmod384_proth1(A,A,A,N);
         
            mulmod384_proth1(A,ONE,A,N);
            
            CheckRedc(A,6,2,n-1,k,n);
            
            return (A[0] == 1 && A[1] == 0 && A[2] == 0 && A[3] == 0 && A[4] == 0 && A[5] == 0);
         }
         else /* 2^384 < k*2^n+1 < 2^447 */
         {
            N[0] = 1; N[1] = 0; N[2] = 0; N[3] = 0; N[4] = 0;
            N[5] = (k << (n-320));
            N[6] = (k >> (384-n));
            
            mpn_tdiv_qr(A+7,A,0L,TWO448,8,N,7);
            
            for (i = n-1-PRE_SQUARE; i > 0; i--)
               mulmod448_proth0(A,A,A,N);
         
            mulmod448_proth0(A,ONE,A,N);
            
            CheckRedc(A,7,2,n-1,k,n);
            
            return (A[0] == 1 && A[1] == 0 && A[2] == 0 && A[3] == 0 && A[4] == 0 && A[5] == 0 && A[6] == 0);
         }

      case 6: /* n < 448 */
         if (n == 384 || (k >> (448-n)) == 0) /* 2^384 < k*2^n+1 < 2^448 */
         {
            N[0] = 1; N[1] = 0; N[2] = 0; N[3] = 0; N[4] = 0; N[5] = 0;
            N[6] = (k << (n-384));
            
            mpn_tdiv_qr(A+7,A,0L,TWO448,8,N,7);
            
            for (i = n-1-PRE_SQUARE; i > 0; i--)
               mulmod448_proth1(A,A,A,N);
         
            mulmod448_proth1(A,ONE,A,N);
            
            CheckRedc(A,7,2,n-1,k,n);
            
            return (A[0] == 1 && A[1] == 0 && A[2] == 0 && A[3] == 0 && A[4] == 0 && A[5] == 0 && A[6] == 0);
         }
         else /* 2^448 < k*2^n+1 < 2^511 */
         {
            N[0] = 1; N[1] = 0; N[2] = 0; N[3] = 0; N[4] = 0; N[5] = 0;
            N[6] = (k << (n-384));
            N[7] = (k >> (448-n));
            
            mpn_tdiv_qr(A+8,A,0L,TWO512,9,N,8);
            
            for (i = n-1-PRE_SQUARE; i > 0; i--)
               mulmod512_proth0(A,A,A,N);
         
            mulmod512_proth0(A,ONE,A,N);
            
            CheckRedc(A,8,2,n-1,k,n);
            
            return (A[0] == 1 && A[1] == 0 && A[2] == 0 && A[3] == 0 && A[4] == 0 &&
                    A[5] == 0 && A[6] == 0 && A[7] == 0);
         }

      case 7: /* n < 512 */
         if (n == 448 || (k >> (512-n)) == 0) /* 2^448 < k*2^n+1 < 2^512 */
         {
            N[0] = 1; N[1] = 0; N[2] = 0; N[3] = 0; N[4] = 0; N[5] = 0;
            N[6] = 0;
            N[7] = (k << (n-448));
            
            mpn_tdiv_qr(A+8,A,0L,TWO512,9,N,8);
            
            for (i = n-1-PRE_SQUARE; i > 0; i--)
               mulmod512_proth1(A,A,A,N);
         
            mulmod512_proth1(A,ONE,A,N);
            
            CheckRedc(A,8,2,n-1,k,n);
            
            return (A[0] == 1 && A[1] == 0 && A[2] == 0 && A[3] == 0 && A[4] == 0 &&
                    A[5] == 0 && A[6] == 0 && A[7] == 0);
         }
         else /* 2^512 < k*2^n+1 < 2^575 */
         {
            N[0] = 1; N[1] = 0; N[2] = 0; N[3] = 0; N[4] = 0; N[5] = 0;
            N[6] = 0;
            N[7] = (k << (n-448));
            N[8] = (k >> (512-n));
            
            mpn_tdiv_qr(A+9,A,0L,TWO576,10,N,9);
            
            for (i = n-1-PRE_SQUARE; i > 0; i--)
               mulmod576_proth0(A,A,A,N);
         
            mulmod576_proth0(A,ONE,A,N);
            
            CheckRedc(A,9,2,n-1,k,n);
            
            return (A[0] == 1 && A[1] == 0 && A[2] == 0 && A[3] == 0 && A[4] == 0 &&
                    A[5] == 0 && A[6] == 0 && A[7] == 0 && A[8] == 0);
         }

      case 8: /* n < 576 */
         if (n == 512 || (k >> (576-n)) == 0) /* 2^512 < k*2^n+1 < 2^576 */
         {
            N[0] = 1; N[1] = 0; N[2] = 0; N[3] = 0; N[4] = 0; N[5] = 0;
            N[6] = 0; N[7] = 0;
            N[8] = (k << (n-512));
            
            mpn_tdiv_qr(A+9,A,0L,TWO576,10,N,9);
            
            for (i = n-1-PRE_SQUARE; i > 0; i--)
               mulmod576_proth1(A,A,A,N);
         
            mulmod576_proth1(A,ONE,A,N);
            
            CheckRedc(A,9,2,n-1,k,n);
            
            return (A[0] == 1 && A[1] == 0 && A[2] == 0 && A[3] == 0 && A[4] == 0 &&
                    A[5] == 0 && A[6] == 0 && A[7] == 0 && A[8] == 0);
         }
         else /* 2^576 < k*2^n+1 < 2^639 */
         {
            N[0] = 1; N[1] = 0; N[2] = 0; N[3] = 0; N[4] = 0; N[5] = 0;
            N[6] = 0; N[7] = 0;
            N[8] = (k << (n-512));
            N[9] = (k >> (576-n));
            
            mpn_tdiv_qr(A+10,A,0L,TWO640,11,N,10);
            
            for (i = n-1-PRE_SQUARE; i > 0; i--)
               mulmod640_proth0(A,A,A,N);
         
            mulmod640_proth0(A,ONE,A,N);
            
            CheckRedc(A,10,2,n-1,k,n);
            
            return (A[0] == 1 && A[1] == 0 && A[2] == 0 && A[3] == 0 && A[4] == 0 &&
                    A[5] == 0 && A[6] == 0 && A[7] == 0 && A[8] == 0 && A[9] == 0);
         }

      case 9: /* n < 640 */

         if (n == 576 || (k >> (640-n)) == 0) /* 2^576 < k*2^n+1 < 2^640 */
         {
            N[0] = 1; N[1] = 0; N[2] = 0; N[3] = 0; N[4] = 0; N[5] = 0;
            N[6] = 0; N[7] = 0; N[8] = 0;
            N[9] = (k << (n-576));

            mpn_tdiv_qr(A+10,A,0L,TWO640,11,N,10);

            for (i = n-1-PRE_SQUARE; i > 0; i--)
               mulmod640_proth1(A,A,A,N);

            mulmod640_proth1(A,ONE,A,N);
            
            CheckRedc(A,10,2,n-1,k,n);
            
            return (A[0] == 1 && A[1] == 0 && A[2] == 0 && A[3] == 0 && A[4] == 0 && 
                    A[5] == 0 && A[6] == 0 && A[7] == 0 && A[8] == 0 && A[9] == 0);
         }
         else /* 2^640 < k*2^n+1 < 2^703 */
         {
            N[0] = 1; N[1] = 0; N[2] = 0; N[3] = 0; N[4] = 0; N[5] = 0;
            N[6] = 0; N[7] = 0; N[8] = 0;
            N[9] = (k << (n-576));
            N[10] = (k >> (640-n));
            
            mpn_tdiv_qr(A+11,A,0L,TWO704,12,N,11);
            
            for (i = n-1-PRE_SQUARE; i > 0; i--)
               mulmod704_proth0(A,A,A,N);
         
            mulmod704_proth0(A,ONE,A,N);
            
            CheckRedc(A,11,2,n-1,k,n);
            
            return (A[0] == 1 && A[1] == 0 && A[2] == 0 && A[3] == 0 && A[4] == 0 &&
                    A[5] == 0 && A[6] == 0 && A[7] == 0 && A[8] == 0 && A[9] == 0 && A[10] == 0);
         }

      case 10: /* n < 704 */
         if (n == 640 || (k >> (704-n)) == 0) /* 2^640 < k*2^n+1 < 2^704 */
         {
            N[0] = 1; N[1] = 0; N[2] = 0; N[3] = 0; N[4] = 0; N[5] = 0;
            N[6] = 0; N[7] = 0; N[8] = 0; N[9] = 0;
            N[10] = (k << (n-640));
            
            mpn_tdiv_qr(A+11,A,0L,TWO704,12,N,11);

            for (i = n-1-PRE_SQUARE; i > 0; i--)
               mulmod704_proth1(A,A,A,N);

            mulmod704_proth1(A,ONE,A,N);

            CheckRedc(A,11,2,n-1,k,n);

            return (A[0] == 1 && A[1] == 0 && A[2] == 0 && A[3] == 0 && A[4] == 0 &&
                    A[5] == 0 && A[6] == 0 && A[7] == 0 && A[8] == 0 &&
                    A[9] == 0 && A[10] == 0);
         }
         else /* 2^704 < k*2^n+1 < 2^767 */
         {
            N[0] = 1; N[1] = 0; N[2] = 0; N[3] = 0; N[4] = 0; N[5] = 0;
            N[6] = 0; N[7] = 0; N[8] = 0; N[9] = 0;
            N[10] = (k << (n-640));
            N[11] = (k >> (704-n));

            mpn_tdiv_qr(A+12,A,0L,TWO768,13,N,12);

            for (i = n-1-PRE_SQUARE; i > 0; i--)
               mulmod768_proth0(A,A,A,N);

            mulmod768_proth0(A,ONE,A,N);

            CheckRedc(A,12,2,n-1,k,n);

            return (A[0] == 1 && A[1] == 0 && A[2] == 0 && A[3] == 0 && A[4] == 0 &&
                    A[5] == 0 && A[6] == 0 && A[7] == 0 && A[8] == 0 &&
                    A[9] == 0 && A[10] == 0 && A[11] == 0);
         }

      case 11: /* n < 768 */
         if (n == 704 || (k >> (768-n)) == 0) /* 2^704 < k*2^n+1 < 2^768 */
         {
            N[0] = 1; N[1] = 0; N[2] = 0; N[3] = 0; N[4] = 0; N[5] = 0; N[6] = 0;
            N[7] = 0; N[8] = 0; N[9] = 0; N[10] = 0;
            N[11] = (k << (n-704));
            
            mpn_tdiv_qr(A+12,A,0L,TWO768,13,N,12);
            
            for (i = n-1-PRE_SQUARE; i > 0; i--)
               mulmod768_proth1(A,A,A,N);
         
            mulmod768_proth1(A,ONE,A,N);
            
            CheckRedc(A,12,2,n-1,k,n);
            
            return (A[0] == 1 && A[1] == 0 && A[2] == 0 && A[3] == 0 && A[4] == 0 &&
                    A[5] == 0 && A[6] == 0 && A[7] == 0 && A[8] == 0 &&
                    A[9] == 0 && A[10] == 0 && A[11] == 0);
         }

      default:
      {
         uint32_t nn = n/64;
         mp_limb_t np[nn+2], t1[2*nn+4], t2[2*nn+4];

         np[0] = 1;
         
         for (i = 1; i < nn; i++)
            np[i] = 0;
         
         np[nn++] = k<<(n%64);
         
         if (n%64)
            if ((np[nn] = k>>(64-n%64)) != 0)
               nn++;

         for (i = 0; i < nn; i++)
            t1[i] = 0;
         
         t1[nn] = (mp_limb_t)1<<(1<<PRE_SQUARE);
         mpn_tdiv_qr(t2+nn,t2,0L,t1,nn+1,np,nn);

         for (i = n-1-PRE_SQUARE; i > 0; i--)
         {
            mpn_mul_n(t1,t2,t2,nn);
            redc_proth0(t1,np,nn);
            if (mpn_add_n(t2,t1+nn,t1,nn))
               mpn_sub_n(t2,t2,np,nn);
         }
         
         for (i = nn; i < 2*nn; i++)
            t2[i] = 0;
         
         redc_proth0(t2,np,nn);
         
         if (mpn_add_n(t1,t2+nn,t2,nn))
            mpn_sub_n(t1,t1,np,nn);
      
#if 1
         // Is this needed?
         while(mpn_cmp(t1,np,nn) >= 0)
            mpn_sub_n(t1,t1,np,nn);
#endif
         
         CheckRedc(t1,nn,2,n-1,k,n);
         
         if (t1[0] != 1)
            return false;
      
         for (i = 1; i < nn; i++)
            if (t1[i] != 0)
               return false;
         
         return true;
      }
   }
}

void  GFNDivisorTester::CheckRedc(mp_limb_t *xp, uint32_t xn, uint32_t b, uint32_t m, uint64_t k, uint32_t n)
{
#ifndef TEST_REDC
   return;
#endif

   mpz_t     B, E, N;
   uint32_t  i;
   
   mpz_init_set_ui(B,b);

   mpz_init_set_ui(E,1);
   mpz_mul_2exp(E,E,m);

   mpz_init(N);
   
#ifdef WIN32
   // Even though build with 64-bit limbs, mpz_set_ui doesn't
   // populate kTemp correctly when k > 32 bits.
   mpz_set_ui(N, k >> 32);
   mpz_mul_2exp(N, N, 32);
   mpz_add_ui(N, N, k & (0xffffffff));
#else
   mpz_set_ui(N, k);
#endif

   mpz_mul_2exp(N,N,n);
   mpz_add_ui(N,N,1);

   mpz_powm(B,B,E,N);

   if (mpz_size(B) > xn)
   {
      ip_App->WriteToConsole(COT_OTHER, "REDC ERROR: %u^2^%u mod %" PRIu64"*2^%u+1\n", b, m, k, n);
      
      FatalError("mpz_size=%" PRIu64", xn=%u\n", (uint64_t) mpz_size(B), xn);
   }

   for (i = 0; i < xn; i++)
      if (xp[i] != mpz_getlimbn(B,i))
      {         
         ip_App->WriteToConsole(COT_OTHER, "REDC ERROR: %u^2^%u mod %" PRIu64"*2^%u+1\n", b, m, k, n);
         
         FatalError("limb %u=%" PRIu64", xp[%u]=%" PRIu64"\n", i, (uint64_t) mpz_getlimbn(B, i), i, (uint64_t) xp[i]);
      }

   mpz_clear(N);
   mpz_clear(E);
   mpz_clear(B);
}

