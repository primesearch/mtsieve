/* DMDivisorApp.cpp -- (C) Mark Rodenkirch, September 2018

   Sieve for 2*k*M+1 where M = is a Mersenne Prime.  The remaining terms,
   if prime, might be divisors of 2^(2^p-1)-1.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <time.h>
#include "../core/inline.h"
#include "../core/Parser.h"
#include "../core/Clock.h"
#include "DMDivisorApp.h"
#include "DMDivisorWorker.h"
#include "../x86_asm/fpu-asm-x86.h"
#include "../x86_asm_ext/asm-ext-x86.h"

#define APP_NAME        "dmdsieve"
#define APP_VERSION     "1.3"

#define BIT(k)          ((k) - il_MinK)

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

// This is declared in App.h, but implemented here.  This means that App.h
// can remain unchanged if using the mtsieve framework for other applications.
App *get_app(void)
{
   return new DMDivisorApp();
}

DMDivisorApp::DMDivisorApp() : FactorApp()
{
   SetBanner(APP_NAME " v" APP_VERSION ", a program to find terms that are potential factors of 2^(2^p-1)-1");
   SetLogFileName("dmdsieve.log");
   
   SetAppMinPrime(3);
   il_MinKOriginal = il_MinK = 0;
   il_MaxKOriginal = il_MaxK = 0;
   ii_N    = 0;
   
   ib_TestTerms = false;
   il_KPerChunk = 10000000000L;
   il_TotalTerms = 0;
   il_TotalTermsEvaluated = 0;
   
   iv_MMPTerms.clear();
}

void DMDivisorApp::Help(void)
{
   FactorApp::ParentHelp();

   printf("-k --kmin=k           Minimum k to search\n");
   printf("-K --kmax=K           Maximum k to search\n");
   printf("-n --exp=n            Exponent to search\n");
   printf("-x --testterms        test remaining terms for DM divisibility\n");
   printf("-X --kperchunk=X      when using -x, number of k to sieve at a time (default 1e10)\n");
}

void  DMDivisorApp::AddCommandLineOptions(std::string &shortOpts, struct option *longOpts)
{
   FactorApp::ParentAddCommandLineOptions(shortOpts, longOpts);

   shortOpts += "k:K:b:n:f:xX:";

   AppendLongOpt(longOpts, "kmin",              required_argument, 0, 'k');
   AppendLongOpt(longOpts, "kmax",              required_argument, 0, 'K');
   AppendLongOpt(longOpts, "exp",               required_argument, 0, 'n');
   AppendLongOpt(longOpts, "testterms",         no_argument,       0, 'x');
   AppendLongOpt(longOpts, "termsperchunk",     required_argument, 0, 'X');
}

parse_t DMDivisorApp::ParseOption(int opt, char *arg, const char *source)
{
   parse_t status = P_UNSUPPORTED;

   status = FactorApp::ParentParseOption(opt, arg, source);
   if (status != P_UNSUPPORTED) return status;

   switch (opt)
   {
      case 'k':
         status = Parser::Parse(arg, 1, KMAX_MAX, il_MinK);
         break;

      case 'K':
         status = Parser::Parse(arg, 1, KMAX_MAX, il_MaxK);
         break;

      case 'n':
         status = Parser::Parse(arg, 1, NMAX_MAX, ii_N);
         break;
         
      case 'x':
         ib_TestTerms = true;
         status = P_SUCCESS;
         break;
         
      case 'X':
         status = Parser::Parse(arg, 1000000000, KMAX_MAX, il_KPerChunk);
         break;
   }

   return status;
}

void DMDivisorApp::ValidateOptions(void)
{
   uint32_t mList[] = { 2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 
                        607, 1279, 2203, 2281, 3217, 4253, 4423, 9689, 9941, 
                        11213, 19937, 21701, 23209, 44497, 86243, 110503, 
                        132049, 216091, 756839, 859433, 1257787, 1398269, 
                        2976221, 3021377, 6972593, 13466917, 20996011, 24036583, 
                        25964951, 30402457, 32582657, 37156667, 42643801, 43112609, 
                        57885161, 74207281, 77232917, 82589933, 0};

#ifdef WIN32
   if (sizeof(unsigned long long) != sizeof(mp_limb_t))
     FatalError("GMP limb size is not 64 bits");
#else
   if (sizeof(unsigned long) != sizeof(mp_limb_t))
     FatalError("GMP limb size is not 64 bits");
#endif

   if (is_InputTermsFileName.length() > 0)
   {
      if (ib_TestTerms)
         FatalError("cannot use -i and -x together");

      ProcessInputTermsFile(false);

      iv_MMPTerms.resize(il_MaxK - il_MinK + 1);
      std::fill(iv_MMPTerms.begin(), iv_MMPTerms.end(), false);
      
      ProcessInputTermsFile(true);
   }
   else
   {
      if (il_MinK == 0)
         FatalError("kmin must be specified");

      if (il_MaxK == 0)
         FatalError("kmax must be specified");
      
      if (il_MaxK <= il_MinK)
         FatalError("kmax must be greater than kmin");
            
      if (ii_N == 0)
         FatalError("exponent must be specified"); 
   }
            
   if (is_OutputTermsFileName.length() == 0)
   {
      char fileName[30];
      
      sprintf(fileName, "mmp_%u.pfgw", ii_N);
      
      is_OutputTermsFileName = fileName;
   }

   int i = 0;
   while (mList[i] != 0 && mList[i] != ii_N)
      i++;
   
   if (mList[i] == 0)
     FatalError("exponent must be for a Mersenne Prime");

   if (ii_N < 13)
     FatalError("MM%u is a known prime", ii_N);
 
   if (ib_TestTerms && il_MaxPrime == il_AppMaxPrime)
      FatalError("must specify -P when testing terms");
   
   FactorApp::ParentValidateOptions();

   // Since the worker wants primes in groups of 4
   while (ii_CpuWorkSize % 4 != 0)
      ii_CpuWorkSize++;
   
   // Allow only one worker to do work when processing small primes.  This allows us to avoid 
   // locking when factors are reported, which significantly hurts performance as most terms 
   // will be removed due to small primes.
   SetMaxPrimeForSingleWorker(10000);
}

void  DMDivisorApp::PreSieveHook(void)
{
   il_TotalTerms = il_MaxK - il_MinK + 1;
   
   if (is_InputTermsFileName.length() > 0)
      return;
   
   if (!ib_TestTerms)
   {      
      il_TermCount = il_TotalTerms;
      
      iv_MMPTerms.resize(il_TermCount);
      std::fill(iv_MMPTerms.begin(), iv_MMPTerms.end(), true);   

      return;
   }
   
   if (il_MinKOriginal == 0)
   {      
      iv_MMPTerms.resize(il_KPerChunk);
            
      il_MinKOriginal = il_MinK;
      il_MaxKOriginal = il_MaxK;
      il_MinKInChunk = il_MinK;
   }
  
   std::fill(iv_MMPTerms.begin(), iv_MMPTerms.end(), true);   

   // We want il_MinK and il_MaxK to be set to the correct range
   // of k before we start sieving.
   il_MinK = il_MinKInChunk;
   il_MaxK = il_MinK + il_KPerChunk + 1;
   
   uint64_t kInChunk = il_KPerChunk;
   
   if (il_MaxK > il_MaxKOriginal)
   {
      il_MaxK = il_MaxKOriginal;
      kInChunk = il_MaxK - il_MinK + 1;
   }
   
   il_TotalTermsInChunk = il_TermCount = kInChunk;
   il_FactorCount = 0;
   
   il_StartSievingUS = Clock::GetCurrentMicrosecond();
}

bool  DMDivisorApp::PostSieveHook(void)
{
   if (!ib_TestTerms)
      return true;
   
   TestRemainingTerms();
   
   // Set the starting k for the next range to be sieved.
   il_MinKInChunk = il_MaxK;

   if (il_MinKInChunk < il_MaxKOriginal)
      return false;

   return true;
}

Worker *DMDivisorApp::CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested)
{
   Worker *theWorker;

   // Note that MMP inherits from Worker.  This will not
   // only create the worker, but also start it.
   theWorker = new DMDivisorWorker(id, this);

   return theWorker;
}

void DMDivisorApp::ProcessInputTermsFile(bool haveBitMap)
{
   FILE    *fPtr = fopen(is_InputTermsFileName.c_str(), "r");
   char     buffer[1000];
   uint32_t n, bit;
   uint64_t k, lastPrime = 0;
   
   if (!fPtr)
      FatalError("Unable to open input file %s", is_InputTermsFileName.c_str());

   if (fgets(buffer, sizeof(buffer), fPtr) == NULL)
      FatalError("No data in input file %s", is_InputTermsFileName.c_str());
   
   if (!haveBitMap)
      il_MinK = il_MaxK = 0;
   
   if (!memcmp(buffer, "ABCD ", 5))
   {
      if (sscanf(buffer, "ABCD 2*$a*(2^%u-1)+1 [%" SCNu64"] // Sieved to %" SCNu64"", &n, &k, &lastPrime) != 3)
         FatalError("Line 1 is not a valid ABCD line in input file %s", is_InputTermsFileName.c_str());
        
      ii_N = n;

      if (haveBitMap)
      {
         iv_MMPTerms[k-il_MinK] = true;
         il_TermCount++;
      }
      else
         il_MinK = il_MaxK = k;
   }
   else if (!memcmp(buffer, "ABC ", 4))
   {
      if (sscanf(buffer, "ABC 2*$a*(2^%u-1)+1 // Sieved to %" SCNu64"", &n, &lastPrime) != 2)
         FatalError("Line 1 is not a valid ABCD line in input file %s", is_InputTermsFileName.c_str());
        
      ii_N = n;
   }
   else
      FatalError("Input file %s has unknown format", is_InputTermsFileName.c_str());
   
   SetMinPrime(lastPrime);
   
   while (fgets(buffer, sizeof(buffer), fPtr) != NULL)
   {
      if (!StripCRLF(buffer))
         continue;

      if (sscanf(buffer, "%" SCNu64"", &k) != 1)
         FatalError("Line %s is malformed", buffer);
  
      if (haveBitMap)
      {
         bit = BIT(k);
      
         iv_MMPTerms[bit] = true;
         il_TermCount++;
      }
      else
      {
         if (il_MinK > k || il_MinK == 0) il_MinK = k;
         if (il_MaxK < k) il_MaxK = k;
      }
   }

   fclose(fPtr);
}

bool DMDivisorApp::ApplyFactor(uint64_t theFactor, const char *term)
{
   uint64_t k;
   uint32_t n;
      
   if (sscanf(term, "2*%" SCNu64"*(2^%u-1)+1", &k, &n) != 2)
      FatalError("Could not parse term %s", term);

   if (n != ii_N)
      FatalError("Expected n %u in factor but found %d", ii_N, n);
   
   if (k < il_MinK || k > il_MaxK)
      return false;

   VerifyFactor(theFactor, k);
      
   uint64_t bit = k - il_MinK;
   
   // No locking is needed because the Workers aren't running yet
   if (iv_MMPTerms[bit])
   {
      iv_MMPTerms[bit] = false;
      il_TermCount--;

      return true;
   }
      
   return false;
}

void DMDivisorApp::WriteOutputTermsFile(uint64_t largestPrime)
{
   uint64_t termsCounted = 0;
   uint64_t k;
   uint64_t bit;

   // With super large ranges, wait until we can lock because without locking
   // the term count can change between opening and closing the file.
   if (IsRunning() && largestPrime < GetMaxPrimeForSingleWorker())
      return;   

   k = il_MinK;
   
   bit = BIT(k);
   for (; k<=il_MaxK; k++)
   {      
      if (iv_MMPTerms[bit])
         break;
      
      bit++;
   }

   if (k > il_MaxK)
      return;
   
   FILE    *termsFile = fopen(is_OutputTermsFileName.c_str(), "w");

   if (!termsFile)
      FatalError("Unable to open output file %s", is_OutputTermsFileName.c_str());
   
   ip_FactorAppLock->Lock();

   fprintf(termsFile, "ABC 2*$a*(2^%u-1)+1 // Sieved to %" SCNu64"\n", ii_N, largestPrime);
   //fprintf(termsFile, "ABCD 2*$a*(2^%u-1)+1 [%" SCNu64"] // Sieved to %" SCNu64"\n", ii_N, k, largestPrime);
   
   for (k=il_MinK; k<=il_MaxK; k++)
   {
      bit = BIT(k);
   
      if (iv_MMPTerms[bit])
      {
         fprintf(termsFile, "%" PRIu64"\n", k);
         termsCounted++;
      }
   }
   
   fclose(termsFile);
   
   if (termsCounted != il_TermCount)
      FatalError("Something is wrong.  Counted terms (%" PRIu64") != expected terms (%" PRIu64")", termsCounted, il_TermCount);

   ip_FactorAppLock->Release();
}

void  DMDivisorApp::GetExtraTextForSieveStartedMessage(char *extraTtext)
{
   sprintf(extraTtext, "%" PRIu64 " < k < %" PRIu64", 2*k*(2^%u-1)+1", il_MinK, il_MaxK, ii_N);
}

bool  DMDivisorApp::ReportFactor(uint64_t theFactor, uint64_t k, bool verifyFactor)
{
   bool     removedTerm = false;
   char     kStr[50];

   if (theFactor > GetMaxPrimeForSingleWorker())
      ip_FactorAppLock->Lock();

   uint64_t bit = BIT(k);

   if (iv_MMPTerms[bit])
   {
      if (ii_N == 13 && theFactor == 8191)
         return false;
      
      if (ii_N == 17 && theFactor == 131071)
         return false;
      
      if (ii_N == 19 && theFactor == 524287)
         return false;
      
      if (ii_N == 31 && theFactor == 2147483647)
         return false;
      
      if (ii_N == 61 && theFactor == 2305843009213693951)
         return false;
         
      if (verifyFactor)
         VerifyFactor(theFactor, k);
      
      iv_MMPTerms[bit] = false;
      removedTerm = true;
      
      sprintf(kStr, "%" PRIu64"", k);
      
      LogFactor(theFactor, "2*%s*(2^%u-1)+1", kStr, ii_N);
      
      il_FactorCount++;
      il_TermCount--;
   }

   if (theFactor > GetMaxPrimeForSingleWorker())
      ip_FactorAppLock->Release();
   
   return removedTerm;
}

void  DMDivisorApp::TestRemainingTerms(void)
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
   mpz_t    rem, mersenne, nTemp, kTemp, factor;

   startTestingUS = Clock::GetCurrentMicrosecond();
   sievingUS = startTestingUS - il_StartSievingUS;
            
   mpz_init(rem);
   mpz_init(mersenne);
   mpz_init(nTemp);
   mpz_init(kTemp);
   mpz_init(factor);
   
   mpz_set_ui(nTemp, 2);
   mpz_pow_ui(nTemp, nTemp, ii_N);
   mpz_sub_ui(nTemp, nTemp, 1);
   mpz_set_ui(mersenne, 2);

   uint64_t bit = BIT(il_MinK);
      
   for (uint64_t k=il_MinK; k<=il_MaxK; k++, bit++)
   {
      il_TotalTermsEvaluated++;
      termsEvaluatedInChunk++;
      
      if (!iv_MMPTerms[bit])
         continue;

      termsTested++;

      if (time(NULL) > lastCheckPointTime + 60)
      {            
         percentSievingTimeSlice = ((double) termsEvaluatedInChunk) / (double) il_TotalTermsInChunk;
         percentChunkTested = 100.0 * percentSievingTimeSlice;
         
         currentUS = Clock::GetCurrentMicrosecond();
         
         // Add the time to test this range to the time to sieve this range.  So if it took 180 seconds
         // to sieve this chunk and we have tested 20 percent of this chunk then "assign" 36 seconds
         // (as 36 is 20% of 180) of sieving time to this chunk.
         calculationUS = (sievingUS * percentSievingTimeSlice) + (currentUS - startTestingUS);
                     
         percentTermsRequiringTest = (100.0 * (double) il_TermCount) / (double) il_TotalTermsInChunk;
         
         termsTestedPerSecond = termsEvaluatedInChunk / (calculationUS / 1000000);
         
         // We really didn't evaluate even k, but we count against the rate anyways.
         termsTestedPerSecond *= 2;

         if (il_TotalTermsInChunk == il_TotalTerms)
            WriteToConsole(COT_SIEVE, "Tested %5.2f pct of range at %" PRIu64" terms per second (%5.2f pct terms passed sieving)", 
                           percentChunkTested, termsTestedPerSecond, percentTermsRequiringTest);
         else
         {
            percentRangeTested = (100.0 * (double) il_TotalTermsEvaluated) / (double) il_TotalTerms;
         
            WriteToConsole(COT_SIEVE, "Tested %5.2f pct of chunk at %" PRIu64" terms per second (%5.2f pct terms passed sieving) (%5.2f pct of range)", 
                           percentChunkTested, termsTestedPerSecond, percentTermsRequiringTest, percentRangeTested);
         }
   
         lastCheckPointTime = time(NULL);
      }
      
      // Due to issues in the IsDoubleMersenneDivisor() function, we won't call it.
      // These issues exist in gmp-double-mersenne as well when the REDC code is enabled.
      
      //if (ii_N < PRE_SQUARE || ii_N >= 128 || IsDoubleMersenneDivisor(k))
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

         mpz_mul(factor, kTemp, nTemp);
         mpz_mul_ui(factor, factor, 2);
         mpz_add_ui(factor, factor, 1);

         mpz_powm(rem, mersenne, nTemp, factor);

         if (mpz_cmp_ui(rem, 1) == 0)
         {
            WriteToConsole(COT_OTHER, "Found factor 2*%" PRIu64"*(2^%u-1)+1 of 2^(2^%u-1)-1", k, ii_N, ii_N);
            
            fPtr = fopen("dm_factors.txt", "a+");
            fprintf(fPtr, "Found factor 2*%" PRIu64"*(2^%u-1)+1 of 2^(2^%u-1)-1\n", k, ii_N, ii_N);
            fclose(fPtr);
         }
      }
   }

   mpz_clear(rem);
   mpz_clear(mersenne);
   mpz_clear(nTemp);
   mpz_clear(kTemp);
   mpz_clear(factor);

   currentUS = Clock::GetCurrentMicrosecond();
   
   // If we took less than 5 seconds to test the range, then the terms per second
   // calculation is rather meaningless so we won't show it.
   if (currentUS - il_StartSievingUS > 50000000)
   {
      percentChunkTested = (100.0 * (double) termsEvaluatedInChunk) / (double) il_TotalTermsInChunk;
      percentTermsRequiringTest = (100.0 * (double) il_TermCount) / (double) il_TotalTermsInChunk;
      
      termsTestedPerSecond = termsEvaluatedInChunk / ((currentUS - il_StartSievingUS) / 1000000);

      // We really didn't evaluate even k, but we count against the rate anyways.
      termsTestedPerSecond *= 2;
               
      if (il_TotalTermsInChunk == il_TotalTerms)
         WriteToConsole(COT_SIEVE, "Tested %5.2f pct of range at %" PRIu64" terms per second (%5.2f pct terms passed sieving)",
                        percentChunkTested, termsTestedPerSecond, percentTermsRequiringTest);
      else
      {
         percentRangeTested = (100.0 * (double) il_TotalTermsEvaluated) / (double) il_TotalTerms;
      
         WriteToConsole(COT_SIEVE, "Tested %5.2f pct of chunk at %" PRIu64" terms per second (%5.2f pct terms passed sieving) (%5.2f pct of range)",
                        percentChunkTested, termsTestedPerSecond, percentTermsRequiringTest, percentRangeTested);
      }
   }
}

// Return 1 iff k*2^n+1 is a Double-Mersenne divisor.
// This code is from fermat_redc.c of GMP-Double-Mersenne.
// Note that this is limited to n < 128.
//
// This method is not called because if PRE_SQUARE > 0, then it doesn't produce
// correct results.  Also this code doesn't prevent overflows when k is large.
// I've left the code here in case someone wants to try to fix it in the future,
// but since mmff and mfac are faster, I don't see much value in that.
bool  DMDivisorApp::IsDoubleMersenneDivisor(uint64_t k)
{
// Handle the possibility that mp_limb_t and uint64_t are different typedefs.
// We don't use mp_limb_t here, but the variables are passed to our assembler
// ext functions and to GMP functions, so theu must be compatible.
#ifdef WIN32
   unsigned long long T[4], A[3], N[3], inv;
#else
   unsigned long T[4], A[3], N[3], inv;
#endif

   int i;

   k *= 2;
   
   // If 2^2^n = 2 (mod k*(2^n-1)+1) then k*(2^n-1)+1) divides 2^(2^n-1)-1.
   if (ii_N < 64) /* k*(2^n-1)+1 < 2^127 */
   {
      /* N <-- k*(2^n-1)+1 */
      N[0] = (k << ii_N) - (k-1);
      N[1] = (k >> (64-ii_N)) - ((k << ii_N) < (k-1));

      /* Encode 2 into Montgomery form, A <-- 2*2^128 mod N */
      A[1] = 0;
      mpn_tdiv_qr(T,A,0L,TWO128,3,N,(N[1]>0)?2:1);
      
#if (PRE_SQUARE==0)
      /* T <-- A */
      T[0] = A[0];
      T[1] = A[1];
#endif

      /* A <-- 2^2^n mod N, in Montgomery form */
      inv = inv_mont(N[0]);
      
      for (i = ii_N-PRE_SQUARE; i > 0; i--)
         sqrmod128_proth0(A,N,inv);

#if (PRE_SQUARE==0)
      /* Return 1 iff 2^2^n = 2 (mod N) */
      return (A[0] == T[0] && A[1] == T[1]);
#else
      /* Convert from Montgomery form */
      mulmod128(A,ONE,A,N,inv);
      CheckRedc(A,2,2,k);

      /* Return 1 iff 2^2^n = 2 (mod N) */
      return (A[0] == 2 && A[1] == 0);
#endif
  }
  
   // n >= 64 and n < 128
   N[0] = 1-k;
   N[1] = (k << (ii_N-64))-1;
   N[2] = (k >> (128-ii_N));

   if ( N[2] == 0) /* k*(2^n-1)+1 < 2^128 */
   {
      mpn_tdiv_qr(T,A,0L,TWO128,3,N,2);
#if (PRE_SQUARE==0)
      T[0] = A[0];
      T[1] = A[1];
#endif
      inv = inv_mont(N[0]);
      for (i = ii_N-PRE_SQUARE; i > 0; i--)
         sqrmod128(A,N,inv);

#if (PRE_SQUARE==0)
      return (A[0] == T[0] && A[1] == T[1]);
#else
      mulmod128(A,ONE,A,N,inv);
      CheckRedc(A,2,2,k);
      return (A[0] == 2 && A[1] == 0);
#endif
   }
   else /* k*(2^n-1)+1 < 2^191 */
   {
      mpn_tdiv_qr(T,A,0L,TWO192,4,N,3);
#if (PRE_SQUARE==0)
      T[0] = A[0];
      T[1] = A[1];
      T[2] = A[2];
#endif
      inv = inv_mont(N[0]);
      for (i = ii_N-PRE_SQUARE; i > 0; i--)
         sqrmod192(A,N,inv);

#if (PRE_SQUARE==0)
      return (A[0] == T[0] && A[1] == T[1] && A[2] == T[2]);
#else
      mulmod192(A,ONE,A,N,inv);
      CheckRedc(A,3,2,k);
      return (A[0] == 0 && A[1] == 0 && A[2] == 0);
#endif
   }
}

void  DMDivisorApp::CheckRedc(mp_limb_t *xp, uint32_t xn, uint32_t b, uint64_t k)
{
#ifndef TEST_REDC
   return;
#endif

   mpz_t     B, E, N;
   uint32_t  i;

   mpz_init_set_ui(B,b);

   mpz_init_set_ui(E,1);
   mpz_mul_2exp(E,E,ii_N);

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

   mpz_mul_2exp(N,N, ii_N);
   mpz_add_ui(N,N,1);

   mpz_powm(B,B,E,N);

   if (mpz_size(B) > xn)
   {
      WriteToConsole(COT_OTHER, "REDC ERROR: %u^2^%u mod %" PRIu64"*2^%u+1\n", b, ii_N, k, ii_N);
      
      FatalError("mpz_size=%" PRIu64", xn=%u\n", (uint64_t) mpz_size(B), xn);
   }

   for (i = 0; i < xn; i++)
      if (xp[i] != mpz_getlimbn(B,i))
      {        
         WriteToConsole(COT_OTHER, "REDC ERROR: %u^2^%u mod %" PRIu64"*2^%u+1\n", b, ii_N, k, ii_N);
         
         FatalError("limb %u=%" PRIu64", xp[%u]=%" PRIu64"\n", i, (uint64_t) mpz_getlimbn(B, i), i, (uint64_t) xp[i]);
      }

   mpz_clear(N);
   mpz_clear(E);
   mpz_clear(B);
}

void  DMDivisorApp::VerifyFactor(uint64_t theFactor, uint64_t k)
{
   uint64_t rem;

   fpu_push_1divp(theFactor);
      
   rem = fpu_powmod(2, ii_N, theFactor);
   
   if (rem == 0)
      rem = theFactor - 1;
   else
      rem--;
   
   rem = fpu_mulmod(rem, k, theFactor);

   rem <<= 1;
   rem++;
   
   if (rem >= theFactor)
      rem -= theFactor;
   
   fpu_pop();
   
   if (rem == 0)
      return;
   
   FatalError("Invalid factor: 2*%" PRIu64"*(2^%u-1)+1 mod %" PRIu64" = %" PRIu64"", k, ii_N, theFactor, rem);
}
