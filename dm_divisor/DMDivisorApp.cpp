/* DMDivisorApp.cpp -- (C) Mark Rodenkirch, September 2018

   Sieve for 2*k*M+1 where M = is a Mersenne Prime.  The remaining terms,
   if prime, might be divisors of 2^(2^n-1)-1.

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
#include "../core/MpArith.h"
#include "DMDivisorApp.h"
#include "DMDivisorWorker.h"

#if defined(USE_OPENCL) || defined(USE_METAL)
#include "DMDivisorGpuWorker.h"
#define APP_NAME        "dmdsievecl"
#else
#define APP_NAME        "dmdsieve"
#endif

#define APP_VERSION     "1.7"

#define MERSENNE_PRIMES 52

#define BIT_EK(k)         (((k) - il_MinK) / 4)
#define BIT_OK(k)         (((k-1) - il_MinK) / 4)

typedef struct {
   uint32_t n;
   uint32_t factorsPerSecond;
   uint32_t secodsPerFactor;
   uint32_t minutesForAverage;
} dmd_t;

// Eventually this program will calculate optimal removal rate on the fly
// rather than being hard-coded in this table.
dmd_t  dmdList[MERSENNE_PRIMES] = {
   {        2,     0,     0,    0},
   {        3,     0,     0,    0},
   {        5,     0,     0,    0},
   {        7,     0,     0,    0},
   {       13,     0,     0,    0},
   {       17,     0,     0,    0},
   {       19,     0,     0,    0},
   {       31,     0,     0,    0},
   {       61,     0,     0,    0},
   {       89,     0,     0,    0},
   {      107,     0,     0,    0},
   {      127,     0,     0,    0}, 
   {      521, 20000,     0,    1},
   {      607,  9000,     0,    1},
   {     1279,  1800,     0,    1},
   {     2203,   500,     0,    1},
   {     2281,   400,     0,    1},
   {     3217,   200,     0,    1},
   {     4253,    80,     0,    1},
   {     4423,    70,     0,    1},
   {     9689,    10,     0,    1},
   {     9941,     9,     0,    1},
   {    11213,     8,     0,    1},
   {    19937,     1,     0,    1},
   {    21701,     1,     0,    1},

// Potential divisors for MMs here and below should be sieved with dmdsieve, but PRP tested
// with external software (such as pfgw) before doing any MM divisibility testing with pfgw.
   {    23209,     0,     4,    5},
   {    44497,     0,     7,    5},
   {    86243,     0,    10,    5},
   {   110503,     0,    30,    6}, 
   {   132049,     0,    50,    7},
   {   216091,     0,    80,   10},
   {   756839,     0,   250,  100},
   {   859433,     0,   400,  200},
   {  1257787,     0,   800,  400},
   {  1398269,     0,   950,  500},
   {  2976221,     0,  2000, 1000},
   {  3021377,     0,     0,    0},
   {  6972593,     0,     0,    0},
   { 13466917,     0,     0,    0},
   { 20996011,     0,     0,    0},
   { 24036583,     0,     0,    0},
   { 25964951,     0,     0,    0},
   { 30402457,     0,     0,    0},
   { 32582657,     0,     0,    0},
   { 37156667,     0,     0,    0},
   { 42643801,     0,     0,    0},
   { 43112609,     0,     0,    0},
   { 57885161,     0,     0,    0},
   { 74207281,     0,     0,    0},
   { 77232917,     0,     0,    0},
   { 82589933,     0,     0,    0},
   {136279841,     0,     0,    0}
};

// This is declared in App.h, but implemented here.  This means that App.h
// can remain unchanged if using the mtsieve framework for other applications.
App *get_app(void)
{
   return new DMDivisorApp();
}

DMDivisorApp::DMDivisorApp() : FactorApp()
{
   SetBanner(APP_NAME " v" APP_VERSION ", a program to find terms that are potential factors of 2^(2^n-1)-1");
   SetLogFileName("dmdsieve.log");
   
   SetAppMinPrime(3);
   il_MinK = 0;
   il_MaxK = 0;
   ii_N    = 0;

   ib_TestTerms = false;
   
   iv_MMPTerms0.clear();
   iv_MMPTerms1.clear();

#if defined(USE_OPENCL) || defined(USE_METAL)
   ii_MaxGpuFactors = GetGpuWorkGroups() * 1000;
#endif
}

void DMDivisorApp::Help(void)
{
   FactorApp::ParentHelp();

   printf("-k --kmin=k           Minimum k to search\n");
   printf("-K --kmax=K           Maximum k to search\n");
   printf("-n --exp=n            Exponent to search\n");
   printf("-x --testterms        test remaining terms for DM divisibility\n");
   
#if defined(USE_OPENCL) || defined(USE_METAL)
   printf("-M --maxfactors=M        max number of factors to support per GPU worker chunk (default %u)\n", ii_MaxGpuFactors);
#endif
}

void  DMDivisorApp::AddCommandLineOptions(std::string &shortOpts, struct option *longOpts)
{
   FactorApp::ParentAddCommandLineOptions(shortOpts, longOpts);

   shortOpts += "k:K:b:n:f:x";

   AppendLongOpt(longOpts, "kmin",              required_argument, 0, 'k');
   AppendLongOpt(longOpts, "kmax",              required_argument, 0, 'K');
   AppendLongOpt(longOpts, "exp",               required_argument, 0, 'n');
   AppendLongOpt(longOpts, "testterms",         no_argument,       0, 'x');
   
#if defined(USE_OPENCL) || defined(USE_METAL)
   shortOpts += "M:";
   
   AppendLongOpt(longOpts, "maxfactors",        required_argument, 0, 'M');
#endif
}

parse_t DMDivisorApp::ParseOption(int opt, char *arg, const char *source)
{
   parse_t status = P_UNSUPPORTED;

   status = FactorApp::ParentParseOption(opt, arg, source);
   if (status != P_UNSUPPORTED) return status;

   switch (opt)
   {
      case 'k':
         status = Parser::Parse(arg, 4, KMAX_MAX, il_MinK);
         break;

      case 'K':
         status = Parser::Parse(arg, 4, KMAX_MAX, il_MaxK);
         break;

      case 'n':
         status = Parser::Parse(arg, 1, NMAX_MAX, ii_N);
         break;

      case 'x':
         ib_TestTerms = true;
         status = P_SUCCESS;
         break;
         
#if defined(USE_OPENCL) || defined(USE_METAL)
      case 'M':
         status = Parser::Parse(arg, 1, 100000000, ii_MaxGpuFactors);
         break;
#endif
   }

   return status;
}

void DMDivisorApp::ValidateOptions(void)
{
#ifdef WIN32
   if (sizeof(unsigned long long) != sizeof(mp_limb_t))
     FatalError("GMP limb size is not 64 bits");
#else
   if (sizeof(unsigned long) != sizeof(mp_limb_t))
     FatalError("GMP limb size is not 64 bits");
#endif

   if (is_InputTermsFileName.length() > 0)
   {
      ProcessInputTermsFile(false);

      iv_MMPTerms0.resize((il_MaxK - il_MinK) / 4 + 1);
      iv_MMPTerms1.resize((il_MaxK - il_MinK) / 4 + 1);
      std::fill(iv_MMPTerms0.begin(), iv_MMPTerms0.end(), false);
      std::fill(iv_MMPTerms1.begin(), iv_MMPTerms1.end(), false);
      
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

      // We want kmin where kmin % 4 == 0
      while (il_MinK % 4 > 1)
         il_MinK++;
      
      while (il_MinK % 4 != 0)
         il_MinK--;
      
      // We want kmin where kmax % 4 == 1
      while (il_MaxK % 4 != 1)
         il_MaxK++;

      iv_MMPTerms0.resize((il_MaxK - il_MinK) / 4 + 1);
      iv_MMPTerms1.resize((il_MaxK - il_MinK) / 4 + 1);
      std::fill(iv_MMPTerms0.begin(), iv_MMPTerms0.end(), false);
      std::fill(iv_MMPTerms1.begin(), iv_MMPTerms1.end(), false);
      
      for (uint64_t k=il_MinK; k<=il_MaxK; k+=4)
      {
         uint32_t bit = BIT_EK(k);
         
         // Only terms where k % 4 == 0 or k % 4 == 1 can be a factor of a double-mersenne
         iv_MMPTerms0[bit] = true;
         iv_MMPTerms1[bit] = true;
         
         il_TermCount += 2;
      }
   }
            
   if (is_OutputTermsFileName.length() == 0)
   {
      char fileName[30];
      
      snprintf(fileName, sizeof(fileName), "mmp_%u.pfgw", ii_N);
      
      is_OutputTermsFileName = fileName;
   }
      
   uint32_t i = 0;
   while (i != MERSENNE_PRIMES && dmdList[i].n != ii_N)
      i++;
   
   if (i == MERSENNE_PRIMES)
   {
      WriteToConsole(COT_OTHER, "Known Mersenne Primes (as of 2025) are:");
      for (i=0; i<MERSENNE_PRIMES-1; i++)
      {
         printf(" %u,", dmdList[i].n);
         if (i > 0 && i % 8 == 0)
            printf("\n");
      }
      printf(" %u\n", dmdList[i].n);
      
      FatalError("exponent must be for a Mersenne Prime");
   }

   if (ii_N <= 7)
      FatalError("MM%u is a known prime", ii_N);
   else if (ii_N < 30)
      WriteToConsole(COT_OTHER, "MM%u should be factored with ECM", ii_N);
   else if (ii_N < 130)
      WriteToConsole(COT_OTHER, "MM%u should be tested with mmff", ii_N);
 
   FactorApp::ParentValidateOptions();

   if (ib_TestTerms && ii_N > 100000)
      WriteToConsole(COT_OTHER, "It is recommended to use PRP test remaining candidates before testing for DM divisibility");

   // Since the worker wants primes in groups of 4
   while (ii_CpuWorkSize % 4 != 0)
      ii_CpuWorkSize++;
   
   // The GPU kernel doesn't have separate small k logic
   SetMinGpuPrime(il_MaxK);

   if (ii_N <= 31)
   {
      double sqrtK = sqrt(il_MaxK);
      SetMaxPrime((uint64_t) sqrtK * (2 << (ii_N+1)), "sqrt of largest term");
   }

   for (uint32_t i=0; i<MERSENNE_PRIMES; i++)
   {
      if (dmdList[i].n == ii_N)
      {
         if (dmdList[i].minutesForAverage > 0)
         {
            if (dmdList[i].factorsPerSecond > 0 && id_FPSTarget == 0.0)
            {
               id_FPSTarget = (double) dmdList[i].factorsPerSecond;
               
               // At this time FactorApp only supports 1 minute averages for factors per second.
               // We won't change it because FactorApp won't work correctly if ii_MinutesForStatus 
               // is modified when computing factors per second.
               
               WriteToConsole(COT_OTHER, "Will stop sieving after averaging fewer than %d factors per second for %u minute",
                  dmdList[i].factorsPerSecond, 1);
            }

            if (dmdList[i].secodsPerFactor > 0 && id_SPFTarget == 0.0)
            {
               id_SPFTarget = (double) dmdList[i].secodsPerFactor;
               
               if (ii_MinutesForStatus == MAX_FACTOR_REPORT_COUNT)
                  ii_MinutesForStatus = dmdList[i].minutesForAverage;
            
               WriteToConsole(COT_OTHER, "Will stop sieving after averaging more than %d second per factor for %u minutes",
                  dmdList[i].secodsPerFactor, ii_MinutesForStatus);
            }
         }
      }
   }

   // Allow only one worker to do work when processing small primes.  This allows us to avoid 
   // locking when factors are reported, which significantly hurts performance as most terms 
   // will be removed due to small primes.
   SetMaxPrimeForSingleWorker(10000);
}

Worker *DMDivisorApp::CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested)
{
   Worker *theWorker;

#if defined(USE_OPENCL) || defined(USE_METAL)
   if (gpuWorker)
      return new DMDivisorGpuWorker(id, this);
#endif

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
         k = il_MinK;
         
         while (il_MinK % 4 != 0)
            il_MinK--;

         if (k % 4 == 0)
         {
            bit = BIT_EK(k);
            iv_MMPTerms0[bit] = true;
         }
         else
         {
            bit = BIT_OK(k);
            iv_MMPTerms1[bit] = true;
         }
         
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
      
         if (k % 4 == 0)
         {
            bit = BIT_EK(k);
            iv_MMPTerms0[bit] = true;
         }
         else
         {
            bit = BIT_OK(k);
            iv_MMPTerms1[bit] = true;
         }
         
         il_TermCount++;
      }
      else
      {
         if (il_MinK > k || il_MinK == 0)
         {
            while (k % 4 != 0)
               k--;
            
            il_MinK = k;
         }
         
         if (il_MaxK < k) il_MaxK = k;
      }
   }

   fclose(fPtr);
}

bool DMDivisorApp::ApplyFactor(uint64_t theFactor, const char *term)
{
   uint64_t k, bit;
   uint32_t n;
      
   if (sscanf(term, "2*%" SCNu64"*(2^%u-1)+1", &k, &n) != 2)
      FatalError("Could not parse term %s", term);

   if (n != ii_N)
      FatalError("Expected n %u in factor but found %d", ii_N, n);
   
   if (k < il_MinK || k > il_MaxK)
      return false;

   VerifyFactor(theFactor, k);
   
   // No locking is needed because the Workers aren't running yet
   if (k % 4 == 0)
   {
      bit = BIT_EK(k);
   
      if (iv_MMPTerms0[bit])
      {
         iv_MMPTerms0[bit] = false;
         il_TermCount--;

         return true;
      }
   }
   else
   {
      bit = BIT_OK(k);
      
      if (iv_MMPTerms1[bit])
      {
         iv_MMPTerms1[bit] = false;
         il_TermCount--;

         return true;
      }
   }
      
   return false;
}

void DMDivisorApp::WriteOutputTermsFile(uint64_t largestPrime)
{
   uint64_t termsCounted = 0;
   uint64_t k;
   uint64_t bit;

   if (ib_TestTerms)
      return;
   
   // With super large ranges, wait until we can lock because without locking
   // the term count can change between opening and closing the file.
   if (IsRunning() && largestPrime < GetMaxPrimeForSingleWorker())
      return;   
   
   for (k=il_MinK; k<=il_MaxK; k+=4)
   {
      bit = BIT_EK(k);
      
      if (iv_MMPTerms0[bit] || iv_MMPTerms1[bit])
         break;
   }

   if (k > il_MaxK)
      return;
   
   FILE    *termsFile = fopen(is_OutputTermsFileName.c_str(), "w");

   if (!termsFile)
      FatalError("Unable to open output file %s", is_OutputTermsFileName.c_str());
   
   ip_FactorAppLock->Lock();

   fprintf(termsFile, "ABC 2*$a*(2^%u-1)+1 // Sieved to %" SCNu64"\n", ii_N, largestPrime);
   // There is no logic for writing ABCD files at this time.
   //fprintf(termsFile, "ABCD 2*$a*(2^%u-1)+1 [%" SCNu64"] // Sieved to %" SCNu64"\n", ii_N, k, largestPrime);
   
   for (k=il_MinK; k<=il_MaxK; k+=4)
   {   
      if (iv_MMPTerms0[BIT_EK(k)])
      {
         fprintf(termsFile, "%" PRIu64"\n", k);
         termsCounted++;
      }
 
      if (iv_MMPTerms1[BIT_OK(k+1)])
      {
         fprintf(termsFile, "%" PRIu64"\n", k+1);
         termsCounted++;
      }
   }
   
   fclose(termsFile);
   
   if (termsCounted != il_TermCount)
      FatalError("Something is wrong.  Counted terms (%" PRIu64") != expected terms (%" PRIu64")", termsCounted, il_TermCount);

   ip_FactorAppLock->Release();
}

void  DMDivisorApp::GetExtraTextForSieveStartedMessage(char *extraText, uint32_t maxTextLength)
{
   snprintf(extraText, maxTextLength, "%" PRIu64 " < k < %" PRIu64", 2*k*(2^%u-1)+1", il_MinK, il_MaxK, ii_N);
}

bool  DMDivisorApp::ReportFactor(uint64_t theFactor, uint64_t k, bool verifyFactor)
{
   bool     removedTerm = false;
   char     kStr[50];

   if (theFactor > GetMaxPrimeForSingleWorker())
      ip_FactorAppLock->Lock();

   uint64_t km4 = k % 4;

   uint64_t bit = ((km4 == 0) ? BIT_EK(k) : BIT_OK(k));

   if ((km4 == 0 && iv_MMPTerms0[bit]) || (km4 == 1 && iv_MMPTerms1[bit]))
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
      
      if (km4 == 0)
         iv_MMPTerms0[bit] = false;
      else
         iv_MMPTerms1[bit] = false;

      removedTerm = true;
      
      snprintf(kStr, sizeof(kStr), "%" PRIu64"", k);
      
      LogFactor(theFactor, "2*%s*(2^%u-1)+1", kStr, ii_N);
      
      il_FactorCount++;
      il_TermCount--;
   }

   if (theFactor > GetMaxPrimeForSingleWorker())
      ip_FactorAppLock->Release();
   
   return removedTerm;
}

void  DMDivisorApp::VerifyFactor(uint64_t theFactor, uint64_t k)
{
   const MpArith mp(theFactor);

   MpRes res = mp.pow(mp.nToRes(2), ii_N);
   
   res = mp.sub(res, mp.one());
   res = mp.mul(res, mp.nToRes(2*k));
   
   uint64_t rem = mp.resToN(res);

   rem++;
   
   // Now we have 2*k*2^(n-1)+1   
   if (rem >= theFactor)
      rem -= theFactor;
   
   if (rem == 0)
      return;
   
   FatalError("Invalid factor: 2*%" PRIu64"*(2^%u-1)+1 mod %" PRIu64" = %" PRIu64"", k, ii_N, theFactor, rem);
}

bool  DMDivisorApp::PostSieveHook(void)
{
   if (!ib_TestTerms)
      return true;
   
   time_t   lastCheckPointTime;
   time_t   startTime, currentTime, estimatedFinishTime;
   uint64_t kEvaluated = 0, kTested = 0;
   uint32_t factorsFound = 0;
   uint64_t kTestedPerSecond;
   uint32_t kSecondsPerTest;
   double   percentCompleted;
   double   percentTermsRequiringTest;
   FILE    *fPtr;
   mpz_t    rem, mersenne, nTemp, kTemp, factor;
   char     finishTimeBuffer[40];
   struct tm   *finish_tm;

   WriteToConsole(COT_OTHER, "Starting Double-Mesenne divisibility checking");

   lastCheckPointTime = startTime = currentTime = time(NULL);
            
   mpz_init(rem);
   mpz_init(mersenne);
   mpz_init(nTemp);
   mpz_init(kTemp);
   mpz_init(factor);
   
   mpz_set_ui(nTemp, 2);
   mpz_pow_ui(nTemp, nTemp, ii_N);
   mpz_sub_ui(nTemp, nTemp, 1);
   mpz_set_ui(mersenne, 2);

   for (uint64_t k=il_MinK; k<=il_MaxK; k++)
   {
      uint64_t km4 = k % 4;

      if (km4 > 1)
         continue;

      kEvaluated++;

      if (km4 == 0 && !iv_MMPTerms0[BIT_EK(k)])
         continue;

      if (km4 == 1 && !iv_MMPTerms1[BIT_OK(k)])
         continue;

      currentTime = time(NULL);
      kTested++;
      
#ifdef WIN32
      // Even though built with 64-bit limbs, mpz_set_ui doesn't
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
         factorsFound++;
         WriteToConsole(COT_OTHER, "Found factor 2*%" PRIu64"*(2^%u-1)+1 of 2^(2^%u-1)-1", k, ii_N, ii_N);
         
         fPtr = fopen("dm_factors.txt", "a+");
         fprintf(fPtr, "Found factor 2*%" PRIu64"*(2^%u-1)+1 of 2^(2^%u-1)-1\n", k, ii_N, ii_N);
         fclose(fPtr);
      }

      // Report once every 10 seconds.
      if (currentTime >= lastCheckPointTime + 10)
      {            
         percentCompleted = ((double) (k - il_MinK)) / ((double) (il_MaxK - il_MinK));

         percentTermsRequiringTest = (100.0 * (double) kTested) / (double) kEvaluated;
         
         kTestedPerSecond = kTested / (currentTime - startTime);
         
         estimatedFinishTime = startTime + (currentTime - startTime)/percentCompleted;
         finish_tm = localtime(&estimatedFinishTime);
         strftime(finishTimeBuffer, sizeof(finishTimeBuffer), "%Y-%m-%d %H:%M", finish_tm);
         
         if (kTestedPerSecond >= 1)
            WriteToConsole(COT_SIEVE, " Tested %5.2f pct of range at %" PRIu64" tests per second (%4.2f pct terms passed sieving) ETC %s", 
                           percentCompleted * 100.0, kTestedPerSecond, percentTermsRequiringTest, finishTimeBuffer);
         else 
         {
            kSecondsPerTest = (currentTime - startTime) / kTested;

            WriteToConsole(COT_SIEVE, " Tested %5.2f pct of range at %u seconds per test (%4.2f pct terms passed sieving) ETC %s", 
                           percentCompleted * 100.0, kSecondsPerTest, percentTermsRequiringTest, finishTimeBuffer);
         }

         lastCheckPointTime = currentTime;
      }
   }

   mpz_clear(rem);
   mpz_clear(mersenne);
   mpz_clear(nTemp);
   mpz_clear(kTemp);
   mpz_clear(factor);
               
   percentTermsRequiringTest = (100.0 * (double) kTested) / (double) kEvaluated;
   
   kTestedPerSecond = kTested / (currentTime - startTime);
   
   if (currentTime - startTime > 0)
   {
      if (kTestedPerSecond >= 1)
         WriteToConsole(COT_SIEVE, "Tested %" PRIu64" terms for divisilibity at %" PRIu64" tests per second (%4.2f pct terms passed sieving)", 
                        kTested, kTestedPerSecond, percentTermsRequiringTest);
      else 
      {
         kSecondsPerTest = (currentTime - startTime) / kTested;
               
         WriteToConsole(COT_SIEVE, "Tested %" PRIu64" terms for divisilibity at %u seconds per test (%4.2f pct terms passed sieving)", 
                        kTested, kSecondsPerTest, percentTermsRequiringTest);
      }
   }

   WriteToConsole(COT_OTHER, "Double-Mersenne divisibility checking completed in %llu seconds.  %u factors found",
                  (currentTime - startTime),   factorsFound);
   
   return true;
}