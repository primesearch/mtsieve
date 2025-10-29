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

#define APP_VERSION     "1.8.5"

// Arrays cannot have more than 2^32 elements on some systems.
// We will limit the vector to have at most 2^30 on all systems.
#define BITS_PER_ARRAY     30
#define TERMS_PER_VECTOR   (1 << BITS_PER_ARRAY)

// But we can have a vector of vectors, so that is what we are going
// to do to get beyond the 2^32 limit for an individual array.

// Only k % 4 = 0 and k % 4 = 1 can be factors
// Even k are tracked in iv_MMPTermsKm40.
// Odd k are tracked in iv_MMPTermsKm41.
#define IDX_K(k)           ((((((k) & 1) ? (k)-1 : (k)) - il_TermsMinK) / 4) >> BITS_PER_ARRAY)
#define BIT_K(k)           ((((((k) & 1) ? (k)-1 : (k)) - il_TermsMinK) / 4) & (TERMS_PER_VECTOR - 1))

#define MERSENNE_PRIMES    52

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
   {    19937,     4,     0,    1},
   {    21701,     2,     0,    1},
   {    23209,     2,     0,    1},

// Potential divisors for MMs here and below should be sieved with dmdsieve, but PRP tested
// with external software (such as pfgw) before doing any MM divisibility testing with pfgw.
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
   
   it_Format = FF_ABCD;
   SetAppMinPrime(3);
   SetAppMaxPrime(1ULL<<63);
   il_MinK = 0;
   il_MaxK = 0;
   ii_N    = 0;

   ib_TestTerms = false;
   
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
   printf("-f --format=f         Format of output file (A=ABC, D=ABCD (default))\n");
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
   AppendLongOpt(longOpts, "format",            required_argument, 0, 'f');
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
      case 'f':
         char value;
         status = Parser::Parse(arg, "AB", value);
         
         it_Format = FF_UNKNOWN;
   
         if (value == 'A')
            it_Format = FF_ABC;
         if (value == 'D')
            it_Format = FF_ABCD;
         break;
         

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
   uint32_t vCount, idx;

#ifdef WIN32
   if (sizeof(unsigned long long) != sizeof(mp_limb_t))
     FatalError("GMP limb size is not 64 bits");
#else
   if (sizeof(unsigned long) != sizeof(mp_limb_t))
     FatalError("GMP limb size is not 64 bits");
#endif

   if (it_Format == FF_UNKNOWN)
      FatalError("the specified file format in not valid, use A (ABC) or D (ABCD)");

   if (is_InputTermsFileName.length() > 0)
   {
      ProcessInputTermsFile(false);

      // We need this to be even
      il_TermsMinK = (il_MinK & 1) ? il_MinK - 1 : il_MinK;
      
      vCount = 1 + (1 + il_MaxK - il_TermsMinK) / TERMS_PER_VECTOR;

      iv_MMPTermsKm40.resize(vCount);
      iv_MMPTermsKm41.resize(vCount);
      for (idx=0; idx<vCount; idx++)
      {
         iv_MMPTermsKm40[idx].resize(20 + TERMS_PER_VECTOR);
         iv_MMPTermsKm41[idx].resize(20 + TERMS_PER_VECTOR);
         
         std::fill(iv_MMPTermsKm40[idx].begin(), iv_MMPTermsKm40[idx].end(), false);
         std::fill(iv_MMPTermsKm41[idx].begin(), iv_MMPTermsKm41[idx].end(), false);
      }
      
      ProcessInputTermsFile(true);
      
      if (GetMinPrime() >= GetMaxPrime() && ib_TestTerms)
      {
         WriteToConsole(COT_OTHER, "Cannot continue sieving due to p_min >= p_max");
         PostSieveHook();
         exit(0);
      }
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

      // We need this to be even
      il_TermsMinK = (il_MinK & 1) ? il_MinK - 1 : il_MinK;

      vCount = 1 + (1 + il_MaxK - il_TermsMinK) / TERMS_PER_VECTOR;

      iv_MMPTermsKm40.resize(vCount);
      iv_MMPTermsKm41.resize(vCount);
      for (uint32_t idx=0; idx<vCount; idx++)
      {
         iv_MMPTermsKm40[idx].resize(20 + TERMS_PER_VECTOR);
         iv_MMPTermsKm41[idx].resize(20 + TERMS_PER_VECTOR);
         
         std::fill(iv_MMPTermsKm40[idx].begin(), iv_MMPTermsKm40[idx].end(), false);
         std::fill(iv_MMPTermsKm41[idx].begin(), iv_MMPTermsKm41[idx].end(), false);
      }

      // Let d be a prime divisor of MMp. We know that
      //
      //    A factor d will always be of the form 2*k*m + 1 where m = M(p) = 2p - 1 is a Mersenne prime.
      //   
      // Note that, since m == 1 (mod 3), factors of M(m) cannot have k == 1 (mod 3),
      // since 2*k*m + 1 == 0 (mod 3) in that case.
      // Chris Nash pointed out on the Mersenne list that k must be 0 or 1 mod 4,
      // since, otherwise, 2 is not a quadratic residue of the supposed factor.
      // The combination of the prior two restrictions limits k to 0, 5, 8, and 9 modulo 12.
      uint64_t k;

      for (k=il_MinK; k%12>0; k++)
      {
         uint64_t km12 = k % 12;
         uint64_t bit = BIT_K(k);
         uint64_t idx = IDX_K(k);

         if (km12 == 0 || km12 == 8)
         {
            iv_MMPTermsKm40[idx][bit] = true;
            il_TermCount++;
         }

         if (km12 == 5 || km12 == 9)
         {
            iv_MMPTermsKm41[idx][bit] = true;
            il_TermCount++;
         }
      }

      for (; k<il_MaxK-24; k+=12)
      {         
         uint64_t bit0 = BIT_K(k+0);
         uint64_t idx0 = IDX_K(k+0);
         uint64_t bit5 = BIT_K(k+5);
         uint64_t idx5 = IDX_K(k+5);
         uint64_t bit8 = BIT_K(k+8);
         uint64_t idx8 = IDX_K(k+8);
         uint64_t bit9 = BIT_K(k+9);
         uint64_t idx9 = IDX_K(k+9);

         // For terms where k % 4 = 0
         iv_MMPTermsKm40[idx0][bit0] = true;      // for k = 0 % 12
         iv_MMPTermsKm40[idx8][bit8] = true;      // for k = 8 % 12
        
         // For terms where k % 4 = 1
         iv_MMPTermsKm41[idx5][bit5] = true;      // for k = 5 % 12
         iv_MMPTermsKm41[idx9][bit9] = true;      // for k = 9 % 12
         
         il_TermCount += 4;
      }
      
      for (; k<=il_MaxK; k++)
      {
         uint64_t km12 = k % 12;
         uint64_t bit = BIT_K(k);
         uint64_t idx = IDX_K(k);

         if (km12 == 0 || km12 == 8)
         {
            iv_MMPTermsKm40[idx][bit] = true;
            il_TermCount++;
         }

         if (km12 == 5 || km12 == 9)
         {
            iv_MMPTermsKm41[idx][bit] = true;
            il_TermCount++;
         }
      }
   }
            
   if (is_OutputTermsFileName.length() == 0)
   {
      char fileName[30];
      
      if (it_Format == FF_ABC)
         snprintf(fileName, sizeof(fileName), "mmp_%u.pfgw", ii_N);
      else
         snprintf(fileName, sizeof(fileName), "mmp_%u.abcd", ii_N);
      
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

   if (ib_TestTerms && ii_N > 10000)
      WriteToConsole(COT_OTHER, "It is recommended to PRP test remaining candidates before testing for DM divisibility with other software");

   // Since the worker wants primes in groups of 4
   while (ii_CpuWorkSize % 4 != 0)
      ii_CpuWorkSize++;
   
   // The GPU kernel doesn't have separate small k logic
   SetMinGpuPrime(100000);

   if (ii_N <= 19)
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
            if (dmdList[i].factorsPerSecond > 0 && id_FPSTarget == 0.0 && id_SPFTarget == 0.0)
            {
               id_FPSTarget = (double) dmdList[i].factorsPerSecond;     
               
               // At this time FactorApp only supports 1 minute averages for factors per second.
               // We won't change it because FactorApp won't work correctly if ii_MinutesForStatus 
               // is modified when computing factors per second.
            }

            if (dmdList[i].secodsPerFactor > 0 && id_SPFTarget == 0.0 && id_FPSTarget == 0.0)
            {
               id_SPFTarget = (double) dmdList[i].secodsPerFactor;
               
               if (ii_MinutesForStatus == MAX_FACTOR_REPORT_COUNT)
                  ii_MinutesForStatus = dmdList[i].minutesForAverage;
            }
         }
      }
   }

   if (id_FPSTarget > 0.0)               
      WriteToConsole(COT_OTHER, "Will stop sieving after averaging fewer than %.2lf factors per second for %u minute",
         id_FPSTarget, 1);
   
   if (id_SPFTarget > 0.0)
      WriteToConsole(COT_OTHER, "Will stop sieving after averaging more than %.2lf second per factor for %u minutes",
         id_SPFTarget, ii_MinutesForStatus);

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
   bool     isAbcd = false;
   uint32_t n;
   uint64_t k, prevk, lastPrime = 0;
   
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
      isAbcd = true;

      if (haveBitMap)
      {
         k = il_MinK;

         if (k % 4 > 1)
            FatalError("Bad k %llu", k);
         
         uint64_t bit = BIT_K(k);
         uint64_t idx = IDX_K(k); 
         
         if (k & 1)
            iv_MMPTermsKm41[idx][bit] = true;
         else
            iv_MMPTermsKm40[idx][bit] = true;
         
         il_TermCount++;
      }
      else
         il_MinK = il_MaxK = k;

      prevk = k;
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
  
      if (isAbcd)
      {
         k = prevk + k;
         prevk = k;
      }
      
      if (k % 4 > 1)
         FatalError("Bad k %llu.  Expecting k mod 4 < 2", k);
      
      if (haveBitMap)
      {
         uint64_t bit = BIT_K(k);
         uint64_t idx = IDX_K(k); 
         
         if (k & 1)
            iv_MMPTermsKm41[idx][bit] = true;
         else
            iv_MMPTermsKm40[idx][bit] = true;
         
         il_TermCount++;
      }
      else
      {
         if (il_MinK == 0)
            il_MinK = il_MaxK = k;
         
         if (k > il_MaxK)
            il_MaxK = k;
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

   if (k % 4 > 1)
      return false;

   uint64_t bit = BIT_K(k);
   uint64_t idx = IDX_K(k); 
         
   // No locking is needed because the Workers aren't running yet
   if (k & 1)
   {
      if (iv_MMPTermsKm41[idx][bit])
      {
         iv_MMPTermsKm41[idx][bit] = false;
         il_TermCount--;

         return true;
      }
   }
   else
   {
      if (iv_MMPTermsKm40[idx][bit])
      {
         iv_MMPTermsKm40[idx][bit] = false;
         il_TermCount--;

         return true;
      }
   }
      
   return false;
}

void DMDivisorApp::WriteOutputTermsFile(uint64_t largestPrime)
{
   uint64_t termsCounted = 0;
   uint64_t prevk = 0;

   if (ib_TestTerms)
      return;
   
   // With super large ranges, wait until we can lock because without locking
   // the term count can change between opening and closing the file.
   if (IsRunning() && largestPrime < GetMaxPrimeForSingleWorker())
      return;   
   
   // If we have no terms, then we have nothing to write
   if (il_TermCount == 0)
      return;
   
   ip_FactorAppLock->Lock();
   
   FILE    *termsFile = fopen(is_OutputTermsFileName.c_str(), "w");

   if (!termsFile)
      FatalError("Unable to open output file %s", is_OutputTermsFileName.c_str());

   if (it_Format == FF_ABC)
      fprintf(termsFile, "ABC 2*$a*(2^%u-1)+1 // Sieved to %" SCNu64"\n", ii_N, largestPrime);

   for (uint64_t k=il_MinK; k<=il_MaxK; k++)
   {    
      uint64_t km4 = k % 4;

      if (km4 > 1)
         continue;

      uint64_t bit = BIT_K(k);
      uint64_t idx = IDX_K(k); 

      if (km4 == 0 && !iv_MMPTermsKm40[idx][bit])
         continue;

      if (km4 == 1 && !iv_MMPTermsKm41[idx][bit])
         continue;
      
      if (it_Format == FF_ABC)
         fprintf(termsFile, "%" PRIu64"\n", k);
      else
      {
         if (termsCounted == 0)
            fprintf(termsFile, "ABCD 2*$a*(2^%u-1)+1 [%" SCNu64"] // Sieved to %" SCNu64"\n", ii_N, k, largestPrime);
         else
            fprintf(termsFile, "%" PRIu64"\n", k - prevk);
         
         prevk = k;
      }

      termsCounted++;
   }
   
   fclose(termsFile);
   
   if (termsCounted != il_TermCount)
      FatalError("Something is wrong.  Counted terms (%" PRIu64") != expected terms (%" PRIu64")", termsCounted, il_TermCount);

   ip_FactorAppLock->Release();
}

void  DMDivisorApp::GetExtraTextForSieveStartedMessage(char *extraText, uint32_t maxTextLength)
{
   snprintf(extraText, maxTextLength, "%" PRIu64 " <= k <= %" PRIu64", 2*k*(2^%u-1)+1", il_MinK, il_MaxK, ii_N);
}

void  DMDivisorApp::ReportFactor(uint64_t theFactor, uint64_t k)
{
   uint32_t verifiedCount = 0;
      
   if (theFactor > GetMaxPrimeForSingleWorker())
      ip_FactorAppLock->Lock();

   if (ii_N < 62)
   {
      if (ii_N == 13 && theFactor == 338193759479) return;
      if (ii_N == 17 && theFactor == 231733529) return;
      if (ii_N == 17 && theFactor == 64296354767) return;
      if (ii_N == 19 && theFactor == 62914441) return;
      if (ii_N == 19 && theFactor == 5746991873407) return;
      if (ii_N == 31 && theFactor == 295257526626031) return;
      if (ii_N == 61 && theFactor == 2305843009213693951) return;
   }

   // Only k where k%4 < 2 are considered.
   while (k <= il_MaxK)
   {
      uint64_t km4 = k % 4;
      
      if (km4 < 2)
      {
         uint64_t bit = BIT_K(k);
         uint64_t idx = IDX_K(k); 

         if ((km4 == 0 && iv_MMPTermsKm40[idx][bit]) || (km4 == 1 && iv_MMPTermsKm41[idx][bit]))
         {
            // We only need to verify the first two
            if (++verifiedCount < 3)
               VerifyFactor(theFactor, k);
            
            if (km4 == 0)
               iv_MMPTermsKm40[idx][bit] = false;
            else
               iv_MMPTermsKm41[idx][bit] = false;

            LogFactor(theFactor, "2*%" PRIu64"*(2^%u-1)+1", k, ii_N);
            
            il_FactorCount++;
            il_TermCount--;
         }
      }
      
		k += theFactor; 
	};

   if (theFactor > GetMaxPrimeForSingleWorker())
      ip_FactorAppLock->Release();
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

   FILE *tPtr = fopen("tested.out", "w");
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
      kEvaluated++;
      
      uint64_t km4 = k % 4;

      if (km4 > 1)
         continue;

      uint64_t bit = BIT_K(k);
      uint64_t idx = IDX_K(k); 

      if (km4 == 0 && !iv_MMPTermsKm40[idx][bit])
         continue;

      if (km4 == 1 && !iv_MMPTermsKm41[idx][bit])
         continue;

      currentTime = time(NULL);
      kTested++;
      
         fprintf(tPtr, "%" PRIu64"\n", k);
         
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
            WriteToConsole(COT_SIEVE, " Tested %5.2f pct of range at %" PRIu64" tests per second (%4.2f pct of all k after sieving) ETC %s", 
                           percentCompleted * 100.0, kTestedPerSecond, percentTermsRequiringTest, finishTimeBuffer);
         else 
         {
            kSecondsPerTest = (currentTime - startTime) / kTested;

            WriteToConsole(COT_SIEVE, " Tested %5.2f pct of range at %u seconds per test (%4.2f pct of all k after sieving) ETC %s", 
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
   
   fclose(tPtr);
   
   if (il_TermCount != kTested)
      FatalError("Expected to test %" PRIu64" terms, but tested %" PRIu64" terms                                                  ", il_TermCount, kTested);
   
   percentTermsRequiringTest = (100.0 * (double) kTested) / (double) kEvaluated;
   
   if (currentTime > startTime)
   {
      kTestedPerSecond = kTested / (currentTime - startTime);
   
      if (kTestedPerSecond >= 1)
         WriteToConsole(COT_SIEVE, "Tested %" PRIu64" terms for divisilibity at %" PRIu64" tests per second (%4.2f pct of all k after sieving)                ", 
                        kTested, kTestedPerSecond, percentTermsRequiringTest);
      else 
      {
         kSecondsPerTest = (currentTime - startTime) / kTested;
               
         WriteToConsole(COT_SIEVE, "Tested %" PRIu64" terms for divisilibity at %u seconds per test (%4.2f pct of all k after sieving)                ", 
                        kTested, kSecondsPerTest, percentTermsRequiringTest);
      }
   }

   WriteToConsole(COT_OTHER, "Double-Mersenne divisibility checking completed in %" PRIu64" seconds.  %u factors found                                ",
                  (currentTime - startTime), factorsFound);
   
   return true;
}