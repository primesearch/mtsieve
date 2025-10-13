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
#include "../core/MpArith.h"
#include "DMDivisorApp.h"
#include "DMDivisorWorker.h"

#if defined(USE_OPENCL) || defined(USE_METAL)
#include "DMDivisorGpuWorker.h"
#define APP_NAME        "dmdsievecl"
#else
#define APP_NAME        "dmdsieve"
#endif

#define APP_VERSION     "1.6"

#define BIT0(k)         (((k) - il_MinK) / 4)
#define BIT1(k)         (((k-1) - il_MinK) / 4)


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
   il_MinK = 0;
   il_MaxK = 0;
   ii_N    = 0;
      
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
   
#if defined(USE_OPENCL) || defined(USE_METAL)
   printf("-M --maxfactors=M        max number of factors to support per GPU worker chunk (default %u)\n", ii_MaxGpuFactors);
#endif
}

void  DMDivisorApp::AddCommandLineOptions(std::string &shortOpts, struct option *longOpts)
{
   FactorApp::ParentAddCommandLineOptions(shortOpts, longOpts);

   shortOpts += "k:K:b:n:f:";

   AppendLongOpt(longOpts, "kmin",              required_argument, 0, 'k');
   AppendLongOpt(longOpts, "kmax",              required_argument, 0, 'K');
   AppendLongOpt(longOpts, "exp",               required_argument, 0, 'n');
   
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
   uint32_t mList[] = { 2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 
                        607, 1279, 2203, 2281, 3217, 4253, 4423, 9689, 9941, 
                        11213, 19937, 21701, 23209, 44497, 86243, 110503, 
                        132049, 216091, 756839, 859433, 1257787, 1398269, 
                        2976221, 3021377, 6972593, 13466917, 20996011, 24036583, 
                        25964951, 30402457, 32582657, 37156667, 42643801, 43112609, 
                        57885161, 74207281, 77232917, 82589933, 136279841, 0};

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
         uint32_t bit = BIT0(k);
         
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
   while (mList[i] != 0 && mList[i] != ii_N)
      i++;
   
   if (mList[i] == 0)
   {
      WriteToConsole(COT_OTHER, "Known Mersenne Primes (as of 2025) are:");
      for (i=0; mList[i+1] > 0; i++)
      {
         printf(" %u,", mList[i]);
         if (i > 0 && i % 8 == 0)
            printf("\n");
      }
      printf(" %u\n", mList[i]);
      
      FatalError("exponent must be for a Mersenne Prime");
   }

   if (ii_N < 13)
      FatalError("MM%u is a known prime", ii_N);
 
   FactorApp::ParentValidateOptions();

   // Since the worker wants primes in groups of 4
   while (ii_CpuWorkSize % 4 != 0)
      ii_CpuWorkSize++;
   
   // The GPU kernel doesn't have separate small k logic
   SetMinGpuPrime(il_MaxK);

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
            bit = BIT0(k);
            iv_MMPTerms0[bit] = true;
         }
         else
         {
            bit = BIT1(k);
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
            bit = BIT0(k);
            iv_MMPTerms0[bit] = true;
         }
         else
         {
            bit = BIT1(k);
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
      bit = BIT0(k);
   
      if (iv_MMPTerms0[bit])
      {
         iv_MMPTerms0[bit] = false;
         il_TermCount--;

         return true;
      }
   }
   else
   {
      bit = BIT1(k);
      
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

   // With super large ranges, wait until we can lock because without locking
   // the term count can change between opening and closing the file.
   if (IsRunning() && largestPrime < GetMaxPrimeForSingleWorker())
      return;   
   
   for (k=il_MinK; k<=il_MaxK; k+=4)
   {
      bit = BIT0(k);
      
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
      if (iv_MMPTerms0[BIT0(k)])
      {
         fprintf(termsFile, "%" PRIu64"\n", k);
         termsCounted++;
      }
 
      if (iv_MMPTerms1[BIT1(k+1)])
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

   uint64_t bit = ((km4 == 0) ? BIT0(k) : BIT1(k));

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
