/* KBBApp.cpp -- (C) Mark Rodenkirch, November 2017

   Sieve for k*2^n+1 for a range of k and a range of n

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <stdarg.h>
#include "../core/Parser.h"
#include "../core/Clock.h"
#include "KBBApp.h"
#include "KBBWorker.h"

#define APP_NAME        "kbbsieve"
#define APP_VERSION     "1.0"

// This is declared in App.h, but implemented here.  This means that App.h
// can remain unchanged if using the mtsieve framework for other applications.
App *get_app(void)
{
   return new KBBApp();
}

KBBApp::KBBApp() : AlgebraicFactorApp()
{
   SetBanner(APP_NAME " v" APP_VERSION ", a program to find factors of k*b^b+c numbers for fixed k and c = +1/-1");
   SetLogFileName("kbb.log");

   il_K = 0;
   ii_MinB = 0;
   ii_MaxB = 0;
   
   ii_CpuWorkSize = 10000;
}

void KBBApp::Help(void)
{
   FactorApp::ParentHelp();

   printf("-k ---k=k             k to search\n");
   printf("-b --bmin=b           Minimum b to search\n");
   printf("-B --bmax=B           Maximum B to search\n");
}

void  KBBApp::AddCommandLineOptions(string &shortOpts, struct option *longOpts)
{
   FactorApp::ParentAddCommandLineOptions(shortOpts, longOpts);

   shortOpts += "k:b:B:";

   AppendLongOpt(longOpts, "k",              required_argument, 0, 'k');
   AppendLongOpt(longOpts, "bmin",           required_argument, 0, 'b');
   AppendLongOpt(longOpts, "bmax",           required_argument, 0, 'B');
}

parse_t KBBApp::ParseOption(int opt, char *arg, const char *source)
{
   parse_t status = P_UNSUPPORTED;

   status = FactorApp::ParentParseOption(opt, arg, source);
   if (status != P_UNSUPPORTED) return status;

   switch (opt)
   {
      case 'k':
         status = Parser::Parse(arg, 2, PMAX_MAX_62BIT, il_K);
         break;

      case 'b':
         status = Parser::Parse(arg, 2, BMAX_MAX, ii_MinB);
         break;

      case 'B':
         status = Parser::Parse(arg, 3, BMAX_MAX, ii_MaxB);
         break;
   }

   return status;
}

void KBBApp::ValidateOptions(void)
{   
   if (is_InputTermsFileName.length() > 0)
   {
      ProcessInputTermsFile(false);
      
      iv_MinusTerms.resize(ii_MaxB - ii_MinB + 1);
      iv_PlusTerms.resize(ii_MaxB - ii_MinB + 1);
      
      std::fill(iv_MinusTerms.begin(), iv_MinusTerms.end(), false);
      std::fill(iv_PlusTerms.begin(), iv_PlusTerms.end(), false);
      
      ProcessInputTermsFile(true);
   }
   else
   {
      if (il_K == 0)
         FatalError("k must be specified");
      
      if (ii_MinB == 0)
         FatalError("bmin must be specified");

      if (ii_MaxB == 0)
         FatalError("bmax must be specified");
      
      if (ii_MaxB <= ii_MinB)
         FatalError("bmax must be greater than bmin");
 
      il_TermCount = 2 * (ii_MaxB - ii_MinB + 1);
      
      iv_MinusTerms.resize(ii_MaxB - ii_MinB + 1);
      iv_PlusTerms.resize(ii_MaxB - ii_MinB + 1);
      
      std::fill(iv_MinusTerms.begin(), iv_MinusTerms.end(), true);
      std::fill(iv_PlusTerms.begin(), iv_PlusTerms.end(), true);
      
      EliminateGfnAndMersenneTerms();
   }
   
   if (is_OutputTermsFileName.length() == 0)
   {
      char  fileName[30];
      
      sprintf(fileName, "k%" PRIu64"_bb.pfgw", il_K);
      is_OutputTermsFileName = fileName;
   } 
   
   FactorApp::ParentValidateOptions();

   // This is only done when starting a new sieve, but has to be done
   // after the call to ParentValidateOptions() as that method will 
   // open the factor file that algebraic factors will be written to.
   if (is_InputTermsFileName.length() == 0)
      EliminateTermsWithAlgebraicFactors();
      
   // Since the worker wants primes in groups of 4
   while (ii_CpuWorkSize % 4 != 0)
      ii_CpuWorkSize++;

   // Allow only one worker to do work when processing small primes.  This allows us to avoid 
   // locking when factors are reported, which significantly hurts performance as most terms 
   // will be removed due to small primes.
   SetMaxPrimeForSingleWorker(10000);
}

void KBBApp::EliminateGfnAndMersenneTerms(void)
{
   uint32_t b, bit;
   uint32_t removedCount = 0;
   
   for (b=ii_MinB; b<=ii_MaxB; b++)
   {
      bit = BIT(b);
            
      if (IsGfnOrMersenneForm(il_K, b, +1))
      {
         il_TermCount--;
         iv_PlusTerms[bit] = false;
         removedCount++;
      }
      
      if (IsGfnOrMersenneForm(il_K, b, -1))
      {
         il_TermCount--;
         iv_MinusTerms[bit] = false;
         removedCount++;
      }
   }
   
   removedCount = EliminateTermsWithSimpleRoots();
   
   if (removedCount > 0)
      WriteToConsole(COT_OTHER, "%u terms removed as they were of GFN or Mersenne form", removedCount);
}

Worker *KBBApp::CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested)
{
   Worker *theWorker;

   // Note that KBB inherits from Worker.  This will not
   // only create the worker, but also start it.
   theWorker = new KBBWorker(id, this);

   return theWorker;
}

void KBBApp::ProcessInputTermsFile(bool haveBitMap)
{
   FILE    *fPtr = fopen(is_InputTermsFileName.c_str(), "r");
   char     buffer[200], *pos;
   uint32_t b;
   int32_t  c;
   uint64_t sieveLimit;

   if (!fPtr)
      FatalError("Unable to open input file %s", is_InputTermsFileName.c_str());

   if (fgets(buffer, 1000, fPtr) == NULL)
      FatalError("No data in input file %s", is_InputTermsFileName.c_str());
   
   if (sscanf(buffer, "ABC %" SCNu64"", &il_K) != 1)
      FatalError("Line 1 is malformed in input file %s", is_InputTermsFileName.c_str());

   pos = strstr(buffer, "//");
   if (pos)
      if (sscanf(pos+13, "%" SCNu64"", &sieveLimit) == 1)
         SetMinPrime(sieveLimit);

   if (!haveBitMap)
      ii_MinB = ii_MaxB = 0;
   
   while (fgets(buffer, 1000, fPtr) != NULL)
   {
      if (!StripCRLF(buffer))
         continue;
   
      if (sscanf(buffer, "%u %d", &b, &c) != 2)
         FatalError("Line %s is malformed", buffer);

      if (!ii_MaxB)
         ii_MinB = ii_MaxB = b;

      if (haveBitMap)
      {
         il_TermCount++;
         if (c == +1) iv_PlusTerms[BIT(b)] = true;
         if (c == -1) iv_MinusTerms[BIT(b)] = true;
      }
      else
      {
         if (ii_MinB > b) ii_MinB = b;
         if (ii_MaxB < b) ii_MaxB = b;
      }
   }

   fclose(fPtr);
}

bool KBBApp::ApplyFactor(const char *term)
{
   uint64_t k;
   uint32_t b1, b2;
   char     sign;
   uint32_t c;
   
   if (sscanf(term, "%" SCNu64"*%u^%u%c%u", &k, &b1, &b2, &sign, &c) != 5)
      FatalError("Could not parse term %s\n", term);

   if (k != il_K)
      FatalError("expecting k of %" PRIu64", but found %" PRIu64"", il_K, k);
   
   if (b1 != b2)
      FatalError("b values for term %s do not match (%u != %u)", term, b1, b2);
   
   if (sign != '+' && sign != '-')
      FatalError("expected sign of + or -, but found %c", sign);
   
   if (c != 1)
      FatalError("expecting +1 or -1 term, but found %c%d", sign, c);

   if (b1 < ii_MinB || b1 > ii_MaxB)
      return false;
   
   uint64_t bit = BIT(b1);
   
   // No locking is needed because the Workers aren't running yet
   if (sign == '+' && iv_PlusTerms[bit])
   {
      iv_PlusTerms[bit] = false;
      il_TermCount--;
      return true;
   }
   
   if (sign == '-' && iv_MinusTerms[bit])
   {
      iv_MinusTerms[bit] = false;
      il_TermCount--;
      return true;
   }

   return false;
}

void KBBApp::WriteOutputTermsFile(uint64_t largestPrime)
{
   FILE    *termsFile = fopen(is_OutputTermsFileName.c_str(), "w");
   uint32_t remaining = 0, bit;

   if (!termsFile)
      FatalError("Unable to open output file %s", is_OutputTermsFileName.c_str());
   
   ip_FactorAppLock->Lock();
   
   fprintf(termsFile, "ABC %" PRIu64"*$a^$a+$b // Sieved to %" PRIu64"\n", il_K, largestPrime);

   bit = BIT(ii_MinB);
   
   for (uint32_t b=ii_MinB; b<=ii_MaxB; b++)
   {
      if (iv_MinusTerms[bit])
      {
         fprintf(termsFile, "%u -1\n", b);
         remaining++;
      }
      
      if (iv_PlusTerms[bit])
      {
         fprintf(termsFile, "%u +1\n", b);
         remaining++;
      }
      
      bit++;
   }

   fclose(termsFile);
   
   if (remaining != il_TermCount)
      FatalError("Terms expected != terms counted (%d != %d)", il_TermCount, remaining);
        
   ip_FactorAppLock->Release();
}

void  KBBApp::GetBases(uint32_t *bases)
{
   uint32_t idx = 0;
   uint32_t bit = BIT(ii_MinB);
   
   for (uint32_t b=ii_MinB; b<=ii_MaxB; b++)
   {
      if (iv_PlusTerms[bit] || iv_MinusTerms[bit])
      {
         bases[idx] = b;
         idx++;
      }
      
      bit++;
   }
}

void KBBApp::GetExtraTextForSieveStartedMessage(char *extraTtext)
{
   sprintf(extraTtext, "k = %" PRIu64", %u <= b <= %u", il_K, ii_MinB, ii_MaxB);
}

bool  KBBApp::ReportFactor(uint64_t p, uint32_t b, int32_t c)
{
   bool     removedTerm = false;
   
   if (b < ii_MinB) return false;
   if (b > ii_MaxB) return false;
   
   if (p > GetMaxPrimeForSingleWorker())
      ip_FactorAppLock->Lock();

   uint32_t bit = BIT(b);
   
   if (c == +1 && iv_PlusTerms[bit])
   {
      iv_PlusTerms[bit] = false;
      removedTerm = true;
      
      LogFactor(p, "%" PRIu64"*%u^%u+1", il_K, b, b);
      
      il_FactorCount++;
      il_TermCount--;
   }
   
   if (c == -1 && iv_MinusTerms[bit])
   {
      iv_MinusTerms[bit] = false;
      removedTerm = true;
      
      LogFactor(p, "%" PRIu64"*%u^%u-1", il_K, b, b);
      
      il_FactorCount++;
      il_TermCount--;
   }
   
   if (p > GetMaxPrimeForSingleWorker())
      ip_FactorAppLock->Release();
   
   return removedTerm;
}

void  KBBApp::EliminateTermsWithAlgebraicFactors(void)
{
   uint32_t remvoed1 = 0;
   
   remvoed1 = EliminateTermsWithSimpleRoots();
   
   if (remvoed1 > 0)
      WriteToConsole(COT_OTHER, "%u terms removed due to algebraic factorizations", remvoed1);
}

// If k=x^f and b=y^g then look for algebraic factors.
// Note that all n will be covered by these factors.
//
// If c=-1, z divides f, z divides n*g, any n*g, then x^(f/z)*y^((n*g)/z)-1 will be a factor.
// If c=+1, z divides f, z divides n*g, odd n*g, then x^(f/z)*y^((n*g)/z)+1 will be a factor.
uint32_t  KBBApp::EliminateTermsWithSimpleRoots(void)
{
   uint32_t  removedCount = 0;
   uint32_t  b, idx;
   uint64_t  broot, kroot, curroot;
   uint32_t  bpower, kpower, curpower;
   
   // If k = 1, then RemoveAlgebraicSimpleTerm() will have handled it
   if (il_K == 1)
      return 0;
   
   GetRoot(il_K, &kroot, &kpower);

   // il_AFk must be x^f for some x and f must be greater than 1
   if (kpower == 1)
      return 0;

   curroot = kroot;
      
   for (idx=2; idx<=kpower; idx++)
   {
      // curroot = kroot^idx
      curroot *= kroot;
            
      // Given k=kroot^kpower, find all idx where kpower%idx = 0.
      // If kpower == 6, then we look for algebraic factors with the
      // forms kroot^(6/2), kroot^(6/3), and kroot^(6/6).
      if (kpower % idx != 0)
         continue;
      
      curpower = kpower / idx;
      
      for (b=ii_MinB; b<=ii_MaxB; b++)
      {
         GetRoot(b, &broot, &bpower);
         
         if ((b*bpower) % curpower != 0)
            continue;
                  
         removedCount += CheckAlgebraicFactor(b, -1, "(%" PRIu64"*%" PRIu64"^%u-1)", curroot, broot, (b*bpower)/curpower);

         // x^1 is not a divisor of x^n+1
         if (curpower == 1)
            continue;
         
         // x^y does not divide x^(y*n)+1 when y*n is even
         if ((b*bpower) % 2 == 0)
            continue;

         removedCount += CheckAlgebraicFactor(b, +1, "(%" PRIu64"*%" PRIu64"^%u+1)", curroot, broot, (b*bpower)/curpower);         
      }
   }
   
   return removedCount;
}

bool  KBBApp::CheckAlgebraicFactor(uint32_t b, int32_t c, const char *fmt, ...)
{   
   va_list args;
   char    factor[200];
   bool     removedTerm = false;
   
   if (b < ii_MinB || b > ii_MaxB)
      return false;
   
   if (c != +1 && c != -1)
      return false;
   
   va_start(args,fmt);
   vsprintf(factor, fmt, args);
   va_end(args);
      
   uint32_t bit = BIT(b);
   
   if (c == +1 && iv_PlusTerms[bit])
   {
      iv_PlusTerms[bit] = false;
      removedTerm = true;
      
      LogFactor(factor, "%" PRIu64"*%u^%u+1", il_K, b, b);
      
      il_FactorCount++;
      il_TermCount--;
   }
   
   if (c == -1 && iv_MinusTerms[bit])
   {
      iv_MinusTerms[bit] = false;
      removedTerm = true;
      
      LogFactor(factor, "%" PRIu64"*%u^%u-1", il_K, b, b);
      
      il_FactorCount++;
      il_TermCount--;
   }
   
   return removedTerm;
}

