/* SierpinskiRieselApp.cpp -- (C) Mark Rodenkirch, October 2018

   Sieve for k*b^n+c for a range of n and fixed k, n, and c

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include "../core/inline.h"
#include "../core/Parser.h"
#include "../core/Clock.h"
#include "../sieve/primesieve.hpp"
#include "SierpinskiRieselApp.h"
#include "AlgebraicFactorHelper.h"
#include "../x86_asm/fpu-asm-x86.h"

#include "GenericSequenceHelper.h"
#include "CisOneWithOneSequenceHelper.h"

#define APP_VERSION     "1.5.1"

#ifdef HAVE_GPU_WORKERS
#define APP_NAME        "srsieve2cl"
#else
#define APP_NAME        "srsieve2"
#endif

#define NBIT(n)         ((n) - ii_MinN)

// This is declared in App.h, but implemented here.  This means that App.h
// can remain unchanged if using the mtsieve framework for other applications.
App *get_app(void)
{
   return new SierpinskiRieselApp();
}

SierpinskiRieselApp::SierpinskiRieselApp() : FactorApp()
{
   SetBanner(APP_NAME " v" APP_VERSION ", a program to find factors of k*b^n+c numbers for fixed b and variable k and n");
   SetLogFileName("srsieve2.log");
   
   SetAppMinPrime(3);
   SetAppMaxPrime(PMAX_MAX_62BIT);

   ii_MinN = 0;
   ii_MaxN = 0;
   it_Format = FF_ABCD;
   ib_HaveNewSequences = false;
   il_MaxK = 0;
   il_MaxAbsC = 0;
   
   ip_FirstSequence = NULL;
   ii_SequenceCount = 0;
   ii_SquareFreeB = 0;
   is_LegendreFileName = "";
   ib_UseLengendreTables = true;
   is_SequencesToRemove = "";
   
#ifdef HAVE_GPU_WORKERS
   ib_UseGPUWorkersUponRebuild = false;
   ii_GpuFactorDensity = 10;
#endif
}

void SierpinskiRieselApp::Help(void)
{
   FactorApp::ParentHelp();

   printf("-n --nmin=n           Minimum n to search\n");
   printf("-N --nmax=N           Maximum n to search\n");
   printf("-s --sequence=s       Sequence in form k*b^n+c where k, b, and c are decimal values\n");
   printf("-f --format=f         Format of output file (A=ABC, D=ABCD (default), B=BOINC, P=ABC with number_primes)\n");
   printf("-l --legendre         Disable use of Legendre tables\n");
   printf("-L --legendrefile=L   Input/output file for Legendre tables (tables kept in memory only if -l used without -L)\n");
   printf("-R --remove=r         Remove sequences\n");
   
#ifdef HAVE_GPU_WORKERS
   printf("-M --maxfactordensity=M   factors per 1e6 terms per GPU worker chunk (default %u)\n", ii_GpuFactorDensity);
#endif
}

void  SierpinskiRieselApp::AddCommandLineOptions(string &shortOpts, struct option *longOpts)
{
   FactorApp::ParentAddCommandLineOptions(shortOpts, longOpts);

   shortOpts += "n:N:s:f:lL:R:";

   AppendLongOpt(longOpts, "nmin",           required_argument, 0, 'n');
   AppendLongOpt(longOpts, "nmax",           required_argument, 0, 'N');
   AppendLongOpt(longOpts, "sequence",       required_argument, 0, 's');
   AppendLongOpt(longOpts, "format",         required_argument, 0, 'f');
   AppendLongOpt(longOpts, "legendre",       no_argument,       0, 'l');
   AppendLongOpt(longOpts, "legendrefile",   required_argument, 0, 'L');
   AppendLongOpt(longOpts, "remove",         required_argument, 0, 'R');
   
#ifdef HAVE_GPU_WORKERS
   shortOpts += "M:";
   
   AppendLongOpt(longOpts, "maxfactordensity",  required_argument, 0, 'M');
#endif
}

parse_t SierpinskiRieselApp::ParseOption(int opt, char *arg, const char *source)
{
   parse_t status = P_UNSUPPORTED;

   status = FactorApp::ParentParseOption(opt, arg, source);
   if (status != P_UNSUPPORTED) return status;

   switch (opt)
   {
      case 'n':
         status = Parser::Parse(arg, 1, NMAX_MAX, ii_MinN);
         break;

      case 'N':
         status = Parser::Parse(arg, 1, NMAX_MAX, ii_MaxN);
         break;
         
      case 'l':
         ib_UseLengendreTables = false;
         status = P_SUCCESS;
         break;

      case 'L':
         is_LegendreFileName = arg;
         status = P_SUCCESS;
         break;
         
      case 'f':
         char value;
         status = Parser::Parse(arg, "ABDP", value);
         
         it_Format = FF_UNKNOWN;
   
         if (value == 'A')
            it_Format = FF_ABC;
         if (value == 'B')
            it_Format = FF_BOINC;
         if (value == 'D')
            it_Format = FF_ABCD;
         if (value == 'P')
            it_Format = FF_NUMBER_PRIMES;
         break;
         
      case 's':
         if (!LoadSequencesFromFile(arg))
            ValidateAndAddNewSequence(arg);

         status = P_SUCCESS;
         break;

      case 'R':
         is_SequencesToRemove = arg;
         status = P_SUCCESS;
         break;

#ifdef HAVE_GPU_WORKERS
      case 'M':
         status = Parser::Parse(arg, 10, 10000000, ii_GpuFactorDensity);
         break;
#endif
   }

   return status;
}

void SierpinskiRieselApp::ValidateOptions(void)
{
   seq_t     *seq;
   ib_CanUseCIsOneLogic = true;
   
   if (it_Format == FF_UNKNOWN)
      FatalError("the specified file format in not valid, use A (ABC), D (ABCD), P (ABC with number_primes), or B (BOINC)");
   
   if (is_InputTermsFileName.length() > 0)
   {
      if (ib_HaveNewSequences)
         FatalError("cannot add new candidate sequences in to an existing sieve");
         
      ProcessInputTermsFile(false);
      
      seq = ip_FirstSequence;
      do
      {
         seq->nTerms.resize(ii_MaxN - ii_MinN + 1);
         std::fill(seq->nTerms.begin(), seq->nTerms.end(), false);

         seq = (seq_t *) seq->next;
      } while (seq != NULL);

      ProcessInputTermsFile(true);

      RemoveSequences();
      
      MakeSubsequences(false, GetMinPrime());
   }
   else
   {
      if (!ib_HaveNewSequences)
         FatalError("must specify one or more sequences or an input file");
         
      if (ii_MinN == 0)
         FatalError("nmin must be specified");

      if (ii_MaxN == 0)
         FatalError("nmax must be specified");
      
      if (ii_MaxN <= ii_MinN)
         FatalError("nmax must be greater than nmin");

      AlgebraicFactorHelper *afh = new AlgebraicFactorHelper(this, ii_Base, ii_MinN, ii_MaxN);
      
      seq = ip_FirstSequence;
      do
      {
         seq->nTerms.resize(ii_MaxN - ii_MinN + 1);
         std::fill(seq->nTerms.begin(), seq->nTerms.end(), true);
         
         il_TermCount += (ii_MaxN - ii_MinN + 1);
         
         il_TermCount -= afh->RemoveTermsWithAlgebraicFactors(seq);
         
         seq = (seq_t *) seq->next;
      } while (seq != NULL);
      
      delete afh;
      
      MakeSubsequences(true, GetMinPrime());  
   }

   if (it_Format == FF_BOINC)
   {
      seq = ip_FirstSequence;
      do
      {
         if (abs(seq->c) != 1)
            FatalError("When using BOINC format, all sequences must have c=+1 or c=-1");
         
         if (seq->c != ip_FirstSequence->c)
            FatalError("When using BOINC format, cannot mix c=+1 and c=-1 sequences");

         seq = (seq_t *) seq->next;
      } while (seq != NULL);
   }

   if (is_OutputTermsFileName.length() == 0)
   {
      char  fileName[30];
      
      if (it_Format == FF_ABCD)
         sprintf(fileName, "b%u_n.abcd", ii_Base);
      if (it_Format == FF_ABC)
         sprintf(fileName, "b%u_n.abc", ii_Base);
      if (it_Format == FF_BOINC)
         sprintf(fileName, "b%u_n.boinc", ii_Base);
      if (it_Format == FF_NUMBER_PRIMES)
         sprintf(fileName, "b%u_n.abcnp", ii_Base);
      
      is_OutputTermsFileName = fileName;
   }

   if (!ib_UseLengendreTables && is_LegendreFileName.length() > 0)
   {
      WriteToConsole(COT_OTHER, "Ingoring -L option since Legendre tables cannot be used");
      is_LegendreFileName = "";
   }
   
   // Allow only one worker to do work when processing small primes.  This allows us to avoid 
   // locking when factors are reported, which significantly hurts performance as most terms 
   // will be removed due to small primes.
   
   // This will sieve beyond the limit, but we want to make sure that at least one prime
   // larger than this limit is passed to the worker even if the worker does not test it.
   SetMaxPrimeForSingleWorker(il_SmallPrimeSieveLimit + 1000);

#ifdef HAVE_GPU_WORKERS
   SetMinGpuPrime(1000000);
   
   double factors = (double) (ii_MaxN - ii_MinN) * (double) (ii_SequenceCount) / 1000.0;

   ii_MaxGpuFactors = GetGpuWorkGroups() * (uint64_t) (factors * (double) ii_GpuFactorDensity);
#endif

   FactorApp::ParentValidateOptions();
}

bool  SierpinskiRieselApp::LoadSequencesFromFile(char *fileName)
{
   FILE    *fPtr = fopen(fileName, "r");
   char     buffer[1000];

   if (!fPtr)
      return false;

   while (fgets(buffer, sizeof(buffer), fPtr) != NULL)
   {
      if (!StripCRLF(buffer))
         continue;
      
      ValidateAndAddNewSequence(buffer);
   }

   fclose(fPtr);
   
   if (!ib_HaveNewSequences)
      FatalError("No data in input file %s", fileName);
   
   return true;
}

void  SierpinskiRieselApp::ValidateAndAddNewSequence(char *arg)
{
   uint64_t k;
   uint32_t b;
   int64_t  c;
   uint32_t d;
   
   if (sscanf(arg, "(%" SCNu64"*%u^n%" SCNd64")/%u", &k, &b, &c, &d) != 4)
   {
      d = 1;
      
      if (sscanf(arg, "%" SCNu64"*%u^n%" SCNd64"", &k, &b, &c) != 3)         
         FatalError("sequence %s must be in form k*b^n+c or (k*b^n+c)/d where you specify values for k, b, c, and d", arg);
   }

   if (ib_HaveNewSequences && b != ii_Base)
      FatalError("only one bsae can be specified");
   
   if (k == 0)
      FatalError("k must be non-zero");
   
   if (c == 0)
      FatalError("c must be non-zero");
   
   if (!ib_HaveNewSequences)
   {
      if (b < 2)
         FatalError("base must be greater than 1");
      
      ii_Base = b;
   }

   if (d != 1 && it_Format == FF_BOINC)
      FatalError("d must be 1 for BOINC format");

   if (c > 0)
   {
      if (gcd64(ii_Base, c) != 1)
         FatalError("b and c must be relatively prime");
      
      if (gcd64(k, c) != 1)
         FatalError("k and c must be relatively prime");
   }
   else
   {
      if (gcd64(ii_Base, -c) != 1)
         FatalError("b and c must be relatively prime");
      
      if (gcd64(k, -c) != 1)
         FatalError("k and c must be relatively prime");
   }

   ib_HaveNewSequences = true;
   
   AddSequence(k, c, d);
}

Worker *SierpinskiRieselApp::CreateWorker(uint32_t id,  bool gpuWorker, uint64_t largestPrimeTested)
{
   return ip_AppHelper->CreateWorker(id, gpuWorker, largestPrimeTested);
}

void SierpinskiRieselApp::ProcessInputTermsFile(bool haveBitMap)
{
   FILE    *fPtr = fopen(is_InputTermsFileName.c_str(), "r");
   char     buffer[1000];
   uint32_t n, diff;
   uint64_t k;
   int64_t  c;
   uint32_t d = 1;
   uint64_t lastPrime = 0;
   uint32_t lineNumber = 0;
   format_t format = FF_UNKNOWN;
   seq_t   *currentSequence = 0;

   if (!fPtr)
      FatalError("Unable to open input file %s", is_InputTermsFileName.c_str());
   
   if (!haveBitMap)
      ii_MinN = ii_MaxN = 0;
   
   while (fgets(buffer, sizeof(buffer), fPtr) != NULL)
   {
      lineNumber++;
      
      if (!StripCRLF(buffer))
         continue;

      if (!memcmp(buffer, "pmin=", 5))
      {
         if (sscanf(buffer, "pmin=%" SCNu64"", &lastPrime) != 1)
            FatalError("Line %u is not a valid pmin line in input file %s", lineNumber, is_InputTermsFileName.c_str());
      }
      else if (strstr(buffer, ":P:") != NULL)
      {
         if (sscanf(buffer, "%" SCNu64":P:1:%u:257", &lastPrime, &ii_Base) != 2)
            FatalError("Line %u is not a valid BOINC line in input file %s", lineNumber, is_InputTermsFileName.c_str());
         
         format = FF_BOINC;
         c = +1;
      }
      else if (strstr(buffer, ":M:") != NULL)
      {
         if (sscanf(buffer, "%" SCNu64":M:1:%u:258", &lastPrime, &ii_Base) != 2)
            FatalError("Line %u is not a valid BOINC line in input file %s", lineNumber, is_InputTermsFileName.c_str());
         
         format = FF_BOINC;
         c = -1;
      }
      else if (!memcmp(buffer, "ABCD ", 5))
      {
         if (strstr(buffer, "Sieved") != NULL)
         {
            if (sscanf(buffer, "ABCD (%" SCNu64"*%u^$a%" SCNd64")/%u [%u] // Sieved to %" SCNu64"", &k, &ii_Base, &c, &d, &n, &lastPrime) != 6)
            {
               d = 1;
               
               if (sscanf(buffer, "ABCD %" SCNu64"*%u^$a%" SCNd64" [%u] // Sieved to %" SCNu64"", &k, &ii_Base, &c, &n, &lastPrime) != 5)
                  FatalError("Line %u is not a valid ABCD line in input file %s", lineNumber, is_InputTermsFileName.c_str());
            }
         }
         else
         {
            if (sscanf(buffer, "ABCD (%" SCNu64"*%u^$a%" SCNd64")/%u [%u]", &k, &ii_Base, &c, &d, &n) != 5)
            {
               d = 1;
               
               if (sscanf(buffer, "ABCD %" SCNu64"*%u^$a%" SCNd64" [%u]", &k, &ii_Base, &c, &n) != 4)
                  FatalError("Line %u is not a valid ABCD line in input file %s", lineNumber, is_InputTermsFileName.c_str());
            }
         }
         
         format = FF_ABCD;

         if (haveBitMap)
         {
            currentSequence = GetSequence(k, c, d);
            currentSequence->nTerms[NBIT(n)] = true;
            il_TermCount++;
            continue;
         }
         else
         {            
            if (ii_MinN == 0)
               ii_MinN = ii_MaxN = n;
            
            if (ii_MinN > n) ii_MinN = n;
            if (ii_MaxN < n) ii_MaxN = n;
            
            AddSequence(k, c, d);
         }
      }
      else if (strstr(buffer, "number_primes") != NULL)
      {
         if (strstr(buffer, "Sieved") != NULL)
         {
            if (sscanf(buffer, "ABC ($a*%u^$b$c)/$d // {number_primes,$a,1} Sieved to %" SCNu64"", &ii_Base, &lastPrime) != 2)
            {              
               if (sscanf(buffer, "ABC $a*%u^$b$c // {number_primes,$a,1} Sieved to %" SCNu64"", &ii_Base, &lastPrime) != 2)
                  FatalError("Line %u is not a valid ABC line in input file %s", lineNumber, is_InputTermsFileName.c_str());
            }
         }
         else
         {
            if (sscanf(buffer, "ABC ($a*%u^$b$c)/$d // {number_primes,$a,1}", &ii_Base) != 1)
            {
               if (sscanf(buffer, "ABC $a*%u^$b$c // {number_primes,$a,1}", &ii_Base) != 1)
                  FatalError("Line %u is not a valid ABC line in input file %s", lineNumber, is_InputTermsFileName.c_str());
            }
         }
         
         format = FF_NUMBER_PRIMES;
      }
      else if (!memcmp(buffer, "ABC ", 4))
      {
         if (strstr(buffer, "Sieved") != NULL)
         {
            if (sscanf(buffer, "ABC (%" SCNu64"*%u^$a%" SCNd64")/%u // Sieved to %" SCNu64"", &k, &ii_Base, &c, &d, &lastPrime) != 5)
            {
               d = 1;
               
               if (sscanf(buffer, "ABC %" SCNu64"*%u^$a%" SCNd64" // Sieved to %" SCNu64"", &k, &ii_Base, &c, &lastPrime) != 4)
                  FatalError("Line %u is not a valid ABC line in input file %s", lineNumber, is_InputTermsFileName.c_str());
            }
         }
         else
         {
            if (sscanf(buffer, "ABC (%" SCNu64"*%u^$a%" SCNd64")/%u", &k, &ii_Base, &c, &d) != 4)
            {
               d = 1;
               
               if (sscanf(buffer, "ABC %" SCNu64"*%u^$a%" SCNd64"", &k, &ii_Base, &c) != 3)
                  FatalError("Line %u is not a valid ABC line in input file %s", lineNumber, is_InputTermsFileName.c_str());
            }
         }
         
         format = FF_ABC;
         
         if (haveBitMap)
            currentSequence = GetSequence(k, c, d);
         else
            AddSequence(k, c, d);
      }
      else
      {
         switch (format)
         {
            case FF_ABCD:
               if (sscanf(buffer, "%u", &diff) != 1)
                  FatalError("Line %s is malformed", buffer);
               
               n += diff;
               break;

            case FF_ABC:
               if (sscanf(buffer, "%u", &n) != 1)
                  FatalError("Line %s is malformed", buffer);
               
               if (ii_MinN == 0)
                  ii_MinN = ii_MaxN = n;
               break;

            case FF_BOINC:
               if (sscanf(buffer, "%" SCNu64" %u", &k, &n) != 2)
                  FatalError("Line %s is malformed", buffer);
               
               if (haveBitMap)
                  currentSequence = GetSequence(k, c, d);
               else
                  AddSequence(k, c, d);
                  
               if (ii_MinN == 0)
                  ii_MinN = ii_MaxN = n;
               break;
               
            case FF_NUMBER_PRIMES:
               if (sscanf(buffer, "%" SCNu64" %u %" SCNd64" %u", &k, &n, &c, &d) != 4)
                  FatalError("Line %s is malformed", buffer);
               
               if (haveBitMap)
                  currentSequence = GetSequence(k, c, d);
               else
                  AddSequence(k, c, d);
               
               if (ii_MinN == 0)
                  ii_MinN = ii_MaxN = n;
               break;
               
            default:
               FatalError("Input file %s has unknown format", is_InputTermsFileName.c_str());
         }

         if (haveBitMap)
         {
            currentSequence->nTerms[NBIT(n)] = true;
            il_TermCount++;
         }
         else
         {
            if (ii_MinN > n) ii_MinN = n;
            if (ii_MaxN < n) ii_MaxN = n;
         }
      }
   }

   fclose(fPtr);

   if (lastPrime > 0)
      SetMinPrime(lastPrime);
      
   if (ii_SequenceCount == 0)
      FatalError("No sequences in input file %s", is_InputTermsFileName.c_str());
}

void SierpinskiRieselApp::RemoveSequences(void)
{
   if (is_SequencesToRemove.length() == 0)
      return;

   FILE       *fPtr = fopen(is_SequencesToRemove.c_str(), "r");
   char        buffer[1000];

   if (!fPtr)
   {
      RemoveSequence(is_SequencesToRemove.c_str());
      return;
   }

   while (fgets(buffer, sizeof(buffer), fPtr) != NULL)
   {
      if (!StripCRLF(buffer))
         continue;

      RemoveSequence(buffer);
   }

   fclose(fPtr);

}

void  SierpinskiRieselApp::RemoveSequence(const char *sequence)
{
   uint64_t    k;
   uint32_t    b, n, d;
   int64_t     c;
   
   if (sscanf(sequence, "(%" SCNu64"*%u^%u%" SCNd64")/%u", &k, &b, &n, &c, &d) == 5)
      return RemoveSequence(k, b, c, d);
      
   if (sscanf(sequence, "(%" SCNu64"*%u^n%" SCNd64")/%u", &k, &b, &c, &d) == 4)
      return RemoveSequence(k, b, c, d);
   
   if (sscanf(sequence, "%" SCNu64"*%u^%u%" SCNd64"", &k, &b, &n, &c) == 4)
      return RemoveSequence(k, b, c, 1);
      
   if (sscanf(sequence, "%" SCNu64"*%u^n%" SCNd64"", &k, &b, &c) == 3)
      return RemoveSequence(k, b, c, 1);
 
   FatalError("Sequence to remove (%s) must be in form k*b^n+c or (k*b^n+c)/d where you specify values for k, b, c, and d", sequence);
}

void  SierpinskiRieselApp::RemoveSequence(uint64_t k, uint32_t b, int64_t c, uint32_t d)
{
   seq_t   *seq, *nextSeq, *prevSeq;
   char     sequence[100];
   bool     found = false;
   uint32_t removedCount = 0;
   
   if (d == 1)
      sprintf(sequence, "%" PRIu64"*%u^n%+" PRId64"", k, b, c);
   else
      sprintf(sequence, "(%" PRIu64"*%u^n%+" PRId64")/%u", k, b, c, d);
   
   if (b != ii_Base)
      WriteToConsole(COT_OTHER, "Sequence %s wasn't removed because it is not base %u", sequence, b);
   
   prevSeq = seq = ip_FirstSequence;
   do
   {      
      nextSeq = (seq_t *) seq->next;
      
      if (seq->k == k && seq->c == c && seq->d == d)
      {
         found = true;
         
         for (uint32_t n=ii_MinN; n<=ii_MaxN; n++)
         {
            if (seq->nTerms[NBIT(n)])
               removedCount++;
         }

         seq->nTerms.clear();
                     
         if (seq == ip_FirstSequence)
            ip_FirstSequence = (seq_t *) seq->next;
         else
            prevSeq->next = seq->next;
            
         WriteToConsole(COT_OTHER, "%u terms for sequence %s have been removed", removedCount, sequence);
                  
         xfree(seq);
         
         il_TermCount -= removedCount;
         ii_SequenceCount--;
      }
      
      prevSeq = seq;
      seq = nextSeq;
   } while (seq != NULL);

   
   if (!found)
      WriteToConsole(COT_OTHER, "Sequence %s wasn't removed because it wasn't found", sequence);
}

bool SierpinskiRieselApp::ApplyFactor(uint64_t thePrime, const char *term)
{
   uint64_t   k;
   uint32_t   b, n, d;
   int64_t    c;
   seq_t     *seq;
   
   if (sscanf(term, "(%" SCNu64"*%u^%u%" SCNd64")/%u", &k, &b, &n, &c, &d) != 5)
   {
      d = 1;
      
      if (sscanf(term, "%" SCNu64"*%u^%u%" SCNd64"", &k, &b, &n, &c) != 4)
         FatalError("Could not parse term %s", term);
   }

   if (b != ii_Base)
      FatalError("Expected base %u in factor but found base %u", ii_Base, b);
        
   if (n < ii_MinN || n > ii_MaxN)
      return false;
   
   seq = ip_FirstSequence;
   do
   {
      if (seq->k == k && seq->c == c && seq->d == d)
      {
         if (seq->nTerms[NBIT(n)])
         {
            seq->nTerms[NBIT(n)] = false;
            il_TermCount--;

            if (!VerifyFactor(false, thePrime, seq, n))
               return false;
         
            return true;
         }
      }

      seq = (seq_t *) seq->next;
   } while (seq != NULL);

   return false;
}

void SierpinskiRieselApp::WriteOutputTermsFile(uint64_t largestPrime)
{
   uint32_t   nCount = 0; 
   bool       allSequencesHaveDEqual1 = true;
   seq_t     *seq;
   
   // With super large ranges, wait until we can lock because without locking
   // the term count can change between opening and closing the file.
   if (IsRunning() && largestPrime < GetMaxPrimeForSingleWorker())
      return;
   
   FILE    *termsFile = fopen(is_OutputTermsFileName.c_str(), "w");

   if (!termsFile)
      FatalError("Unable to open output file %s", is_OutputTermsFileName.c_str());
   
   ip_FactorAppLock->Lock();
   
   nCount = 0;

   if (it_Format == FF_NUMBER_PRIMES)
   {
      seq = ip_FirstSequence;
      do
      {
         if (seq->d != 1)
            allSequencesHaveDEqual1 = false;
            
         seq = (seq_t *) seq->next;
      } while (seq != NULL);
         
      if (allSequencesHaveDEqual1)
         fprintf(termsFile, "ABC $a*%u^$b$c // {number_primes,$a,1} Sieved to %" PRIu64"\n", ii_Base, largestPrime);
      else
         fprintf(termsFile, "ABC ($a*%u^$b$c)/$d // {number_primes,$a,1} Sieved to %" PRIu64"\n", ii_Base, largestPrime);
   }

   if (it_Format == FF_BOINC)
   {
      if (ip_FirstSequence->c == +1)
         fprintf(termsFile, "%" PRIu64":P:1:%u:257\n", largestPrime, ii_Base);
      else
         fprintf(termsFile, "%" PRIu64":M:1:%u:258\n", largestPrime, ii_Base);
   }
   
   seq = ip_FirstSequence;
   do
   {
      if (it_Format == FF_ABCD)
         nCount += WriteABCDTermsFile(seq, largestPrime, termsFile);
      
      if (it_Format == FF_ABC)
         nCount += WriteABCTermsFile(seq, largestPrime, termsFile);
      
      if (it_Format == FF_BOINC)
         nCount += WriteBoincTermsFile(seq, largestPrime, termsFile);
      
      if (it_Format == FF_NUMBER_PRIMES)
         nCount += WriteABCNumberPrimesTermsFile(seq, largestPrime, termsFile, allSequencesHaveDEqual1);
            
      seq = (seq_t *) seq->next;
   } while (seq != NULL);
   
   fclose(termsFile);
   
   if (nCount != il_TermCount)
      FatalError("Something is wrong.  Counted terms (%u) != expected terms (%u)", nCount, il_TermCount);
   
   ip_FactorAppLock->Release();
}

uint32_t SierpinskiRieselApp::WriteABCDTermsFile(seq_t *seq, uint64_t maxPrime, FILE *termsFile)
{
   uint32_t n, nCount = 0, previousN;
   uint32_t bit;

   n = ii_MinN;
   
   bit = NBIT(n);
   for (; n<=ii_MaxN; n++)
   {      
      if (seq->nTerms[bit])
         break;
      
      bit++;
   }

   if (n > ii_MaxN)
      return 0;
   
   if (seq->d == 1)
      fprintf(termsFile, "ABCD %" PRIu64"*%u^$a%+" PRId64" [%u] // Sieved to %" PRIu64"\n", seq->k, ii_Base, seq->c, n, maxPrime);
   else
      fprintf(termsFile, "ABCD (%" PRIu64"*%u^$a%+" PRId64")/%u [%u] // Sieved to %" PRIu64"\n", seq->k, ii_Base, seq->c, seq->d, n, maxPrime);
   
   previousN = n;
   nCount = 1;
   n++;
   
   bit = NBIT(n);
   for (; n<=ii_MaxN; n++)
   {
      if (seq->nTerms[bit])
      {
         fprintf(termsFile, "%u\n", n - previousN);
         previousN = n;
         nCount++;
      }
      
      bit++;
   }

   return nCount;
}

uint32_t SierpinskiRieselApp::WriteABCTermsFile(seq_t *seq, uint64_t maxPrime, FILE *termsFile)
{
   uint32_t n, nCount = 0;
   uint32_t bit;

   if (seq->d == 1)
      fprintf(termsFile, "ABC %" PRIu64"*%u^$a%+" PRId64" // Sieved to %" PRIu64"\n", seq->k, ii_Base, seq->c, maxPrime);
   else
      fprintf(termsFile, "ABC (%" PRIu64"*%u^$a%+" PRId64")/%u // Sieved to %" PRIu64"\n", seq->k, ii_Base, seq->c, seq->d, maxPrime);
      
   n = ii_MinN;
   bit = NBIT(n);

   for (; n<=ii_MaxN; n++)
   {
      if (seq->nTerms[bit])
      {
         fprintf(termsFile, "%u\n", n);
         nCount++;
      }
      
      bit++;
   }
   
   return nCount;
}

uint32_t SierpinskiRieselApp::WriteBoincTermsFile(seq_t *seq, uint64_t maxPrime, FILE *termsFile)
{
   uint32_t n, nCount = 0;
   uint32_t bit;
      
   n = ii_MinN;
   bit = NBIT(n);

   for (; n<=ii_MaxN; n++)
   {
      if (seq->nTerms[bit])
      {
         fprintf(termsFile, "%" PRIu64" %u\n", seq->k, n);
         nCount++;
      }
      
      bit++;
   }
   
   return nCount;
}

uint32_t SierpinskiRieselApp::WriteABCNumberPrimesTermsFile(seq_t *seq, uint64_t maxPrime, FILE *termsFile, bool allSequencesHaveDEqual1)
{
   uint32_t n, nCount = 0;
   uint32_t bit;

   n = ii_MinN;
   bit = NBIT(n);

   for (; n<=ii_MaxN; n++)
   {
      if (seq->nTerms[bit])
      {
         if (allSequencesHaveDEqual1)
            fprintf(termsFile, "%" PRIu64" %u %+" PRId64"\n", seq->k, n, seq->c);
         else
            fprintf(termsFile, "%" PRIu64" %u %+" PRId64" %u\n", seq->k, n, seq->c, seq->d);
         
         nCount++;
      }
      
      bit++;
   }
   
   return nCount;
}

void  SierpinskiRieselApp::GetExtraTextForSieveStartedMessage(char *extraTtext)
{
   sprintf(extraTtext, "%u < n < %u, k*%u^n+c", ii_MinN, ii_MaxN, ii_Base);
}

void  SierpinskiRieselApp::AddSequence(uint64_t k, int64_t c, uint32_t d)
{
   seq_t  *seq;

   // If base, k, and c are odd then all terms are even
   if ((ii_Base % 2) && (k % 2) && (c % 2) && (d == 1))
   {
      WriteToConsole(COT_OTHER, "Sequence %" PRIu64"*%u^n%+d not added because all terms are divisible by 2", k, ii_Base, c);
      return;
   }

   // If either the base or k is even and c is even then all terms are even
   if ((!(ii_Base % 2) || !(k % 2)) && !(c % 2) && (d == 1))
   {
      WriteToConsole(COT_OTHER, "Sequence %" PRIu64"*%u^n%+d not added because all terms are divisible by 2", k, ii_Base, c);
      return;
   }
   
   if (ip_FirstSequence != NULL)
   {
      // If the sequence already exists, then nothing to do.
      seq = ip_FirstSequence;
      do
      {
         if (seq->k == k && seq->c == c && seq->d == d)
            return;

         seq = (seq_t *) seq->next;
      } while (seq != NULL);
   }
       
   seq_t *newPtr = (seq_t *) xmalloc(sizeof(seq_t));
   
   uint64_t absc = abs(c);

   if (absc != 1)
      ib_CanUseCIsOneLogic = false;
   
   if (il_MaxK < k) il_MaxK = k;
   if (il_MaxAbsC < absc) il_MaxAbsC = absc;

   newPtr->k = k;
   newPtr->c = c;
   newPtr->d = d;
   newPtr->next = NULL;
   
   if (ii_SequenceCount == 0)
      ip_FirstSequence = newPtr;
   else
   {
      seq = ip_FirstSequence;
      do
      {
         if (seq->next == NULL)
         {
            seq->next = newPtr;
            break;
         }
   
         seq = (seq_t *) seq->next;
      } while (seq != NULL);
   }
   
   ii_SequenceCount++;
}

seq_t    *SierpinskiRieselApp::GetSequence(uint64_t k, int64_t c, uint32_t d) 
{
   seq_t     *seq;

   seq = ip_FirstSequence;
   do
   {
      if (seq->k == k && seq->c == c && seq->d == d)
         return seq;

      seq = (seq_t *) seq->next;
   } while (seq != NULL);
   
   FatalError("Sequence for %" PRIu64" and c=%" PRId64" was not found", k, c);
   
   return 0;
}

// Note that this is only called if all workers are paused and waiting for work
void  SierpinskiRieselApp::NotifyAppToRebuild(uint64_t largestPrimeTested)
{
#ifdef HAVE_GPU_WORKERS
   if (ib_UseGPUWorkersUponRebuild)
      return;
#endif

   ip_AppHelper->CleanUp();
   
   delete ip_AppHelper;
   
   // This allows us to choose the best AbstractSequenceHelper based upon the current status
   MakeSubsequences(false, largestPrimeTested);
}

void  SierpinskiRieselApp::MakeSubsequences(bool newSieve, uint64_t largestPrimeTested)
{
   RemoveSequencesWithNoTerms();

   if (ii_SequenceCount == 0)
      FatalError("All sequences have been removed");
      
   // 65536 is from srsieve.  I don't understand the limit, but if this
   // value is too large, then factors are missed.
   il_SmallPrimeSieveLimit = MAX(257, ii_Base);
   
   if (il_MaxK < 65536 && il_MaxAbsC < 65536)
   {
      il_SmallPrimeSieveLimit = MAX(il_SmallPrimeSieveLimit, il_MaxK);
      il_SmallPrimeSieveLimit = MAX(il_SmallPrimeSieveLimit, il_MaxAbsC);
   }
   
   CheckForLegendreSupport();
   
   // The constructors will make a copy of the sequences.
   if (newSieve || il_MaxK > largestPrimeTested || !ib_CanUseCIsOneLogic)
      ip_AppHelper = new GenericSequenceHelper(this, largestPrimeTested);
   else
   {
      if (ii_SequenceCount == 1)
         ip_AppHelper = new CisOneWithOneSequenceHelper(this, largestPrimeTested);
      else
      {
         WriteToConsole(COT_OTHER, "Must use generic sieving logic because there is more than one sequence");
         
         ip_AppHelper = new GenericSequenceHelper(this, largestPrimeTested);
      }
   }

   if (newSieve)
   {
      uint64_t termCount = ip_AppHelper->MakeSubsequencesForNewSieve();
      
      if (termCount != il_TermCount)
         FatalError("Expected %" PRIu64" terms, but only set %" PRIu64"", il_TermCount, termCount);
   }
   else
      ip_AppHelper->MakeSubsequencesForOldSieve(il_TermCount);
   
   ip_AppHelper->LastChanceLogicBeforeSieving();
}

void  SierpinskiRieselApp::RemoveSequencesWithNoTerms(void)
{
   seq_t     *seq, *nextSeq, *prevSeq;
   
   prevSeq = seq = ip_FirstSequence;
   do
   {
      bool haveTerm = false;
      
      for (uint32_t n=ii_MinN; n<=ii_MaxN; n++)
      {
         if (seq->nTerms[NBIT(n)])
         {
            haveTerm = true;
            break;
         }
      }
      
      if (!haveTerm)
      {
         seq->nTerms.clear();
         
         nextSeq = (seq_t *) seq->next;
         
         if (seq == ip_FirstSequence)
            ip_FirstSequence = (seq_t *) seq->next;
         else
            prevSeq->next = seq->next;
            
         if (seq->d == 1)
            WriteToConsole(COT_OTHER, "Sequence %" PRIu64"*%u^n%+" PRId64" removed as all terms have a factor", 
               seq->k, ii_Base, seq->c);
         else
            WriteToConsole(COT_OTHER, "Sequence (%" PRIu64"*%u^n%+" PRId64")/%u removed as all terms have a factor", 
               seq->k, ii_Base, seq->c, seq->d);
                  
         xfree(seq);
         
         seq = nextSeq;
         ii_SequenceCount--;
      }
      else
      {
         prevSeq = seq;
         seq = (seq_t *) seq->next;
      }
   } while (seq != NULL);
}

void  SierpinskiRieselApp::CheckForLegendreSupport(void)
{
   seq_t  *seq;

   // If we have computed this value, then we don't need to do this again.
   if (ii_SquareFreeB > 0)
      return;

   if (!ib_CanUseCIsOneLogic)
   {
      WriteToConsole(COT_OTHER, "Must use generic sieving logic because abs(c) != 1 for at least one sequence");
         
      if (ib_UseLengendreTables)
      {
         WriteToConsole(COT_OTHER, "Cannot use Legendre tables because abs(c) != 1 for at least one sequence");
         ib_UseLengendreTables = false;
      }
      
      ii_SquareFreeB = 1;
      return;
   }

   // If we can use CisOne logic, then we need to compute squareFreeK and squareFreeB even if
   // Legendre tables cannot be used.
   vector<uint64_t> primes;
   
   primesieve::generate_n_primes(sqrt(il_MaxK) + 1, &primes);

   // Note that GetSquareFreeFactor() must return a value >= 1
   ii_SquareFreeB = GetSquareFreeFactor(ii_Base, primes);
   
   seq = ip_FirstSequence;
   do
   {
      seq->squareFreeK = GetSquareFreeFactor(seq->k, primes);
               
      seq->kcCore = -1 * seq->c * seq->squareFreeK;

      seq = (seq_t *) seq->next;
   } while (seq != NULL);
   
   seq = ip_FirstSequence;
   do
   {
      // This can only happen if k is prime and k > INT64_MAX
      if (seq->squareFreeK > INT64_MAX)
      {
         WriteToConsole(COT_OTHER, "Must use generic sieving logic because square-free part of k is too large");
         ib_CanUseCIsOneLogic = false;
         ib_UseLengendreTables = false;
         break;
      }
      
      if (ib_UseLengendreTables && seq->squareFreeK > INT32_MAX/ii_SquareFreeB)
      {
         WriteToConsole(COT_OTHER, "Cannot use Legendre tables because square-free part of k is too large");
         ib_UseLengendreTables = false;
         break;
      }

      seq = (seq_t *) seq->next;
   } while (seq != NULL);
   
   primes.clear();
}

void     SierpinskiRieselApp::ReportFactor(uint64_t thePrime, seq_t *seq, uint32_t n, bool verifyFactor)
{
   uint32_t nbit;
   bool     wasRemoved = false;
               
   if (n < ii_MinN || n > ii_MaxN)
      return;

   nbit = NBIT(n);
   
   if (thePrime > GetMaxPrimeForSingleWorker())
      ip_FactorAppLock->Lock();
      
   if (seq->nTerms[nbit])
   {
      il_TermCount--;
      il_FactorCount++;
      seq->nTerms[nbit] = false;
      wasRemoved = true;
   }
         
   if (thePrime > GetMaxPrimeForSingleWorker())
      ip_FactorAppLock->Release();

   if (verifyFactor)
      VerifyFactor(true, thePrime, seq, n);

   if (!wasRemoved)
      return;   
 
   char buffer[200];
   
   if (seq->d > 1)
      sprintf(buffer, "(%" PRIu64"*%u^%u%+" PRId64")/%u", seq->k, ii_Base, n, seq->c, seq->d);
   else
      sprintf(buffer, "%" PRIu64"*%u^%u%+" PRId64"", seq->k, ii_Base, n, seq->c);

   if (IsPrime(thePrime, seq, n))
   {
      WriteToConsole(COT_OTHER, "%s is prime!", buffer);
      WriteToLog("%s is prime!", buffer);
   }
   else
      LogFactor(thePrime, buffer);
}

bool  SierpinskiRieselApp::VerifyFactor(bool badFactorIsFatal, uint64_t thePrime, seq_t *seq, uint32_t n)
{
   uint64_t  rem;
   bool      isValid;

   fpu_push_1divp(thePrime);
   
   rem = fpu_powmod(ii_Base, n, thePrime);
   rem = fpu_mulmod(rem, seq->k, thePrime);

   fpu_pop();
   
   if (seq->c > 0)
      rem += seq->c;
   else
      rem += (thePrime + seq->c);
      
   if (rem >= thePrime)
      rem -= thePrime;

   // At some point need logic if gcd(d, thePrime) != 1
   isValid = (rem == 0);
   
   if (isValid)
      return isValid;
      
   char buffer[200];
   
   if (seq->d > 1)
      sprintf(buffer, "Invalid factor: (%" PRIu64"*%u^%u%+" PRId64")/%u mod %" PRIu64" = %" PRIu64"", seq->k, ii_Base, n, seq->c, seq->d, thePrime, rem);
   else
      sprintf(buffer, "Invalid factor: %" PRIu64"*%u^%u%+" PRId64" mod %" PRIu64" = %" PRIu64"", seq->k, ii_Base, n, seq->c, thePrime, rem);
   
   if (badFactorIsFatal)
      FatalError(buffer);
   else
      WriteToConsole(COT_OTHER, buffer);
   
   return isValid;
}

bool  SierpinskiRieselApp::IsPrime(uint64_t p, seq_t *seq, uint32_t n)
{
   __uint128_t bigP;
   
   // Won't try if n is too large
   if (n > 64)
      return false;

   // Just to make sure that our term isn't too big for this logic
   if (log10(ii_Base)*n > 18.0)
      return false;

   bigP = __uint128_t(p);
      
   bigP *= seq->d;
   
   // p = k*b^n+c --> p = k*b^n
   bigP -= seq->c;

   if (bigP % seq->k > 0)
      return false;
   
   // p = k*b^n --> p = b^n
   bigP /= seq->k;

   // At this if bigP = b^n, then (k*b^n+c)/d is prime
   do
   {
      // If b doesn't divide bigP, then the number is not prime
      if (bigP % ii_Base != 0)
         return false;
         
      bigP /= ii_Base;
      n--;
   } while (n > 0);
   
   // If we divided by b^n and bigP == 1, then the number is prime
   return (bigP == 1);
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
uint64_t    SierpinskiRieselApp::GetSquareFreeFactor(uint64_t n, vector<uint64_t> primes)
{
   uint64_t c = 1, q, r;
   
   vector<uint64_t>::iterator it = primes.begin();
   
   r = sqrt(n);
   
   while (it != primes.end())
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
