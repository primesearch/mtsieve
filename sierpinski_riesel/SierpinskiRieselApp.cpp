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
#include "SierpinskiRieselApp.h"
#include "AlgebraicFactorHelper.h"
#include "../x86_asm/fpu-asm-x86.h"

#include "GenericSubsequenceHelper.h"
#include "CisOneSubsequenceHelper.h"

#define APP_NAME        "srsieve2"
#define APP_VERSION     "1.2.1"

#define NBIT(n)         ((n) - ii_MinN)
#define MBIT(m)         ((m) - ii_MinM)

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
   
   SetAppMinPrime(2);
   SetAppMaxPrime(PMAX_MAX_52BIT);

   ii_MinN = 0;
   ii_MaxN = 0;
   it_Format = FF_ABCD;
   ib_HaveNewSequences = false;
   ib_UseAvx = false;
   
   ip_Sequences = 0;
   ii_SequenceCount = 0;
   ii_SequenceCapacity = 0;
   ib_UseLengendreTables = false;
   is_LegendreFileName = "";
}

void SierpinskiRieselApp::Help(void)
{
   FactorApp::ParentHelp();

   printf("-n --nmin=n           Minimum n to search\n");
   printf("-N --nmax=N           Maximum n to search\n");
   printf("-s --sequence=s       Sequence to find factors of in form k*b^n+c where b, n, and c are decimal values\n");
   printf("-f --format=f         Format of output file (A=ABC, D=ABCD (default), B=BOINC, P=ABC with number_primes)\n");
   printf("-l --legendre         Use Legendre tables\n");
   printf("-L --legendrefile=L   Input/output file for Legendre tables (tables kept in memory only if -l used without -L)\n");
//   printf("-D --disableavx       disableavx\n");
}

void  SierpinskiRieselApp::AddCommandLineOptions(string &shortOpts, struct option *longOpts)
{
   FactorApp::ParentAddCommandLineOptions(shortOpts, longOpts);

   shortOpts += "n:N:s:f:lL:";

   AppendLongOpt(longOpts, "nmin",           required_argument, 0, 'n');
   AppendLongOpt(longOpts, "nmax",           required_argument, 0, 'N');
   AppendLongOpt(longOpts, "sequence",       required_argument, 0, 's');
   AppendLongOpt(longOpts, "format",         required_argument, 0, 'f');
   AppendLongOpt(longOpts, "legendre",       no_argument, 0, 'l');
   AppendLongOpt(longOpts, "legendrefile",   required_argument, 0, 'L');
//   AppendLongOpt(longOpts, "disableavx",     no_argument, 0, 'D');
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
         ib_UseLengendreTables = true;
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

      case 'D':
         ib_UseAvx = false;
         status = P_SUCCESS;
         break;         
   }

   return status;
}

void SierpinskiRieselApp::ValidateOptions(void)
{
   if (it_Format == FF_UNKNOWN)
      FatalError("the specified file format in not valid, use A (ABC), D (ABCD), or P (ABC with number_primes))");
   
   if (is_InputTermsFileName.length() > 0)
   {
      if (ib_HaveNewSequences)
         FatalError("cannot add new candidate sequences in to an existing sieve");
         
      ProcessInputTermsFile(false);
      
      for (uint32_t seqIdx=0; seqIdx<ii_SequenceCount; seqIdx++)
      {
         ip_Sequences[seqIdx].nTerms.resize(ii_MaxN - ii_MinN + 1);
         std::fill(ip_Sequences[seqIdx].nTerms.begin(), ip_Sequences[seqIdx].nTerms.end(), false);
      }

      ProcessInputTermsFile(true);
      
      MakeSubsequences(false);
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
      
      for (uint32_t seqIdx=0; seqIdx<ii_SequenceCount; seqIdx++)
      {
         ip_Sequences[seqIdx].nTerms.resize(ii_MaxN - ii_MinN + 1);
         std::fill(ip_Sequences[seqIdx].nTerms.begin(), ip_Sequences[seqIdx].nTerms.end(), true);
         
         il_TermCount += (ii_MaxN - ii_MinN + 1);
         
         il_TermCount -= afh->RemoveTermsWithAlgebraicFactors(&ip_Sequences[seqIdx]);
      }
      
      delete afh;
      
      MakeSubsequences(true);  
   }

   if (it_Format == FF_BOINC)
   {
      for (uint32_t seqIdx=0; seqIdx<ii_SequenceCount; seqIdx++)
      {
         if (abs(ip_Sequences[seqIdx].c) != 1)
            FatalError("When using BOINC format, all sequences must have c=+1 or c=-1");
         
         if (ip_Sequences[0].c != ip_Sequences[seqIdx].c)
            FatalError("When using BOINC format, cannot mix c=+1 and c=-1 sequences");
      }
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
      WriteToConsole(COT_OTHER, "Ingoring -L option since -l was not specified");
      is_LegendreFileName = "";
   }
   
   FactorApp::ParentValidateOptions();
   
   // Allow only one worker to do work when processing small primes.  This allows us to avoid 
   // locking when factors are reported, which significantly hurts performance as most terms 
   // will be removed due to small primes.
   SetMaxPrimeForSingleWorker(il_SmallPrimeSieveLimit);
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
   
   if (sscanf(arg, "%" SCNu64"*%u^n%" SCNd64"", &k, &b, &c) != 3)         
      FatalError("sequence must be in form k*b^n+c where you specify values for k, b and c");

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
   
   AddSequence(k, c);
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

      if (strstr(buffer, ":P:") != NULL)
      {
         if (sscanf(buffer, "%" SCNu64":P:1:%u:257", &lastPrime, &ii_Base) != 2)
            FatalError("Line %u is not a valid BOINC line in intput file %s", lineNumber, is_InputTermsFileName.c_str());
         
         format = FF_BOINC;
         c = +1;
      }
      else if (strstr(buffer, ":M:") != NULL)
      {
         if (sscanf(buffer, "%" SCNu64":M:1:%u:258", &lastPrime, &ii_Base) != 2)
            FatalError("Line %u is not a valid BOINC line in intput file %s", lineNumber, is_InputTermsFileName.c_str());
         
         format = FF_BOINC;
         c = -1;
      }
      else if (!memcmp(buffer, "ABCD ", 5))
      {
         if (strstr(buffer, "Sieved") != NULL)
         {
            if (sscanf(buffer, "ABCD %" SCNu64"*%u^$a%" SCNd64" [%u] // Sieved to %" SCNu64"", &k, &ii_Base, &c, &n, &lastPrime) != 5)
               FatalError("Line %u is not a valid ABCD line in input file %s", lineNumber, is_InputTermsFileName.c_str());
         }
         else
         {
            if (sscanf(buffer, "ABCD %" SCNu64"*%u^$a%" SCNd64" [%u]", &k, &ii_Base, &c, &n) != 4)
               FatalError("Line %u is not a valid ABCD line in input file %s", lineNumber, is_InputTermsFileName.c_str());
         }
         
         format = FF_ABCD;

         if (haveBitMap)
         {
            currentSequence = GetSequence(k, c);
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
            
            AddSequence(k, c);
         }
      }
      else if (strstr(buffer, "number_primes") != NULL)
      {
         if (strstr(buffer, "Sieved") != NULL)
         {
            if (sscanf(buffer, "ABC $a*%u^$b$c // {number_primes,$a,1} Sieved to %" SCNu64"", &ii_Base, &lastPrime) != 2)
               FatalError("Line %u is not a valid ABC line in input file %s", lineNumber, is_InputTermsFileName.c_str());
         }
         else
         {
            if (sscanf(buffer, "ABC $a*%u^$b$c // {number_primes,$a,1}", &ii_Base) != 1)
               FatalError("Line %u is not a valid ABC line in input file %s", lineNumber, is_InputTermsFileName.c_str());
         }
         
         format = FF_NUMBER_PRIMES;
      }
      else if (!memcmp(buffer, "ABC ", 4))
      {
         if (strstr(buffer, "Sieved") != NULL)
         {
            if (sscanf(buffer, "ABC %" SCNu64"*%u^$a%" SCNd64" // Sieved to %" SCNu64"", &k, &ii_Base, &c, &lastPrime) != 4)
               FatalError("Line %u is not a valid ABC line in input file %s", lineNumber, is_InputTermsFileName.c_str());
         }
         else
         {
            if (sscanf(buffer, "ABC %" SCNu64"*%u^$a%" SCNd64"", &k, &ii_Base, &c) != 3)
               FatalError("Line %u is not a valid ABC line in input file %s", lineNumber, is_InputTermsFileName.c_str());
         }
         
         format = FF_ABC;
         
         if (haveBitMap)
            currentSequence = GetSequence(k, c);
         else
            AddSequence(k, c);
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
                  currentSequence = GetSequence(k, c);
               else
                  AddSequence(k, c);
               break;
               
            case FF_NUMBER_PRIMES:
               if (sscanf(buffer, "%" SCNu64" %u %" SCNd64"", &k, &n, &c) != 3)
                  FatalError("Line %s is malformed", buffer);
               
               if (haveBitMap)
                  currentSequence = GetSequence(k, c);
               else
                  AddSequence(k, c);
               
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

bool SierpinskiRieselApp::ApplyFactor(const char *term)
{
   uint64_t k;
   uint32_t b, n;
   int64_t  c;
   seq_t   *seq;
   
   if (sscanf(term, "%" SCNu64"*%u^%u%" SCNd64"", &k, &b, &n, &c) != 4)
      FatalError("Could not parse term %s", term);

   if (b != ii_Base)
      FatalError("Expected base %u in factor but found base %u", ii_Base, b);
        
   if (n < ii_MinN || n > ii_MaxN)
      return false;
   
   // No locking is needed because the Workers aren't running yet
   seq = GetSequence(k, c);
   
   if (seq->nTerms[NBIT(n)])
   {
      seq->nTerms[NBIT(n)] = false;
      il_TermCount--;

      return true;
   }
      
   return false;
}

void SierpinskiRieselApp::WriteOutputTermsFile(uint64_t largestPrime)
{
   uint32_t nCount = 0;
   uint32_t seqIdx;
   
   FILE    *termsFile = fopen(is_OutputTermsFileName.c_str(), "w");

   if (!termsFile)
      FatalError("Unable to open output file %s", is_OutputTermsFileName.c_str());
   
   ip_FactorAppLock->Lock();
   
   nCount = 0;
   
   if (it_Format == FF_NUMBER_PRIMES)
      fprintf(termsFile, "ABC $a*%u^$b$c // {number_primes,$a,1} Sieved to %" PRIu64"\n", ii_Base, largestPrime);
   
   if (it_Format == FF_BOINC)
   {
      if (ip_Sequences[0].c == +1)
         fprintf(termsFile, "%" PRIu64":P:1:%u:257\n", largestPrime, ii_Base);
      else
         fprintf(termsFile, "%" PRIu64":M:1:%u:258\n", largestPrime, ii_Base);
   }
   
   for (seqIdx=0; seqIdx<ii_SequenceCount; seqIdx++)
   {
      if (it_Format == FF_ABCD)
         nCount += WriteABCDTermsFile(&ip_Sequences[seqIdx], largestPrime, termsFile);
      
      if (it_Format == FF_ABC)
         nCount += WriteABCTermsFile(&ip_Sequences[seqIdx], largestPrime, termsFile);
      
      if (it_Format == FF_BOINC)
         nCount += WriteBoincTermsFile(&ip_Sequences[seqIdx], largestPrime, termsFile);
      
      if (it_Format == FF_NUMBER_PRIMES)
         nCount += WriteABCNumberPrimesTermsFile(&ip_Sequences[seqIdx], largestPrime, termsFile);
   }
   
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
   
   fprintf(termsFile, "ABCD %" PRIu64"*%u^$a%+" PRId64" [%u] // Sieved to %" PRIu64"\n", seq->k, ii_Base, seq->c, n, maxPrime);
   
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

   fprintf(termsFile, "ABC %" PRIu64"*%u^$a%+" PRId64" // Sieved to %" PRIu64"\n", seq->k, ii_Base, seq->c, maxPrime);
      
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

uint32_t SierpinskiRieselApp::WriteABCNumberPrimesTermsFile(seq_t *seq, uint64_t maxPrime, FILE *termsFile)
{
   uint32_t n, nCount = 0;
   uint32_t bit;

      
   n = ii_MinN;
   bit = NBIT(n);

   for (; n<=ii_MaxN; n++)
   {
      if (seq->nTerms[bit])
      {
         fprintf(termsFile, "%" PRIu64" %u %+" PRId64"\n", seq->k, n, seq->c);
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

void  SierpinskiRieselApp::AddSequence(uint64_t k, int64_t c)
{
   uint32_t seqIdx;

   // If base, k, and c are odd then all terms are even
   if ((ii_Base % 2) && (k % 2) && (c % 2))
   {
      WriteToConsole(COT_OTHER, "All terms for sequence %" PRIu64"*%u^n%+d are divisible by 2", k, ii_Base, c);
      return;
   }

   // If either the base or k is even and c is even then all terms are even
   if ((!(ii_Base % 2) || !(k % 2)) && !(c % 2))
   {
      WriteToConsole(COT_OTHER, "All terms for sequence %" PRIu64"*%u^n%+d are divisible by 2", k, ii_Base, c);
      return;
   }
   
   // If the sequence already exists, then nothing to do.
   for (seqIdx=0; seqIdx<ii_SequenceCount; seqIdx++)
   {
      if (ip_Sequences[seqIdx].k == k && ip_Sequences[seqIdx].c == c)
         return;
   }
   
   if (ii_SequenceCapacity == ii_SequenceCount)
   {
      ii_SequenceCapacity += 100;
      
      seq_t *newPtr = (seq_t *) xmalloc(ii_SequenceCapacity * sizeof(seq_t));
      
      if (ip_Sequences)
      {
         // We only need to copy these fields as the other fields will be
         // populated after we know all of the sequences.
         for (seqIdx=0; seqIdx<ii_SequenceCount; seqIdx++)
         {
            newPtr[seqIdx].k = ip_Sequences[seqIdx].k;
            newPtr[seqIdx].c = ip_Sequences[seqIdx].c;
         }
         
         xfree(ip_Sequences);
      }
      
      ip_Sequences = newPtr;
   }

   seq_t *seq = &ip_Sequences[ii_SequenceCount];
   
   seq->k = k;
   seq->c = c;
   
   ii_SequenceCount++;
}

seq_t    *SierpinskiRieselApp::GetSequence(uint64_t k, int64_t c) 
{
   uint32_t seqIdx;
   
   for (seqIdx=0; seqIdx<ii_SequenceCount; seqIdx++)
   {      
      if (ip_Sequences[seqIdx].k == k && ip_Sequences[seqIdx].c == c)
         return &ip_Sequences[seqIdx];
   }
   
   FatalError("Sequence for %" PRIu64" and c=%" PRId64" was not found", k, c);
   
   return 0;
}

// Note that this is only called if all workers are paused and waiting for work
void  SierpinskiRieselApp::NotifyAppToRebuild(void)
{
   ip_AppHelper->CleanUp();
   
   delete ip_AppHelper;
   
   // This allows us to choose the best AbstractSubsequenceHelper based upon the current status
   MakeSubsequences(false);
}

void  SierpinskiRieselApp::MakeSubsequences(bool newSieve)
{
   uint64_t max_k = 0;
   uint64_t max_c = 0;
   uint64_t c;
   uint64_t largestPrimeTested = GetLargestPrimeTested(false);
   uint32_t seqIdx;
   
   seqIdx = 0;
   while (seqIdx < ii_SequenceCount)
   {
      bool haveTerm = false;
      
      for (uint32_t n=ii_MinN; n<=ii_MaxN; n++)
      {
         if (ip_Sequences[seqIdx].nTerms[NBIT(n)])
         {
            haveTerm = true;
            break;
         }
      }
      
      if (haveTerm)
      {
         seqIdx++;
         continue;
      }

      ip_Sequences[seqIdx].nTerms.clear();
      
      WriteToConsole(COT_OTHER, "Sequence %" PRIu64"*%u^n%+" PRId64" removed as all terms have a factor", 
                     ip_Sequences[seqIdx].k, ii_Base, ip_Sequences[seqIdx].c);
            
      for (uint32_t seqIdx2=seqIdx; seqIdx2+1<ii_SequenceCount; seqIdx2++)
         memcpy(&ip_Sequences[seqIdx2], &ip_Sequences[seqIdx2+1], sizeof(seq_t));
      
      ii_SequenceCount--;
   }

   if (ii_SequenceCount == 0)
      FatalError("All sequences have been removed");
   
   for (seqIdx=0; seqIdx<ii_SequenceCount; seqIdx++)
   {
      c = abs(ip_Sequences[seqIdx].c);
      
      max_k = MAX(max_k, ip_Sequences[seqIdx].k);
      max_c = MAX(max_c, c);
   }

   if (max_c > 1 && ib_UseLengendreTables)
   {
      WriteToConsole(COT_OTHER, "Cannot use Legendre tables because at least one sequence has abs(c) > 1");
      ib_UseLengendreTables = false;
      is_LegendreFileName = "";
   }
   
   il_SmallPrimeSieveLimit = 0;
   
   // 65536 is from srsieve.  I don't understand the limit, but if this
   // value is too large, then factors are missed.
   il_SmallPrimeSieveLimit = MAX(257, ii_Base);
   
   if (max_k < 65536 && max_c < 65536)
   {
      il_SmallPrimeSieveLimit = MAX(il_SmallPrimeSieveLimit, max_k);
      il_SmallPrimeSieveLimit = MAX(il_SmallPrimeSieveLimit, max_c);
   }
   
   // The constructors will make a copy of the sequences.
   if (newSieve || max_c > 1 || max_k > largestPrimeTested)
      ;
   
      ip_AppHelper = new GenericSubsequenceHelper(this, ip_Sequences, ii_SequenceCount);
   //else
   //{
   //   bool builtLegendreTables = false;
   //  
   //   if (ib_UseLengendreTables)
   //   {
   //      CIsOneSubsequenceHelper *lah = new CIsOneSubsequenceHelper(this, ip_Sequences, ii_SequenceCount, ib_UseLengendreTables, is_LegendreFileName);
   //    
   //      builtLegendreTables = lah->BuildLegendreTables();
   //}
   
   if (newSieve)
   {
      uint64_t termCount = ip_AppHelper->MakeSubsequencesForNewSieve(ip_Sequences);
      
      if (termCount != il_TermCount)
         FatalError("Expected %" PRIu64" terms, but only set %" PRIu64"", il_TermCount, termCount);
   }
   else
      ip_AppHelper->MakeSubsequencesForOldSieve(ip_Sequences, il_TermCount);
}

void     SierpinskiRieselApp::ReportFactor(uint64_t thePrime, uint32_t seqIdx, uint32_t n)
{
   uint32_t nbit;
               
   if (n < ii_MinN || n > ii_MaxN)
      return;

   if (thePrime > GetMaxPrimeForSingleWorker())
      ip_FactorAppLock->Lock();

   nbit = NBIT(n);
      
   if (ip_Sequences[seqIdx].nTerms[nbit])
   {
      VerifyFactor(thePrime, seqIdx, n);
      
      il_TermCount--;
      il_FactorCount++;
      ip_Sequences[seqIdx].nTerms[nbit] = false;
            
      if (IsPrime(thePrime, ip_Sequences[seqIdx].k, n, ip_Sequences[seqIdx].c))
      {
         WriteToConsole(COT_OTHER, "%" PRIu64"*%u^%u%+" PRId64" is prime!", ip_Sequences[seqIdx].k, ii_Base, n, ip_Sequences[seqIdx].c);

         WriteToLog("%" PRIu64"*%u^%u%+" PRId64" is prime!", ip_Sequences[seqIdx].k, ii_Base, n, ip_Sequences[seqIdx].c);
      }
      else
         LogFactor(thePrime, "%" PRIu64"*%u^%u%+" PRId64"", ip_Sequences[seqIdx].k, ii_Base, n, ip_Sequences[seqIdx].c);
   }

   if (thePrime > GetMaxPrimeForSingleWorker())
      ip_FactorAppLock->Release();
}

bool  SierpinskiRieselApp::IsPrime(uint64_t p, uint64_t k, uint32_t n, int64_t c)
{
   if (n > 64 || (c > 0 && p <= (uint64_t) c))
      return false;

   // p = k*b^n+c --> p = k*b^n
   p -= c;

   if (p % k > 0)
      return false;
   
   // p = k*b^n --> p = b^n
   p /= k;

   // keep dividing by b until p%b != 0
   for ( ; p % ii_Base == 0; p /= ii_Base)
      n--;

   // if p != 1 and n != 0, then p != b^n
   return (n == 0 && p == 1);
}

void  SierpinskiRieselApp::VerifyFactor(uint64_t thePrime, uint32_t seqIdx, uint32_t n)
{
   uint64_t  rem;

   fpu_push_1divp(thePrime);
   
   rem = fpu_powmod(ii_Base, n, thePrime);
   rem = fpu_mulmod(rem, ip_Sequences[seqIdx].k, thePrime);

   if (ip_Sequences[seqIdx].c < 0 && rem < (uint64_t) -ip_Sequences[seqIdx].c)
      rem += thePrime;
   
   rem += ip_Sequences[seqIdx].c;
   
   if (rem >= thePrime)
      rem -= thePrime;

   if (rem != 0)
      FatalError("%" PRIu64"*%u^%u%+" PRId64" mod %" PRIu64" = %" PRIu64"", ip_Sequences[seqIdx].k, ii_Base, n, ip_Sequences[seqIdx].c, thePrime, rem);
   
   fpu_pop();
}