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
#include "../core/BigHashTable.h"
#include "../primesieve/include/primesieve.hpp"
#include "SierpinskiRieselApp.h"
#include "AlgebraicFactorHelper.h"
#include "../core/MpArith.h"

#include "GenericSequenceHelper.h"
#include "CisOneWithOneSequenceHelper.h"
#include "CisOneWithMultipleSequencesHelper.h"

#define APP_VERSION     "1.8.8"

#if defined(USE_OPENCL)
#define APP_NAME        "srsieve2cl"
#elif defined(USE_METAL)
#define APP_NAME        "srsieve2mtl"
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
   ii_MaxD = 0;
   
   ip_FirstSequence = NULL;
   ii_SequenceCount = 0;
   ii_SquareFreeB = 0;
   is_LegendreDirectoryName = "";
   is_SequencesToRemove = "";
   il_LegendreTableBytes = PMAX_MAX_62BIT;

   ii_BaseMultipleMultiplier = 0;
   ii_PowerResidueLcmMulitplier = 0;
   ii_LimitBaseMultiplier = 0;
   id_GiantStepFactor = 1.0;
   ib_ShowQEffort = false;
   ii_UserBestQ = 0;
   ib_SplitByBestQ = false;
   ib_HaveGenericWorkers = false;
   ib_RemoveN = false;
   ib_Algebraic = false;
   ib_OnlyPrimeNs = false;
   ib_UseGenericLogic = false;
   
#if defined(USE_OPENCL) || defined(USE_METAL)
   ib_UseGPUWorkersUponRebuild = false;
   ii_GpuFactorDensity = 100;
   ii_KernelCount = 1;
#endif
}

SierpinskiRieselApp::~SierpinskiRieselApp(void)
{
    if (ip_AppHelper != NULL)
      delete ip_AppHelper;
}

void SierpinskiRieselApp::Help(void)
{
   FactorApp::ParentHelp();

   printf("-n --nmin=n           Minimum n to search (append value with 'p' to limit to n that are prime)\n");
   printf("-N --nmax=N           Maximum n to search\n");
   printf("-s --sequence=s       Sequence in form k*b^n+c where k, b, and c are decimal values\n");
   printf("-f --format=f         Format of output file (A=ABC, D=ABCD (default), B=BOINC, P=ABC with number_primes, M=mfakt)\n");
   printf("-l --legendrebytes=l  Bytes to use for Legendre tables (only used if abs(c)=1 for all sequences)\n");
   printf("-L --legendrefile=L   Input/output diretory for Legendre tables (no files if -L not specified or -l0 is used)\n");
   printf("-Q --showqcost        Output estimated effort for each q\n");
   printf("-q --useq=q           q to use for discrete log\n");
   printf("-r --removen          For sequences with d > 1, remove n where k*b^n+/-c mod d != 0\n");
   printf("-R --remove=r         Remove single sequence r or sequences specified in file r\n");
   printf("-S --splitbybestq     Split sequences into a file based upon best q for each sequence\n");
   printf("-a --algebraic        Exit after generating algebraic factors\n");
   
#if defined(USE_OPENCL) || defined(USE_METAL)
   printf("-M --maxfactordensity=M   factors per 1e6 terms per GPU worker chunk (default %u)\n", ii_GpuFactorDensity);
   printf("-K --kernelcount=K        the number of kernels when splitting large numbers of sequences for the GPU (default %u)\n", ii_KernelCount);
#endif
   
   printf("-c --genericlogic     use generic logic even if abs(c) = 1 for all sequences\n");

   printf("-F --giantstepfactor=F a multiplier used in the calculation giant steps\n");
   printf("                      As F increases, so do the number of giant steps.  default %lf\n", id_GiantStepFactor);

   printf("-U --bmmulitplier=U   multiplied by 2 to compute BASE_MULTIPLE (default %u for single %u for multi\n", 
            DEFAULT_BM_MULTIPLIER_SINGLE, DEFAULT_BM_MULTIPLIER_MULTI);
   printf("                      default BASE_MULTIPLE=%u, BASE_MULTIPLE=%u for multi)\n", 
            DEFAULT_BM_MULTIPLIER_SINGLE * 2, DEFAULT_BM_MULTIPLIER_MULTI * 2);

   printf("-V --prmmultiplier=V  multiplied by BASE_MULTIPLE to compute POWER_RESIDUE_LCM (default %u for single %u for multi\n",
            DEFAULT_PRL_MULTIPLIER_SINGLE, DEFAULT_PRL_MULTIPLIER_MULTI);
   printf("                      default POWER_RESIDUE_LCM=%u, POWER_RESIDUE_LCM=%u for multi)\n", 
            DEFAULT_PRL_MULTIPLIER_SINGLE * DEFAULT_BM_MULTIPLIER_SINGLE, DEFAULT_PRL_MULTIPLIER_MULTI * DEFAULT_BM_MULTIPLIER_MULTI);
            
   printf("-X --lbmultipler=X    multiplied by POWER_RESIDUE_LCM to compute LIMIT_BASE  (default %u for single %u for multi\n",
            DEFAULT_LB_MULTIPLIER_SINGLE, DEFAULT_LB_MULTIPLIER_MULTI);
   printf("                      default LIMIT_BASE=%u, LIMIT_BASE=%u for multi)\n", 
            DEFAULT_LB_MULTIPLIER_SINGLE * DEFAULT_PRL_MULTIPLIER_SINGLE, DEFAULT_LB_MULTIPLIER_MULTI * DEFAULT_PRL_MULTIPLIER_MULTI);
}

void  SierpinskiRieselApp::AddCommandLineOptions(std::string &shortOpts, struct option *longOpts)
{
   FactorApp::ParentAddCommandLineOptions(shortOpts, longOpts);

   shortOpts += "acn:N:s:f:l:L:Qq:rR:U:V:X:F:S";

   AppendLongOpt(longOpts, "nmin",            required_argument, 0, 'n');
   AppendLongOpt(longOpts, "nmax",            required_argument, 0, 'N');
   AppendLongOpt(longOpts, "sequence",        required_argument, 0, 's');
   AppendLongOpt(longOpts, "format",          required_argument, 0, 'f');
   AppendLongOpt(longOpts, "legendrebytes",   required_argument, 0, 'l');
   AppendLongOpt(longOpts, "legendrefile",    required_argument, 0, 'L');
   AppendLongOpt(longOpts, "useq",            required_argument, 0, 'q');
   AppendLongOpt(longOpts, "showqcost",       required_argument, 0, 'Q');
   AppendLongOpt(longOpts, "remove",          required_argument, 0, 'R');
   AppendLongOpt(longOpts, "giantstepfactor", required_argument, 0, 'F');
   AppendLongOpt(longOpts, "basemultiple",    required_argument, 0, 'U');
   AppendLongOpt(longOpts, "limitbase",       required_argument, 0, 'V');
   AppendLongOpt(longOpts, "powerresidue",    required_argument, 0, 'X');
   AppendLongOpt(longOpts, "splitbybestq",    no_argument,       0, 'S');
   AppendLongOpt(longOpts, "algebraic",       no_argument,       0, 'a');
   AppendLongOpt(longOpts, "genericlogic",    no_argument,       0, 'c');

#if defined(USE_OPENCL) || defined(USE_METAL)
   shortOpts += "M:K:";
   
   AppendLongOpt(longOpts, "maxfactordensity",  required_argument, 0, 'M');
   AppendLongOpt(longOpts, "kernelcount",       required_argument, 0, 'K');
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
         if (arg[strlen(arg)-1] == 'P' || arg[strlen(arg)-1] == 'p')
         {
            arg[strlen(arg)-1] = 0;
            ib_OnlyPrimeNs = true;
         }
         
         status = Parser::Parse(arg, 1, NMAX_MAX, ii_MinN);
         break;

      case 'N':
         status = Parser::Parse(arg, 1, NMAX_MAX, ii_MaxN);
         break;

      case 'l':
         status = Parser::Parse(arg, 0, PMAX_MAX_62BIT, il_LegendreTableBytes);
         break;
         
      case 'L':
         is_LegendreDirectoryName = arg;
         status = P_SUCCESS;
         break;
         
      case 'f':
         char value;
         status = Parser::Parse(arg, "ABDPM", value);
         
         it_Format = FF_UNKNOWN;
   
         if (value == 'A')
            it_Format = FF_ABC;
         if (value == 'B')
            it_Format = FF_BOINC;
         if (value == 'D')
            it_Format = FF_ABCD;
         if (value == 'P')
            it_Format = FF_NUMBER_PRIMES;
         if (value == 'M')
            it_Format = FF_MFAKT;
         break;
         
      case 's':
         if (!LoadSequencesFromFile(arg))
            ValidateAndAddNewSequence(arg);

         status = P_SUCCESS;
         break;

      case 'Q':
         ib_ShowQEffort = true;
         status = P_SUCCESS;
         break;
         
      case 'q':
         status = Parser::Parse(arg, 1, 1<<15, ii_UserBestQ);
         status = P_SUCCESS;
         break;
         
      case 'R':
         is_SequencesToRemove = arg;
         status = P_SUCCESS;
         break;
         
      case 'r':
         ib_RemoveN = true;
         status = P_SUCCESS;
         break;

      case 'F':
         sscanf(arg, "%lf", &id_GiantStepFactor);
         status = P_SUCCESS;
         break;

      case 'a':
         ib_Algebraic = true;
         status = P_SUCCESS;
         break;

      case 'c':
         ib_UseGenericLogic = true;
         status = P_SUCCESS;
         break;
      
      case 'S':
         ib_SplitByBestQ = true;
         status = P_SUCCESS;
         break;
         
      case 'U':
         status = Parser::Parse(arg, 1, 50, ii_BaseMultipleMultiplier);
         break;

      case 'V':
         status = Parser::Parse(arg, 1, 400, ii_PowerResidueLcmMulitplier);
         break;

      case 'X':
         status = Parser::Parse(arg, 1, 10, ii_LimitBaseMultiplier);
         break;

#if defined(USE_OPENCL) || defined(USE_METAL)
      case 'M':
         status = Parser::Parse(arg, 1, 1000000, ii_GpuFactorDensity);
         break;
         
      case 'K':
         status = Parser::Parse(arg, 1, 1000000, ii_KernelCount);
         break;
#endif
   }

   return status;
}

void SierpinskiRieselApp::ValidateOptions(void)
{
   seq_t     *seqPtr;

   if (it_Format == FF_UNKNOWN)
      FatalError("the specified file format in not valid, use A (ABC), D (ABCD), P (ABC with number_primes), or B (BOINC)");
   
   if (is_InputTermsFileName.length() > 0)
   {
      if (ib_HaveNewSequences)
         FatalError("cannot add new candidate sequences in to an existing sieve");
         
      ProcessInputTermsFile(false);
         
      seqPtr = ip_FirstSequence;
      do
      {
         seqPtr->nTerms.resize(ii_MaxN - ii_MinN + 1);
         std::fill(seqPtr->nTerms.begin(), seqPtr->nTerms.end(), false);

         seqPtr = (seq_t *) seqPtr->next;
      } while (seqPtr != NULL);

      ProcessInputTermsFile(true);
   
      RemoveSequences();
      
      RemoveN();
      
      if (!ib_ApplyAndExit)
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

      if (GetMinPrime() != 3)
         FatalError("pmin cannot be overridden when starting a new sieve");
         
      AlgebraicFactorHelper *afh = new AlgebraicFactorHelper(this, ii_Base, ii_MinN, ii_MaxN);
      
      seqPtr = ip_FirstSequence;
      if ((seq_t *) seqPtr->next == NULL)
      {
         if (seqPtr->k == 1 && seqPtr->d > 1 && ii_Base == seqPtr->d+1)
         {
            if (ii_MinN == 2)
               ii_MinN = 3;
            
            ib_OnlyPrimeNs = true;
            WriteToConsole(COT_OTHER, "With GRUs, composite n will have factors, so removing them");
         }
      }
      
      std::vector<uint64_t> primes;
   
      if (ib_OnlyPrimeNs)
         primesieve::generate_primes(ii_MinN, ii_MaxN, &primes);
      
      do
      {
         seqPtr->nTerms.resize(ii_MaxN - ii_MinN + 1);
         
         if (ib_OnlyPrimeNs)
         {
            std::fill(seqPtr->nTerms.begin(), seqPtr->nTerms.end(), false);
   
            std::vector<uint64_t>::iterator it = primes.begin();
                        
            while (it != primes.end())
            {
               uint64_t prime = *it;
               it++;
            
               seqPtr->nTerms[NBIT(prime)] = true;
               il_TermCount++;
            }
         }
         else
         {
            std::fill(seqPtr->nTerms.begin(), seqPtr->nTerms.end(), true);
            
            il_TermCount += (ii_MaxN - ii_MinN + 1);
         }
            
         il_TermCount -= afh->RemoveTermsWithAlgebraicFactors(seqPtr);
         
         seqPtr = (seq_t *) seqPtr->next;
      } while (seqPtr != NULL);
      
      delete afh;
         
      if (ib_Algebraic)
         FatalError("Exiting without sieving after generating algebraic factors");
   
      RemoveN();
   
      MakeSubsequences(true, GetMinPrime());  
   }

   seqPtr = ip_FirstSequence;
   int64_t firstC = seqPtr->c;
   ib_HaveSingleC = true;
   
   do
   {
      if (firstC != seqPtr->c)
         ib_HaveSingleC = false;      
      
      if (it_Format == FF_BOINC)
      {
         if (seqPtr->d != 1)
            FatalError("When using BOINC format, all sequences must have d=1");
         
         if (abs(seqPtr->c) != 1)
            FatalError("When using BOINC format, all sequences must have c=+1 or c=-1");
         
         if (seqPtr->c != ip_FirstSequence->c)
            FatalError("When using BOINC format, cannot mix c=+1 and c=-1 sequences");
      }

      if (it_Format == FF_MFAKT)
      {
         if (seqPtr->k != 1 || seqPtr->d == 1 || ii_Base != seqPtr->d+1 || seqPtr != ip_FirstSequence)
            FatalError("MFAKT format is only support for GRUs");

         if (ii_MinN < 50000)
            FatalError("MAKTC format requires minN >= 50000");
      }

      seqPtr = (seq_t *) seqPtr->next;
   } while (seqPtr != NULL);

   if (is_OutputTermsFileName.length() == 0)
   {
      char  fileName[30];
      
      if (it_Format == FF_ABCD)
         snprintf(fileName, sizeof(fileName), "b%u_n.abcd", ii_Base);
      if (it_Format == FF_ABC)
         snprintf(fileName, sizeof(fileName), "b%u_n.abc", ii_Base);
      if (it_Format == FF_BOINC)
         snprintf(fileName, sizeof(fileName), "b%u_n.boinc", ii_Base);
      if (it_Format == FF_NUMBER_PRIMES)
         snprintf(fileName, sizeof(fileName), "b%u_n.abcnp", ii_Base);
      if (it_Format == FF_MFAKT)
         snprintf(fileName, sizeof(fileName), "worktodo_b%u.txt", ii_Base);
      
      is_OutputTermsFileName = fileName;
   }
   
   if (il_LegendreTableBytes == 0 && is_LegendreDirectoryName.length() > 0)
   {
      WriteToConsole(COT_OTHER, "Ingoring -L option since Legendre tables cannot be used");
      is_LegendreDirectoryName = "";
   }
   
   // Allow only one worker to do work when processing small primes.  This allows us to avoid 
   // locking when factors are reported, which significantly hurts performance as most terms 
   // will be removed due to small primes.
   
   // This will sieve beyond the limit, but we want to make sure that at least one prime
   // larger than this limit is passed to the worker even if the worker does not test it.
   SetMaxPrimeForSingleWorker(il_SmallPrimeSieveLimit + 1000);

#if defined(USE_OPENCL) || defined(USE_METAL)
   SetMinGpuPrime(1000000);
   
   double factors = (double) (ii_MaxN - ii_MinN) * (double) (ii_SequenceCount) / 1000000.0;

   ii_MaxGpuFactors = GetGpuWorkGroups() * (uint64_t) (factors * (double) ii_GpuFactorDensity);
   
   if (ii_MaxGpuFactors < 10) ii_MaxGpuFactors = 10;
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
   char     temp[500];

   if (!StripCRLF(arg))
      return;

   if (sscanf(arg, "(%" SCNu64"*%u^n%" SCNd64")/%u", &k, &b, &c, &d) == 4)
   {
      snprintf(temp, sizeof(temp), "(%" PRIu64"*%u^n%+" PRId64")/%u", k, b, c, d);
   }
   else if (sscanf(arg, "(%u^n%" SCNd64")/%u", &b, &c, &d) == 3)
   {
      k = 1;
      
      snprintf(temp, sizeof(temp), "(%u^n%+" PRId64")/%u", b, c, d);
   }
   else if (sscanf(arg, "%" SCNu64"*%u^n%" SCNd64"", &k, &b, &c) == 3)
   {
      d = 1;
      snprintf(temp, sizeof(temp), "%" PRIu64"*%u^n%+" PRId64"", k, b, c);
   }
   else
      FatalError("sequence %s must be in form k*b^n+c, (k*b^n+c)/d, or (b^n+c)/d where you specify values for k, b, c, and d", arg);

   if (strcmp(arg, temp))
      FatalError("Input sequence %s has extraneous trailing characters", arg);
      
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
   char     buffer[1000], *pos;
   uint32_t n, diff;
   uint64_t k, prevK = 0;
   int64_t  c, prevC = 0;
   uint32_t d = 1, prevD = 1;
   uint64_t lastPrime = 0;
   uint32_t lineNumber = 0;
   format_t format = FF_UNKNOWN;
   seq_t   *currentSequence = 0;
   bool     haveMinN = false;
   
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
         pos = strstr(buffer, "// Sieved to");
         
         if (pos != NULL)
         {
            if (sscanf(pos, "// Sieved to %" SCNu64"", &lastPrime) != 1)
               FatalError("Line %u is not a valid ABCD line in input file %s", lineNumber, is_InputTermsFileName.c_str());
            *pos = 0;
         }
         
         d = 1;
         k = 1;
         
         // Either (k*b^$a+c)/d [n] or k*b^$a+c [n] or (b^a$+c)/d [n]
         if (sscanf(buffer, "ABCD (%" SCNu64"*%u^$a%" SCNd64")/%u [%u]", &k, &ii_Base, &c, &d, &n) != 5)
         {
            d = 1;
            k = 1;
            
            if (sscanf(buffer, "ABCD %" SCNu64"*%u^$a%" SCNd64" [%u]", &k, &ii_Base, &c, &n) != 4)
            {
               d = 1;
               k = 1;
               
               if (sscanf(buffer, "ABCD (%u^$a%" SCNd64")/%u [%u]", &ii_Base, &c, &d, &n) != 4)
                  FatalError("Line %u is not a valid ABCD line in input file %s", lineNumber, is_InputTermsFileName.c_str());
            }
         }
         
         format = FF_ABCD;

         if (haveBitMap)
         {
            currentSequence = GetSequence(k, c, d, currentSequence);
            currentSequence->nTerms[NBIT(n)] = true;
            il_TermCount++;
            continue;
         }
         else
         {            
            if (!haveMinN)
            {
               ii_MinN = ii_MaxN = n;
               haveMinN = true;
            }
            
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
         pos = strstr(buffer, "// Sieved to");
         
         if (pos != NULL)
         {
            if (sscanf(pos, "// Sieved to %" SCNu64"", &lastPrime) != 1)
               FatalError("Line %u is not a valid ABCD line in input file %s", lineNumber, is_InputTermsFileName.c_str());
            *pos = 0;
         }
         
         d = 1;
         k = 1;
         
         // Either (k*b^$a+c)/d or k*b^$a+c or (b^a$+c)/d
         if (sscanf(buffer, "ABC (%" SCNu64"*%u^$a%" SCNd64")/%u", &k, &ii_Base, &c, &d) != 4)
         {
            d = 1;
            k = 1;
            
            if (sscanf(buffer, "ABC %" SCNu64"*%u^$a%" SCNd64"", &k, &ii_Base, &c) != 3)
            {
               d = 1;
               k = 1;

               if (sscanf(buffer, "ABC (%u^$a%" SCNd64")/%u", &ii_Base, &c, &d) != 3)
                  FatalError("Line %u is not a valid ABC line in input file %s", lineNumber, is_InputTermsFileName.c_str());
            }
         }

         format = FF_ABC;
         
         if (k != prevK || c != prevC || d != prevD)
         {
            if (haveBitMap)
               currentSequence = GetSequence(k, c, d, currentSequence);
            else
               AddSequence(k, c, d);
         }
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
               
               break;

            case FF_BOINC:
               if (sscanf(buffer, "%" SCNu64" %u", &k, &n) != 2)
                  FatalError("Line %s is malformed", buffer);
               
               if (k != prevK || c != prevC || d != prevD)
               {
                  if (haveBitMap)
                     currentSequence = GetSequence(k, c, d, currentSequence);
                  else
                     AddSequence(k, c, d);
               }
               
               break;
               
            case FF_NUMBER_PRIMES:
               if (sscanf(buffer, "%" SCNu64" %u %" SCNd64"", &k, &n, &c) != 3)
                  FatalError("Line %s is malformed", buffer);
               
               if (k != prevK || c != prevC || d != prevD)
               {
                  if (haveBitMap)
                     currentSequence = GetSequence(k, c, d, currentSequence);
                  else
                     AddSequence(k, c, d);
               }
               
               break;
               
            default:
               FatalError("Input file %s has unknown format", is_InputTermsFileName.c_str());
         }

         prevK = k;
         prevC = c;
         prevD = d;

         if (haveBitMap)
         {
            currentSequence->nTerms[NBIT(n)] = true;
            il_TermCount++;
         }
         else
         {
            if (!haveMinN)
            {
               ii_MinN = ii_MaxN = n;
               haveMinN = true;
            }
         
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
   seq_t   *seqPtr, *nextSeq, *prevSeq;
   char     sequence[100];
   bool     found = false;
   uint32_t removedCount = 0;
   
   if (d == 1)
      snprintf(sequence, sizeof(sequence), "%" PRIu64"*%u^n%+" PRId64"", k, b, c);
   else
      snprintf(sequence, sizeof(sequence), "(%" PRIu64"*%u^n%+" PRId64")/%u", k, b, c, d);
   
   if (b != ii_Base)
      WriteToConsole(COT_OTHER, "Sequence %s wasn't removed because it is not base %u", sequence, b);
   
   prevSeq = seqPtr = ip_FirstSequence;
   do
   {      
      nextSeq = (seq_t *) seqPtr->next;
      
      if (seqPtr->k == k && seqPtr->c == c && seqPtr->d == d)
      {
         found = true;
         
         for (uint32_t n=ii_MinN; n<=ii_MaxN; n++)
         {
            if (seqPtr->nTerms[NBIT(n)])
               removedCount++;
         }

         seqPtr->nTerms.clear();
                     
         if (seqPtr == ip_FirstSequence)
            ip_FirstSequence = (seq_t *) seqPtr->next;
         else
            prevSeq->next = seqPtr->next;
            
         WriteToConsole(COT_OTHER, "%u terms for sequence %s have been removed", removedCount, sequence);
                  
         xfree(seqPtr);
         
         il_TermCount -= removedCount;
         ii_SequenceCount--;
      }
      
      prevSeq = seqPtr;
      seqPtr = nextSeq;
   } while (seqPtr != NULL);

   UpdateSequenceIndexes();
   
   if (!found)
      WriteToConsole(COT_OTHER, "Sequence %s wasn't removed because it wasn't found", sequence);
}

void SierpinskiRieselApp::UpdateSequenceIndexes(void)
{
   uint32_t  ii_SequenceCount = 0;
   seq_t    *seqPtr = ip_FirstSequence;
   
   while (seqPtr != NULL)
   {
      ii_SequenceCount++;
      seqPtr->seqIdx = ii_SequenceCount;
      
      seqPtr = (seq_t *) seqPtr->next;
   };
}

void SierpinskiRieselApp::RemoveN(void)
{
   if (!ib_RemoveN)
      return;

   seq_t *seqPtr = ip_FirstSequence;
   do
   {
      if (seqPtr->d == 1)
      {
         seqPtr = (seq_t *) seqPtr->next;
         continue;
      }
    
      uint32_t d = seqPtr->d;

      MpArith mp(d);

      int64_t c = seqPtr->c;
      
      while (c < 0)
         c += d;
      
      c %= d;
      
      MpRes mpK = mp.nToRes(seqPtr->k % d);
      MpRes mpB = mp.nToRes(ii_Base % d);
      MpRes mpRem = mp.pow(mpB, ii_MinN);
      
      mpRem = mp.mul(mpRem, mpK);
      
      uint64_t remKBN = mp.resToN(mpRem);
      uint64_t rem;      
      uint32_t removed = 0;
      uint32_t n = ii_MinN;
      
      while (n <= ii_MaxN && n < ii_MinN+(d-1))
      {
         rem = remKBN + c;
         
         if (rem >= d)
            rem -= d;
         
         if (rem != 0)
         {
            // if k*b^n+/-c % d != 0, then we know that k*b^(n+m*(d-1))+/-c % d != 0
            for (uint32_t tempN=n; tempN<=ii_MaxN; tempN+=(d-1))
            {
               if (seqPtr->nTerms[NBIT(tempN)])
                  {
                     seqPtr->nTerms[NBIT(tempN)] = false;
                     il_TermCount--;
                     removed++;
                  }
            }
         }
         
         remKBN = mp.mul(mpRem, mpB);
         n++;
      }

      if (removed > 0)
         WriteToConsole(COT_OTHER, "Removed %u terms for (%" PRIu64"*%u^n%+" PRId64")/%u because gcd(%" PRIu64"*%u^n%+" PRId64", %u) != 0",
                        removed, seqPtr->k, ii_Base, seqPtr->c, seqPtr->d, seqPtr->k, ii_Base, seqPtr->c, seqPtr->d);

      seqPtr = (seq_t *) seqPtr->next;
   } while (seqPtr != NULL);
}

bool SierpinskiRieselApp::ApplyFactor(uint64_t theFactor, const char *term)
{
   uint64_t   k;
   uint32_t   b, n, d;
   int64_t    c;
   seq_t     *seqPtr;
   
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
   
   seqPtr = ip_FirstSequence;
   do
   {
      if (seqPtr->k == k && seqPtr->c == c && seqPtr->d == d)
      {
         VerifyFactor(theFactor, seqPtr, n);
            
         if (seqPtr->nTerms[NBIT(n)])
         {
            seqPtr->nTerms[NBIT(n)] = false;
            il_TermCount--;
         
            return true;
         }
      }

      seqPtr = (seq_t *) seqPtr->next;
   } while (seqPtr != NULL);

   return false;
}

void SierpinskiRieselApp::WriteOutputTermsFile(uint64_t largestPrime)
{
   id_EstimatedPrimes = 0;
   
   if (ib_SplitByBestQ)
      return;
   
   WriteOutputTermsFile(largestPrime, 0);
}

void SierpinskiRieselApp::OuptutAdditionalConsoleMessagesUponFinish(void)
{
   WriteToConsole(COT_OTHER, "%.02lf estimated primes", id_EstimatedPrimes);
}

void SierpinskiRieselApp::WriteOutputTermsFilesByQ(void)
{
   uint64_t termsCounted = 0;
   uint32_t maxQ = 0;
   
   seq_t *seqPtr = ip_FirstSequence;
   do
   {
      if (maxQ < seqPtr->bestQ)
         maxQ = seqPtr->bestQ; 
      seqPtr = (seq_t *) seqPtr->next;
   } while (seqPtr != NULL);
   
   if (is_InputTermsFileName.length() == 0)
   {
      for (uint32_t q=1; q<=maxQ; q++)
      {
         seqPtr = ip_FirstSequence;
         
         do
         {
            if (seqPtr->bestQ == q)
            {
               WriteSequenceFilesByQ(q);
               break;
            }
            
            seqPtr = (seq_t *) seqPtr->next;
         } while (seqPtr != NULL);
      }
   }
   else 
   {
      for (uint32_t q=1; q<=maxQ; q++)
      {
         seqPtr = ip_FirstSequence;
         
         do
         {
            if (seqPtr->bestQ == q)
            {
               termsCounted += WriteOutputTermsFile(GetMinPrime(), q);
               break;
            }
            
            seqPtr = (seq_t *) seqPtr->next;
         } while (seqPtr != NULL);
      }
      
      if (termsCounted != il_TermCount)
         FatalError("Something is wrong.  Counted terms (%" PRIu64") != expected terms (%" PRIu64")", termsCounted, il_TermCount);
   }
}

void  SierpinskiRieselApp::WriteSequenceFilesByQ(uint32_t q)
{
   seq_t   *seqPtr;
   char     fileName[200];
   uint32_t sequenceCount = 0;

   snprintf(fileName, sizeof(fileName), "q%03u_b%u.in", q, ii_Base);
   
   FILE    *termsFile = fopen(fileName, "w");

   if (!termsFile)
      FatalError("Unable to open output file %s", is_OutputTermsFileName.c_str());
   
   seqPtr = ip_FirstSequence;
   do
   {
      if (seqPtr->bestQ == q)
      {
         sequenceCount++;

         if (seqPtr->d > 1)
            fprintf(termsFile, "(%" PRIu64"*%u^n%+" PRId64")/%u\n", seqPtr->k, ii_Base, seqPtr->c, seqPtr->d);
         else
            fprintf(termsFile, "%" PRIu64"*%u^n%+" PRId64"\n", seqPtr->k, ii_Base, seqPtr->c);
      }
            
      seqPtr = (seq_t *) seqPtr->next;
   } while (seqPtr != NULL);
   
   fclose(termsFile);
 
   WriteToConsole(COT_OTHER, "%5u sequences for q %3u written to %s", sequenceCount, q, fileName); 
}

uint32_t SierpinskiRieselApp::WriteOutputTermsFile(uint64_t largestPrime, uint32_t q)
{
   uint64_t termsCounted = 0;
   bool     allSequencesHaveDEqual1 = true;
   seq_t   *seqPtr;
   char     fileName[200];
   uint32_t sequenceCount = 0;
   
   // With super large ranges, wait until we can lock because without locking
   // the term count can change between opening and closing the file.
   if (IsRunning() && largestPrime < GetMaxPrimeForSingleWorker())
      return 0;
   
   if (q == 0)
      snprintf(fileName, sizeof(fileName), "%s", is_OutputTermsFileName.c_str());
   else
      snprintf(fileName, sizeof(fileName), "q%03u_%s", q, is_OutputTermsFileName.c_str());
   
   FILE    *termsFile = fopen(fileName, "w");

   if (!termsFile)
      FatalError("Unable to open output file %s", is_OutputTermsFileName.c_str());
   
   ip_FactorAppLock->Lock();
   
   termsCounted = 0;

   if (it_Format == FF_NUMBER_PRIMES)
   {
      seqPtr = ip_FirstSequence;
      do
      {
         if (q == 0 || seqPtr->bestQ == q)
         {
            if (seqPtr->d != 1)
               allSequencesHaveDEqual1 = false;
         }

         seqPtr = (seq_t *) seqPtr->next;
      } while (seqPtr != NULL);
         
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
   
   seqPtr = ip_FirstSequence;
   do
   {
      if (q == 0 || seqPtr->bestQ == q)
      {
         sequenceCount++;

         if (it_Format == FF_ABCD)
            termsCounted += WriteABCDTermsFile(seqPtr, largestPrime, termsFile);
         
         if (it_Format == FF_ABC)
            termsCounted += WriteABCTermsFile(seqPtr, largestPrime, termsFile);
         
         if (it_Format == FF_BOINC)
            termsCounted += WriteBoincTermsFile(seqPtr, largestPrime, termsFile);
         
         if (it_Format == FF_NUMBER_PRIMES)
            termsCounted += WriteABCNumberPrimesTermsFile(seqPtr, largestPrime, termsFile, allSequencesHaveDEqual1);

         if (it_Format == FF_MFAKT)
            termsCounted += WriteMfaktTermsFile(seqPtr, largestPrime, termsFile);
      }
            
      seqPtr = (seq_t *) seqPtr->next;
   } while (seqPtr != NULL);
   
   fclose(termsFile);
   
   if (q == 0 && termsCounted != il_TermCount)
      FatalError("Something is wrong.  Counted terms (%" PRIu64") != expected terms (%" PRIu64")", termsCounted, il_TermCount);

   if (q > 0 && sequenceCount > 0)
      WriteToConsole(COT_OTHER, "%5u sequences with %8" PRIu64" terms written to %s", sequenceCount, termsCounted, fileName);

   ip_FactorAppLock->Release();
   
   return termsCounted;
}

uint32_t SierpinskiRieselApp::WriteABCDTermsFile(seq_t *seqPtr, uint64_t maxPrime, FILE *termsFile)
{
   uint32_t n, nCount = 0, previousN;
   uint32_t bit;

   double dCalc = 1.781 * log((double) maxPrime);
   double dBase = (double) log((double) ii_Base);
   double dK = (double) log((double) seqPtr->k);
   double dD = (double) log((double) seqPtr->d);
   double dLength, dTemp;
  
   n = ii_MinN;
   
   bit = NBIT(n);
   for (; n<=ii_MaxN; n++)
   {      
      if (seqPtr->nTerms[bit])
         break;
      
      bit++;
   }

   if (n > ii_MaxN)
      return 0;
   
   if (seqPtr->d == 1)
      fprintf(termsFile, "ABCD %" PRIu64"*%u^$a%+" PRId64" [%u] // Sieved to %" PRIu64"\n", seqPtr->k, ii_Base, seqPtr->c, n, maxPrime);
   else
      fprintf(termsFile, "ABCD (%" PRIu64"*%u^$a%+" PRId64")/%u [%u] // Sieved to %" PRIu64"\n", seqPtr->k, ii_Base, seqPtr->c, seqPtr->d, n, maxPrime);
   
   if (n > 1)
   {
      dLength = dBase * (double) n + dK - dD;
      dTemp = dCalc / dLength;
      id_EstimatedPrimes += (dTemp > 1.0 ? 1.0 : dTemp);
   }
   
   previousN = n;
   nCount = 1;
   n++;
   
   bit = NBIT(n);
   for (; n<=ii_MaxN; n++)
   {
      if (seqPtr->nTerms[bit])
      {
         fprintf(termsFile, "%u\n", n - previousN);
         previousN = n;
         nCount++;
         
         dLength = dBase * (double) n + dK - dD;
         dTemp = dCalc / dLength;
         id_EstimatedPrimes += (dTemp > 1.0 ? 1.0 : dTemp);
      }
      
      bit++;
   }

   return nCount;
}

uint32_t SierpinskiRieselApp::WriteABCTermsFile(seq_t *seqPtr, uint64_t maxPrime, FILE *termsFile)
{
   uint32_t n, nCount = 0;
   uint32_t bit;

   double dCalc = 1.781 * log((double) maxPrime);
   double dBase = (double) log((double) ii_Base);
   double dK = (double) log((double) seqPtr->k);
   double dD = (double) log((double) seqPtr->d);
   double dLength, dTemp;

   if (seqPtr->d == 1)
      fprintf(termsFile, "ABC %" PRIu64"*%u^$a%+" PRId64" // Sieved to %" PRIu64"\n", seqPtr->k, ii_Base, seqPtr->c, maxPrime);
   else
      fprintf(termsFile, "ABC (%" PRIu64"*%u^$a%+" PRId64")/%u // Sieved to %" PRIu64"\n", seqPtr->k, ii_Base, seqPtr->c, seqPtr->d, maxPrime);
   
   n = ii_MinN;
   bit = NBIT(n);

   for (; n<=ii_MaxN; n++)
   {
      if (seqPtr->nTerms[bit])
      {
         fprintf(termsFile, "%u\n", n);
         nCount++;

         if (n > 1)
         {
            dLength = dBase * (double) n + dK - dD;
            dTemp = dCalc / dLength;
            id_EstimatedPrimes += (dTemp > 1.0 ? 1.0 : dTemp);
         }
      }
      
      bit++;
   }
   
   return nCount;
}

uint32_t SierpinskiRieselApp::WriteBoincTermsFile(seq_t *seqPtr, uint64_t maxPrime, FILE *termsFile)
{
   uint32_t n, nCount = 0;
   uint32_t bit;

   double dCalc = 1.781 * log((double) maxPrime);
   double dBase = (double) log((double) ii_Base);
   double dK = (double) log((double) seqPtr->k);
   double dD = (double) log((double) seqPtr->d);
   double dLength, dTemp;
   
   n = ii_MinN;
   bit = NBIT(n);

   for (; n<=ii_MaxN; n++)
   {
      if (seqPtr->nTerms[bit])
      {
         fprintf(termsFile, "%" PRIu64" %u\n", seqPtr->k, n);
         nCount++;

         if (n > 1)
         {
            dLength = dBase * (double) n + dK - dD;
            dTemp = dCalc / dLength;
            id_EstimatedPrimes += (dTemp > 1.0 ? 1.0 : dTemp);
         }      }
      
      bit++;
   }
   
   return nCount;
}

uint32_t SierpinskiRieselApp::WriteABCNumberPrimesTermsFile(seq_t *seqPtr, uint64_t maxPrime, FILE *termsFile, bool allSequencesHaveDEqual1)
{
   uint32_t n, nCount = 0;
   uint32_t bit;

   double dCalc = 1.781 * log(maxPrime);
   double dBase = (double) log(ii_Base);
   double dK = (double) log(seqPtr->k);
   double dD = (double) log(seqPtr->d);
   double dLength, dTemp;

   n = ii_MinN;
   bit = NBIT(n);

   for (; n<=ii_MaxN; n++)
   {
      if (seqPtr->nTerms[bit])
      {
         if (allSequencesHaveDEqual1)
            fprintf(termsFile, "%" PRIu64" %u %+" PRId64"\n", seqPtr->k, n, seqPtr->c);
         else
            fprintf(termsFile, "%" PRIu64" %u %+" PRId64" %u\n", seqPtr->k, n, seqPtr->c, seqPtr->d);
         
         if (n > 1)
         {
            dLength = dBase * (double) n + dK - dD;
            dTemp = dCalc / dLength;
            id_EstimatedPrimes += (dTemp > 1.0 ? 1.0 : dTemp);
         }
         
         nCount++;
      }
      
      bit++;
   }
   
   return nCount;
}

uint32_t  SierpinskiRieselApp::WriteMfaktTermsFile(seq_t *seqPtr, uint64_t maxPrime, FILE *termsFile)
{
   uint32_t n, nCount = 0;
   uint32_t bit, bits;
   
   n = ii_MinN;
   bit = NBIT(n);
   
   double dCalc = 1.781 * log(maxPrime);
   double dBase = (double) log(ii_Base);
   double dK = (double) log(seqPtr->k);
   double dD = (double) log(seqPtr->d);
   double dLength, dTemp;

   bits = 0;
   while (maxPrime > 1)
   {
      bits++;
      maxPrime >>= 1;
   }

   for (; n<=ii_MaxN; n++)
   {
      if (seqPtr->nTerms[bit])
      {
         fprintf(termsFile, "Factor=base=%u,%u,%u,%u\n", ii_Base, n, bits, 50);
         nCount++;

         if (n > 1)
         {
            dLength = dBase * (double) n + dK - dD;
            dTemp = dCalc / dLength;
            id_EstimatedPrimes += (dTemp > 1.0 ? 1.0 : dTemp);
         }
      }
      
      bit++;
   }
   
   return nCount;
}
void  SierpinskiRieselApp::GetExtraTextForSieveStartedMessage(char *extraText, uint32_t maxTextLength)
{
   seq_t  *seqPtr;
   uint64_t minK = 0, maxK = 0;
   int64_t  minC = 0, maxC = 0;
   uint32_t minD = 0, maxD = 0;
   char theC[20], theK[20], theD[20];
   
   seqPtr = ip_FirstSequence;
   do
   {
      if (minK == 0 || minK > seqPtr->k) minK = seqPtr->k;
      if (maxK == 0 || maxK < seqPtr->k) maxK = seqPtr->k;
      
      if (minC == 0 || minC > seqPtr->c) minC = seqPtr->c;
      if (maxC == 0 || maxC < seqPtr->c) maxC = seqPtr->c;
      
      if (minD == 0 || minD > seqPtr->d) minD = seqPtr->d;
      if (maxD == 0 || maxD < seqPtr->d) maxD = seqPtr->d;
         
      seqPtr = (seq_t *) seqPtr->next;
   } while (seqPtr != NULL);
   
   if (minK == maxK)
      snprintf(theK, sizeof(theK), "%" PRIu64"", minK);
   else
      snprintf(theK, sizeof(theK), "k");

   if (minC == maxC)
      snprintf(theC, sizeof(theC), "%+" PRId64"", minC);
   else
   {
      if (minC > 0)
         snprintf(theC, sizeof(theC), "+c");
      else if (maxC < 0)
         snprintf(theC, sizeof(theC), "-c");
      else
         snprintf(theC, sizeof(theC), "+/-c");
   }

   if (minD == maxD)
      snprintf(theD, sizeof(theD), "%u", minD);
   else
      snprintf(theK, sizeof(theK), "d");
   
   if (maxD > 1)
   {            
      if (maxK > 1)
         snprintf(extraText, maxTextLength, "%u <= n <= %u, (%s*%u^n%s)/%s", ii_MinN, ii_MaxN, theK, ii_Base, theC, theD);
      else
         snprintf(extraText, maxTextLength, "%u <= n <= %u, (%u^n%s)/%s", ii_MinN, ii_MaxN, ii_Base, theC, theD);
   }
   else
   {
      if (maxK > 1)
         snprintf(extraText, maxTextLength, "%u <= n <= %u, %s*%u^n%s", ii_MinN, ii_MaxN, theK, ii_Base, theC);
      else
         snprintf(extraText, maxTextLength, "%u <= n <= %u, %u^n%s", ii_MinN, ii_MaxN, ii_Base, theC);
   }
}

void  SierpinskiRieselApp::AddSequence(uint64_t k, int64_t c, uint32_t d)
{
   seq_t   *seqPtr;

   // If base, k, and c are odd then all terms are even
   if ((ii_Base % 2) && (k % 2) && (c % 2) && (d == 1))
   {
      WriteToConsole(COT_OTHER, "Sequence %" PRIu64"*%u^n%+" PRId64" not added because all terms are divisible by 2", k, ii_Base, c);
      return;
   }

   // If either the base or k is even and c is even then all terms are even
   if ((!(ii_Base % 2) || !(k % 2)) && !(c % 2) && (d == 1))
   {
      WriteToConsole(COT_OTHER, "Sequence %" PRIu64"*%u^n%+" PRId64" not added because all terms are divisible by 2", k, ii_Base, c);
      return;
   }
   
   char sequence[50];
      
   // If the given k is greater than any other k, then we can bypas the next check.
   if (k <= il_MaxK && is_SequenceKs.find(k) != is_SequenceKs.end())
   {
      // We have seen the k before, so if is it for the same c and d, then this is
      // a known sequence and we can just return.
      snprintf(sequence, sizeof(sequence), "%" PRIu64"*%u^n%+" PRId64"/%u", k, ii_Base, c, d);

      if (is_SequenceSet.find(sequence) != is_SequenceSet.end())
         return;
   }

   snprintf(sequence, sizeof(sequence), "%" PRIu64"*%u^n%+" PRId64"/%u", k, ii_Base, c, d);

   is_SequenceKs.insert(k);
   is_SequenceSet.insert(sequence);
   
   seq_t *newPtr = (seq_t *) xmalloc(1, sizeof(seq_t), "seq");
   
   uint64_t absc = abs(c);

   if (absc != 1)
   {
      if (!ib_UseGenericLogic)
         WriteToConsole(COT_OTHER, "Must use generic sieving logic because abs(c) != 1 for at least one sequence");

      ib_UseGenericLogic = true;
   }
   
   if (il_MaxK < k) il_MaxK = k;
   if (ii_MaxD < d) ii_MaxD = d;
   if (il_MaxAbsC < absc) il_MaxAbsC = absc;

   newPtr->k = k;
   newPtr->c = c;
   newPtr->d = d;
   newPtr->next = NULL;
   
   if (ii_SequenceCount == 0)
   {
      ip_FirstSequence = newPtr;
      ii_SequenceCount = 1;
      newPtr->seqIdx = ii_SequenceCount;
      
      ip_LastSequence = newPtr;
      return;
   }
   
   ii_SequenceCount++;
   
   if (ip_LastSequence->k <= newPtr->k)
   {
      ip_LastSequence->next = newPtr;
      ip_LastSequence = newPtr;
      newPtr->seqIdx = ii_SequenceCount;
               
      return;
   }

   seqPtr = ip_FirstSequence;
   seq_t *prevSeqPtr = NULL;
   do
   {
      // If the new sequence has a smaller k then the current sequence
      // then insert the new sequence here.
      if (newPtr->k < seqPtr->k)
      {
         if (prevSeqPtr == NULL)
            ip_FirstSequence = newPtr;
         else
            prevSeqPtr->next = newPtr;
         
         newPtr->next = seqPtr;
                     
         break;
      }
      
      prevSeqPtr = seqPtr;
      
      seqPtr = (seq_t *) seqPtr->next;
   } while (seqPtr != NULL);

   UpdateSequenceIndexes();
}

seq_t    *SierpinskiRieselApp::GetSequence(uint64_t k, int64_t c, uint32_t d, seq_t *startSeqPtr) 
{
   seq_t  *seqPtr;

   if (startSeqPtr != NULL)
   {
      seqPtr = startSeqPtr;
      do
      {
         if (seqPtr->k == k && seqPtr->c == c && seqPtr->d == d)
            return seqPtr;

         seqPtr = (seq_t *) seqPtr->next;
      } while (seqPtr != NULL);
   }
   
   seqPtr = ip_FirstSequence;
   do
   {
      if (seqPtr->k == k && seqPtr->c == c && seqPtr->d == d)
         return seqPtr;

      seqPtr = (seq_t *) seqPtr->next;
   } while (seqPtr != startSeqPtr);
   
   FatalError("Sequence for k=%" PRIu64" and c=%+" PRId64" was not found", k, c);
   
   return 0;
}

// Note that this is only called if all workers are paused and waiting for work
void  SierpinskiRieselApp::NotifyAppToRebuild(uint64_t largestPrimeTested)
{
#if defined(USE_OPENCL) || defined(USE_METAL)
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
   if (newSieve || il_MaxK > largestPrimeTested || ib_UseGenericLogic)
   {
      ip_AppHelper = new GenericSequenceHelper(this, largestPrimeTested);
      ib_HaveGenericWorkers = true;
   }
   else
   {
      if (ii_SequenceCount == 1)
         ip_AppHelper = new CisOneWithOneSequenceHelper(this, largestPrimeTested);
      else
         ip_AppHelper = new CisOneWithMultipleSequencesHelper(this, largestPrimeTested);
   }

   if (newSieve)
   {
      uint64_t termsCounted = ip_AppHelper->MakeSubsequencesForNewSieve();
      
      if (termsCounted != il_TermCount)
         FatalError("Expected (%" PRIu64") terms, but only set (%" PRIu64")", il_TermCount, termsCounted);
   }
   else
      ip_AppHelper->MakeSubsequencesForOldSieve(il_TermCount);
   
   ip_AppHelper->LastChanceLogicBeforeSieving();
   
   if (!newSieve && ib_SplitByBestQ)
   {
      if (ib_HaveGenericWorkers)
         WriteOutputTermsFilesByQ();
      else
         WriteToConsole(COT_OTHER, "-S not supported with CisOne sieving logic");
      
      exit(0);
   }
}

void  SierpinskiRieselApp::RemoveSequencesWithNoTerms(void)
{
   seq_t  *seqPtr, *nextSeq, *prevSeq;
   
   prevSeq = seqPtr = ip_FirstSequence;
   do
   {
      bool haveTerm = false;
      
      for (uint32_t n=ii_MinN; n<=ii_MaxN; n++)
      {
         if (seqPtr->nTerms[NBIT(n)])
         {
            haveTerm = true;
            break;
         }
      }
      
      if (!haveTerm)
      {
         seqPtr->nTerms.clear();
         
         nextSeq = (seq_t *) seqPtr->next;
         
         if (seqPtr == ip_FirstSequence)
            ip_FirstSequence = (seq_t *) seqPtr->next;
         else
            prevSeq->next = seqPtr->next;
            
         if (seqPtr->d == 1)
            WriteToConsole(COT_OTHER, "Sequence %" PRIu64"*%u^n%+" PRId64" removed as all terms have a factor", 
               seqPtr->k, ii_Base, seqPtr->c);
         else
            WriteToConsole(COT_OTHER, "Sequence (%" PRIu64"*%u^n%+" PRId64")/%u removed as all terms have a factor", 
               seqPtr->k, ii_Base, seqPtr->c, seqPtr->d);
                  
         xfree(seqPtr);
         
         seqPtr = nextSeq;
         ii_SequenceCount--;
      }
      else
      {
         prevSeq = seqPtr;
         seqPtr = (seq_t *) seqPtr->next;
      }
   } while (seqPtr != NULL);
   
   UpdateSequenceIndexes();
}

void  SierpinskiRieselApp::CheckForLegendreSupport(void)
{
   seq_t  *seqPtr;

   // If we have computed this value, then we don't need to do this again.
   if (ii_SquareFreeB > 0)
      return;

   if (ib_UseGenericLogic)
   {
      ii_SquareFreeB = 1;
      return;
   }

   // If not using generic logic, then we need to compute squareFreeK and squareFreeB even if
   // Legendre tables cannot be used.
   std::vector<uint64_t> primes;
   
   primesieve::generate_primes(sqrt(ii_Base) + 1, &primes);

   // Note that GetSquareFreeFactor() must return a value >= 1
   ii_SquareFreeB = GetSquareFreeFactor(ii_Base, primes);
   
   primesieve::generate_primes(sqrt(il_MaxK) + 1, &primes);
   
   seqPtr = ip_FirstSequence;
   do
   {
      seqPtr->squareFreeK = GetSquareFreeFactor(seqPtr->k, primes);

      if (seqPtr->squareFreeK > PMAX_MAX_62BIT)
      {
         WriteToConsole(COT_OTHER, "Must use generic sieving logic because square free value of %" PRIu64" > 2^62", seqPtr->k);     
         ib_UseGenericLogic = true;
         ii_SquareFreeB = 1;
         break;
      }

      seqPtr->kcCore = -1 * seqPtr->c * seqPtr->squareFreeK;

      seqPtr = (seq_t *) seqPtr->next;
   } while (seqPtr != NULL);
   
   primes.clear();
}

void     SierpinskiRieselApp::ReportFactor(uint64_t theFactor, seq_t *seqPtr, uint32_t n, bool verifyFactor)
{
   uint32_t nbit;
   bool     isPrime = false;
   char     buffer[200];
   bool     needToLock = (theFactor > GetMaxPrimeForSingleWorker());
               
   if (n < ii_MinN || n > ii_MaxN)
      return;

   // We cannot verify these factors, so ignore them
   if (seqPtr->d > 1)
   {
      if (gcd64(seqPtr->d, theFactor) > 1) 
      {      
         uint64_t fPower = theFactor;
         
         while (fPower < seqPtr->d)
            fPower *= theFactor;

         // If d = theFactor^n for some integer n, then we can verify the factor.
         // If not, then there is no way to verify the factor.
         // For example if d = 9 then we can verify if theFactor = 3, but
         // if d = 15, then we cannot verify if theFactor is 3 or 5.
         if (fPower >= seqPtr->d)
           return;
      }
   }

   if (seqPtr->d > 1)
   {
      if (seqPtr->k > 1)
         snprintf(buffer, sizeof(buffer), "(%" PRIu64"*%u^%u%+" PRId64")/%u", seqPtr->k, ii_Base, n, seqPtr->c, seqPtr->d);
      else
         snprintf(buffer, sizeof(buffer), "(%u^%u%+" PRId64")/%u", ii_Base, n, seqPtr->c, seqPtr->d);
   }
   else
      snprintf(buffer, sizeof(buffer), "%" PRIu64"*%u^%u%+" PRId64"", seqPtr->k, ii_Base, n, seqPtr->c);
   
   nbit = NBIT(n);
   
   if (needToLock)
      ip_FactorAppLock->Lock();

   if (seqPtr->nTerms[nbit])
   {
      // Do not remove terms where k*b^n+c is prime.  This means that PRP testing program
      // should identify this term as prime and stop testing other terms of this sequence.
      if (IsPrime(theFactor, seqPtr, n))
      {
         WriteToConsole(COT_OTHER, "%s is prime!", buffer);
         WriteToLog("%s is prime!", buffer);
         isPrime = true;
      }
   }
   
   if (isPrime) {
      if (needToLock)
         ip_FactorAppLock->Release();

      return;
   }
   
   if (seqPtr->nTerms[nbit])
   {
      if (verifyFactor)
         VerifyFactor(theFactor, seqPtr, n);

      il_TermCount--;
      il_FactorCount++;
      seqPtr->nTerms[nbit] = false;
      
      LogFactor(theFactor, "%s", buffer);
   }

   if (needToLock)
      ip_FactorAppLock->Release();
}

void  SierpinskiRieselApp::VerifyFactor(uint64_t theFactor, seq_t *seqPtr, uint32_t n)
{
   MpArith mp(theFactor);
   
   MpRes mmRem = mp.pow(mp.nToRes(ii_Base), n);
   mmRem = mp.mul(mmRem, mp.nToRes(seqPtr->k));

   uint64_t rem = mp.resToN(mmRem);
      
   if (seqPtr->c > 0)
      rem += seqPtr->c;
   else
   {
      int64_t adj = seqPtr->c + theFactor;
      
      while (adj < 0)
         adj += theFactor;
      
      rem += adj;
   }
      
   if (rem >= theFactor)
      rem -= theFactor;

   if (rem == 0 || rem % theFactor == 0)
      return;
         
   char buffer[200];
   
   if (seqPtr->d > 1)
      snprintf(buffer, sizeof(buffer), "Invalid factor: (%" PRIu64"*%u^%u%+" PRId64")/%u mod %" PRIu64" = %" PRIu64"", seqPtr->k, ii_Base, n, seqPtr->c, seqPtr->d, theFactor, rem);
   else
      snprintf(buffer, sizeof(buffer), "Invalid factor: %" PRIu64"*%u^%u%+" PRId64" mod %" PRIu64" = %" PRIu64"", seqPtr->k, ii_Base, n, seqPtr->c, theFactor, rem);
   
   FatalError("%s", buffer);
}

bool  SierpinskiRieselApp::IsPrime(uint64_t p, seq_t *seqPtr, uint32_t n)
{
   __uint128_t bigP;
   
   // Won't try if n is too large
   if (n > 64)
      return false;

   // Just to make sure that our term isn't too big for this logic
   if (log10(ii_Base)*n > 18.0)
      return false;

   bigP = __uint128_t(p);
      
   bigP *= seqPtr->d;
   
   // p = k*b^n+c --> p = k*b^n
   bigP -= seqPtr->c;

   if (bigP % seqPtr->k > 0)
      return false;
   
   // p = k*b^n --> p = b^n
   bigP /= seqPtr->k;

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
uint64_t    SierpinskiRieselApp::GetSquareFreeFactor(uint64_t n, std::vector<uint64_t> primes)
{
   uint64_t c = 1, q, r;
   
   std::vector<uint64_t>::iterator it = primes.begin();
      
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
