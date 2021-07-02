/* PrimesInXApp.cpp -- (C) Mark Rodenkirch, October 2012

   PrimesInX/Wall-Sun-Sun Search OpenCL application

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <time.h>
#include "../core/main.h"
#include "../core/Parser.h"
#include "../core/Clock.h"
#include "PrimesInXApp.h"
#include "PrimesInXWorker.h"
#ifdef HAVE_GPU_WORKERS  
#include "PrimesInXGpuWorker.h"
#endif

#define APP_NAME        "pixsieve"
#define APP_VERSION     "2.5"

#define BIT(l)          ((l) - ii_MinLength)

// This is declared in App.h, but implemented here.  This means that App.h
// can remain unchanged if using the CPUSieve framework for other applications.
App *get_app(void)
{
   return new PrimesInXApp();
}

PrimesInXApp::PrimesInXApp() : FactorApp()
{
   SetBanner(APP_NAME " v" APP_VERSION ", a program to find factors of substrings of a decimal string");
   SetLogFileName("pixsieve.log");

   // This is because the assembly code is using SSE to do the mulmods
   SetAppMaxPrime(PMAX_MAX_52BIT);
   
   ii_MinLength = 0;
   ii_MinLengthRemaining = 0;
   ii_MaxLength = 0;

   is_SearchString = "";
   is_FullTerm = "";
   is_StringFileName = "";
   
   ii_e1TermCount = ii_e3TermCount = ii_e6TermCount = ii_e9TermCount = 0;
   ii_e1Terms = ii_e3Terms = ii_e6Terms = ii_e9Terms = 0;

#ifdef HAVE_GPU_WORKERS
   ii_StepL = 5000;
#endif
}

PrimesInXApp::~PrimesInXApp()
{
   if (ii_e1Terms) xfree(ii_e1Terms);
   if (ii_e3Terms) xfree(ii_e3Terms);
   if (ii_e6Terms) xfree(ii_e6Terms);
   if (ii_e9Terms) xfree(ii_e9Terms);
}
   
void PrimesInXApp::Help(void)
{
   FactorApp::ParentHelp();

   printf("-l --minlength=l      minimum length to search\n");
   printf("-L --maxlength=L      maximum length to search\n");
   printf("-s --stringfile=s     file containing a decimal representation of any number\n");
   printf("-S --searchstring=S   starting point of substring to start factoring\n");
   
#ifdef HAVE_GPU_WORKERS
   printf("-N --step=N           N iterated per call to GPU (default %u)\n", ii_StepL);
#endif
}

void  PrimesInXApp::AddCommandLineOptions(string &shortOpts, struct option *longOpts)
{
   FactorApp::ParentAddCommandLineOptions(shortOpts, longOpts);

   shortOpts += "l:L:s:S:i:o:N:";

   AppendLongOpt(longOpts, "minlength",     required_argument, 0, 'l');
   AppendLongOpt(longOpts, "maxlength",     required_argument, 0, 'L');
   AppendLongOpt(longOpts, "stringfile",    required_argument, 0, 's');
   AppendLongOpt(longOpts, "searchstring",  required_argument, 0, 'S');
   AppendLongOpt(longOpts, "inputfile",     required_argument, 0, 'i');
   AppendLongOpt(longOpts, "outputfile",    required_argument, 0, 'o');
   
#ifdef HAVE_GPU_WORKERS
   AppendLongOpt(longOpts, "steps",         required_argument, 0, 'N');
#endif
}


parse_t PrimesInXApp::ParseOption(int opt, char *arg, const char *source)
{
   parse_t status = P_UNSUPPORTED;

   status = FactorApp::ParentParseOption(opt, arg, source);
   if (status != P_UNSUPPORTED) return status;

   switch (opt)
   {
      case 'l':
         status = Parser::Parse(arg, 1, 5000000, ii_MinLength);
         break;

      case 'L':
         status = Parser::Parse(arg, 1, 5000000, ii_MaxLength);
         break;

      case 's':
         is_StringFileName = arg;
         status = P_SUCCESS;
         break;
         
      case 'S':
         is_SearchString = arg;
         status = P_SUCCESS;
         break;
         
#ifdef HAVE_GPU_WORKERS
      case 'N':
         status = Parser::Parse(arg, 1, 1000000000, ii_StepL);
         break;
#endif
   }

   return status;
}

void PrimesInXApp::ValidateOptions(void)
{ 
   if (is_OutputTermsFileName.length() == 0)
      FatalError("An output file name must be specified");

   if (is_InputTermsFileName.length() > 0)
   {
      ProcessInputTermsFile(false);

      il_TermCount = ii_MaxLength - ii_MinLength + 1;
      
      iv_Terms.resize(il_TermCount);
      std::fill(iv_Terms.begin(), iv_Terms.end(), false);
      
      ProcessInputTermsFile(true);
   }
   else
   {
      if (is_StringFileName.length() == 0)
         FatalError("The string file name must be specified");  
  
      ProcessInputStringFile();
   
      if (ii_MinLength > ii_MaxLength)
         FatalError("Min length must be less than max length.");
      
      if (ii_MinLength < is_SearchString.length())
         ii_MinLength = is_SearchString.length();
      
      il_TermCount = ii_MaxLength - ii_MinLength + 1;
      
      iv_Terms.resize(il_TermCount);
      std::fill(iv_Terms.begin(), iv_Terms.end(), true);
   }

   FactorApp::ParentValidateOptions();

   // We want worksize to be divisible by 4
   while (ii_CpuWorkSize % 4 != 0)
      ii_CpuWorkSize++;
}

Worker *PrimesInXApp::CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested)
{
   Worker *theWorker;

#ifdef HAVE_GPU_WORKERS  
   if (gpuWorker)
      theWorker = new PrimesInXGpuWorker(id, this);
   else
#endif
      theWorker = new PrimesInXWorker(id, this);

   return theWorker;
}

void PrimesInXApp::ProcessInputStringFile(void)
{
   FILE  *fPtr;
   char  *buffer, *pos;
   
   // Assume that nobody will ever specify and lmax greater than 100000000
   // because those PRP tests would take many days.
   buffer = (char *) xmalloc(100000000);

   fPtr = fopen(is_StringFileName.c_str(), "r");

   if (fPtr == NULL)
      FatalError("Error reading input file `%s'", is_StringFileName.c_str());
     
   fgets(buffer, 100000000, fPtr);
   fclose(fPtr);

   pos = buffer;
   if (is_SearchString.length() != 0)
   {
      pos = strstr(buffer, is_SearchString.c_str());
      if (pos == 0)
         FatalError("Start string %s not found in input string.", is_SearchString.c_str());
   }

   pos[ii_MaxLength] = 0;
   BuildTerms(pos);
}

void PrimesInXApp::ProcessInputTermsFile(bool haveBitMap)
{
   FILE    *fPtr = fopen(is_InputTermsFileName.c_str(), "r");
   char    *buffer, *pos;
   uint32_t c;
   uint64_t minPrime;

   if (!fPtr)
      FatalError("Unable to open input file %s", is_InputTermsFileName.c_str());

   il_TermCount = 0;
   buffer = (char *) xmalloc(10000100);

   while (fgets(buffer, 10000100, fPtr) != NULL)
   {
      if (!memcmp(buffer, "DECIMAL ", 8))
      {
         pos = strstr(buffer, "//");
         if (pos)
         {
            if (sscanf(pos+13, "%" SCNu64"", &minPrime) == 1)      
               SetMinPrime(minPrime);
         }
 
         if (haveBitMap)
            BuildTerms(buffer+8);
         continue;
      }
      
      if (sscanf(buffer, "%u", &c) != 1)
         FatalError("Line %s is malformed", buffer);

      if (ii_MinLength == 0)
         ii_MinLength = c-1;

      if (c > ii_MaxLength)
         ii_MaxLength = c;
      
      if (haveBitMap)
      {
         iv_Terms[BIT(c)] = true;
         il_TermCount++;
      }
   }

   fclose(fPtr);

   xfree(buffer);
}

bool  PrimesInXApp::ApplyFactor(uint64_t thePrime, const char *term)
{
   uint32_t c;
   
   if (sscanf(term, "pix(%u)", &c) != 1)
      FatalError("Could not parse term %s\n", term);

   if (c < ii_MinLength || c > ii_MaxLength)
      return false;
   
   uint64_t bit = BIT(c);
   
   // No locking is needed because the Workers aren't running yet
   if (iv_Terms[bit])
   {
      iv_Terms[bit] = false;
      il_TermCount--;

      return true;
   }
      
   return false;
}

void  PrimesInXApp::BuildTerms(char *inputTerm)
{
   char  *charTerms, *pos = inputTerm;
   uint32_t index = 0;

   ii_e1Terms = (uint32_t *) xmalloc((size_t) (ii_MaxLength/1+32)*sizeof(uint32_t));
   ii_e3Terms = (uint32_t *) xmalloc((size_t) (ii_MaxLength/3+32)*sizeof(uint32_t));
   ii_e6Terms = (uint32_t *) xmalloc((size_t) (ii_MaxLength/6+32)*sizeof(uint32_t));
   ii_e9Terms = (uint32_t *) xmalloc((size_t) (ii_MaxLength/9+32)*sizeof(uint32_t));
   charTerms  = (char *) xmalloc((size_t) ii_MaxLength+10);

   memset(ii_e3Terms, 0x00, (size_t) (ii_MaxLength/3+32)*sizeof(uint32_t));
   memset(ii_e3Terms, 0x00, (size_t) (ii_MaxLength/6+32)*sizeof(uint32_t));
   memset(ii_e9Terms, 0x00, (size_t) (ii_MaxLength/9+32)*sizeof(uint32_t));

   // Build our array of terms
   while (*pos)
   {
      char x = *pos;

      charTerms[index] = x;

      ii_e1Terms[ii_e1TermCount] = x - '0';
      
      ii_e3Terms[ii_e3TermCount] *= 10;
      ii_e3Terms[ii_e3TermCount] += (x - '0');
      
      ii_e6Terms[ii_e6TermCount] *= 10;
      ii_e6Terms[ii_e6TermCount] += (x - '0');
      
      ii_e9Terms[ii_e9TermCount] *= 10;
      ii_e9Terms[ii_e9TermCount] += (x - '0');
      
      ii_e1TermCount++;

      index++;
      pos++;

      if ((index % 3) == 0) ii_e3TermCount++;
      if ((index % 6) == 0) ii_e6TermCount++;
      if ((index % 9) == 0) ii_e9TermCount++;
         
      if (index == ii_MaxLength)
         break;
   }
  
   if (ii_MinLength > index)
      FatalError("Starting length is longer than the nubmer of digits in the string.");
   
   if (index < ii_MaxLength) {
      ii_MaxLength = index;
      WriteToConsole(COT_OTHER, "Adjusting maxlength to %u, the length of the longest term to factor\n", ii_MaxLength);
   }
   
   charTerms[index] = 0;
   is_FullTerm = charTerms;
   xfree(charTerms);

   // This is a signal to the pixsieve() function that there are no more terms
   ii_e1Terms[ii_e1TermCount] = 10;
   ii_e3Terms[ii_e3TermCount] = 1000;
   ii_e6Terms[ii_e6TermCount] = 1000000;
   ii_e9Terms[ii_e9TermCount] = 1000000000;
}

uint32_t	*PrimesInXApp::Get3DigitTermsCopy(void)
{
   uint32_t *ptr;
   
   ptr = (uint32_t *) xmalloc((size_t) (ii_MaxLength/3+32)*sizeof(uint32_t));
   memcpy(ptr, ii_e3Terms, (size_t) (ii_MaxLength/3+32)*sizeof(uint32_t));
   
   return ptr;
}

uint32_t	*PrimesInXApp::Get6DigitTermsCopy(void)
{
   uint32_t *ptr;
   
   ptr = (uint32_t *) xmalloc((size_t) (ii_MaxLength/6+32)*sizeof(uint32_t));
   memcpy(ptr, ii_e6Terms, (size_t) (ii_MaxLength/6+32)*sizeof(uint32_t));
   
   return ptr;
}

uint32_t	*PrimesInXApp::Get9DigitTermsCopy(void)
{
   uint32_t *ptr;
   
   ptr = (uint32_t *) xmalloc((size_t) (ii_MaxLength/9+32)*sizeof(uint32_t));
   memcpy(ptr, ii_e9Terms, (size_t) (ii_MaxLength/9+32)*sizeof(uint32_t));
   
   return ptr;
}

void PrimesInXApp::WriteOutputTermsFile(uint64_t largestPrime)
{
   FILE    *termsFile = fopen(is_OutputTermsFileName.c_str(), "w");
   uint64_t termsCounted = 0;
   uint32_t maxLength = 0;
   char     termCountStr[50];
   char     termsCountedStr[50];

   if (!termsFile)
      FatalError("Unable to open input file %s", is_OutputTermsFileName.c_str());
      
   ip_FactorAppLock->Lock();
   
   for (uint32_t l=ii_MinLength; l<=ii_MaxLength; l++)
   {
      if (iv_Terms[BIT(l)])
         maxLength = l;
   }
   
   // Set the terminator so that we don't write unnecessary data
   is_FullTerm[maxLength] = 0;

   fprintf(termsFile, "DECIMAL %s // Sieved to %" PRIu64"\n", is_FullTerm.c_str(), largestPrime);
      
   for (uint32_t l=ii_MinLength; l<=ii_MaxLength; l++)
   {
      if (iv_Terms[BIT(l)])
      {
         fprintf(termsFile, "%u\n", l);
         termsCounted++;
      }
   }

   fclose(termsFile);

   if (termsCounted != il_TermCount)
   {
      sprintf(termCountStr, "%" PRIu64"", il_TermCount);
      sprintf(termsCountedStr, "%" PRIu64"", termsCounted);
      
      FatalError("Something is wrong.  Counted terms (%s) != expected terms (%s)", termsCountedStr, termCountStr);
   }
      
   ip_FactorAppLock->Release();
}

void PrimesInXApp::GetExtraTextForSieveStartedMessage(char *extraText)
{                 
   if (is_SearchString.length() > 0)
      sprintf(extraText, "%d <= length <= %d, terms starting with %s", ii_MinLength, ii_MaxLength, is_SearchString.c_str());
   else
      sprintf(extraText, "%d <= length <= %d", ii_MinLength, ii_MaxLength);
}

bool PrimesInXApp::ReportFactor(uint64_t p, uint32_t n)
{
   bool newFactor = false;
   
   if (n < ii_MinLength || n > ii_MaxLength)
      return false;
   
   ip_FactorAppLock->Lock();
      
   uint64_t bit = BIT(n);

   if (iv_Terms[bit])
   {
      newFactor = true;
      iv_Terms[bit] = false;
      il_FactorCount++;
      il_TermCount--;
      
      LogFactor(p, "pix(%u)", n);
   }
   
   ip_FactorAppLock->Release();
   
   return newFactor;
}

void PrimesInXApp::ReportPrime(uint64_t p, uint32_t n)
{
   char  pStr[50];
   
   if (n < ii_MinLength || n > ii_MaxLength)
      return;
   
   ip_FactorAppLock->Lock();
      
   uint64_t bit = BIT(n);

   if (iv_Terms[bit])
   {
      iv_Terms[bit] = false;
      il_FactorCount++;
      il_TermCount--;

      sprintf(pStr, "%" PRIu64"", p);
 
      LogFactor(p, "pix(%u)", n);
         
      WriteToConsole(COT_OTHER, "pix(%u) is prime! (%s)", n, pStr);

      WriteToLog("pix(%u) is prime! (%s)", n, pStr);
   }
   
   ip_FactorAppLock->Release();
}
