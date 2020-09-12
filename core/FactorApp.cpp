/* FactorApp.cpp -- (C) Mark Rodenkirch, January 2018

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <inttypes.h>
#include <memory.h>
#include <time.h>
#include <stdarg.h>
#include "Clock.h"
#include "FactorApp.h"

#define CHECKPOINT_SECONDS    3600

FactorApp::FactorApp(void)
{
   ip_FactorAppLock = new SharedMemoryItem("factorapp");
   
   ir_ReportStatus[0].reportTimeUS = Clock::GetCurrentMicrosecond();
   ir_ReportStatus[0].factorsFound = 0;
   ii_NextStatusEntry = 1;

   is_InputTermsFileName = "";
   is_InputFactorsFileName = "";
   is_OutputTermsFileName = "";
   is_OutputFactorsFileName = "";
   
   it_CheckpointTime = time(NULL) + CHECKPOINT_SECONDS;
   il_FactorCount = 0;
   il_TermCount = 0;
   if_FactorFile = 0;
   
   ib_ApplyAndExit = false;
}

FactorApp::~FactorApp(void)
{
   if (if_FactorFile)
      fclose(if_FactorFile);
}

void FactorApp::ParentHelp(void)
{
   App::ParentHelp();

   printf("-A --applyandexit     apply factors and exit (used with -I)\n");
   printf("-i --inputterms=i     input file of remaining candidates\n");
   printf("-I --inputfactors=I   input file with factors (used with -A)\n");
   printf("-o --outputterms=o    output file of remaining candidates\n");
   printf("-O --outputfactors=O  output file with new factors\n");
}

void  FactorApp::ParentAddCommandLineOptions(string &shortOpts, struct option *longOpts)
{
   App::ParentAddCommandLineOptions(shortOpts, longOpts);

   shortOpts += "Ai:o:I:O:";

   AppendLongOpt(longOpts, "applyandexit",   no_argument, 0, 'A');
   AppendLongOpt(longOpts, "inputterms",     required_argument, 0, 'i');
   AppendLongOpt(longOpts, "inputfactors",   required_argument, 0, 'I');
   AppendLongOpt(longOpts, "outputterms",    required_argument, 0, 'o');
   AppendLongOpt(longOpts, "outputfactors",  required_argument, 0, 'O');
}

parse_t FactorApp::ParentParseOption(int opt, char *arg, const char *source)
{
   parse_t status = P_UNSUPPORTED;

   status = App::ParentParseOption(opt, arg, source);
   
   if (status != P_UNSUPPORTED)
      return status;

   switch (opt)
   {
      case 'A':
         ib_ApplyAndExit = true;
         status = P_SUCCESS;
         break;

      case 'i':
         is_InputTermsFileName = arg;
         status = P_SUCCESS;
         break;
         
      case 'I':
         is_InputFactorsFileName = arg;
         status = P_SUCCESS;
         break;

      case 'o':
         is_OutputTermsFileName = arg;
         status = P_SUCCESS;
         break;
         
      case 'O':
         is_OutputFactorsFileName = arg;
         status = P_SUCCESS;
         break;
   }

   return status;
}

void  FactorApp::ParentValidateOptions(void)
{
   char     buffer[1000];
   char    *pos;
   uint32_t factors = 0, applied = 0;
   
   if (is_OutputTermsFileName.length() == 0)
   {
      FatalError("An output terms file name must be specified");
      
      FILE *termsFile = fopen(is_OutputTermsFileName.c_str(), "w");
      
      if (termsFile == NULL)
         FatalError("Could not open terms file %s for output", is_OutputTermsFileName.c_str());
         
      fclose(termsFile);
   }
   
   App::ParentValidateOptions();
   
   if (is_InputFactorsFileName.length() > 0)
   {
      FILE *factorFile = fopen(is_InputFactorsFileName.c_str(), "r");
      
      if (factorFile == NULL)
         FatalError("Could not open factor file %s for input", is_InputFactorsFileName.c_str());
         
      while (fgets(buffer, 1000, factorFile) != NULL)
      {
         if (!StripCRLF(buffer))
            continue;

         // All factors are of the form "p | term"
         pos = strchr(buffer, '|');

         if (pos == 0)
            FatalError("Could not parse factor %s", buffer);
         
         factors++;
         if (ApplyFactor(pos + 2))
            applied++;
      }

      WriteToConsole(COT_OTHER, "Read %u factors from %s.  Removed %u terms", factors, is_InputFactorsFileName.c_str(), applied);
      
      fclose(factorFile);
   }
   
   // I know this is dirty, but it is much easier than other options.
   if (ib_ApplyAndExit)
   {
      WriteOutputTermsFile(il_MinPrime);
      exit(0);
   }
      
   if (is_OutputFactorsFileName.length() > 0)
   {
      if_FactorFile = fopen(is_OutputFactorsFileName.c_str(), "a");
      
      if (if_FactorFile == NULL)
         FatalError("Could not open factor file %s for output", is_OutputFactorsFileName.c_str());
   }
}

void  FactorApp::ResetFactorStats(void)
{
   ir_ReportStatus[0].reportTimeUS = Clock::GetCurrentMicrosecond();
   ir_ReportStatus[0].factorsFound = 0;
   ii_NextStatusEntry = 1;

   il_FactorCount = 0;
}
      
bool  FactorApp::StripCRLF(char *line)
{
   uint32_t n;
   
   n = strlen(line) - 1;
   while (n >= 0 && (line[n] == '\n' || line[n] == '\r'))
   {
      line[n] = 0;
      
      // If the length of the line is 0, return false
      if (n == 0)
         return false;
      
      n--;
   }

   return true;
}

void  FactorApp::LogStartSievingMessage(void)
{
   char  minPrime[30];
   char  maxPrime[30];
   char  startOfMessage[100];
   char  endOfMessage[100];
   char  extraText[200];
   char  fullMessage[500];
   
   *extraText = 0;
   *endOfMessage = 0;
   
   ConvertNumberToShortString(il_MinPrime, (char *) minPrime);
   ConvertNumberToShortString(il_MaxPrime, (char *) maxPrime);
      
   sprintf(startOfMessage, "Sieve started: %s < p < %s with %" PRIu64" terms", minPrime, maxPrime, il_TermCount);

   GetExtraTextForSieveStartedMessage(extraText);
   
   if (il_MaxPrime != il_AppMaxPrime)
   {
      uint64_t minPrime = ((il_MinPrime == 1) ? 2 : il_MinPrime);
      double expectedFactors = ((double) il_TermCount) * (1.0 - log(minPrime) / log(il_MaxPrime));
   
      sprintf(endOfMessage, "expecting %.f factors", expectedFactors);
   }

   if (*extraText != 0)
      if (*endOfMessage != 0)
         sprintf(fullMessage, "%s (%s) (%s)", startOfMessage, extraText, endOfMessage);
      else
         sprintf(fullMessage, "%s (%s)", startOfMessage, extraText);
   else
      if (*endOfMessage != 0)
         sprintf(fullMessage, "%s (%s)", startOfMessage, endOfMessage);
      else
         sprintf(fullMessage, "%s", startOfMessage);
      
   WriteToConsole(COT_OTHER, fullMessage);
           
   WriteToLog(fullMessage);

}

void  FactorApp::Finish(const char *finishMethod, uint64_t elapsedTimeUS, uint64_t largestPrimeTested, uint64_t primesTested)
{
   double   elapsedSeconds = ((double) elapsedTimeUS) / 1000000.0;
   
   if (IsWritingOutputTermsFile())
   {
      WriteOutputTermsFile(largestPrimeTested);
   
      WriteToConsole(COT_OTHER, "%" PRIu64" terms written to %s", il_TermCount, is_OutputTermsFileName.c_str());
   }
   
   WriteToConsole(COT_OTHER, "Primes tested: %" PRIu64".  Factors found: %" PRIu64".  Remaining terms: %" PRIu64".  Time: %.2f seconds.",
           primesTested, il_FactorCount, il_TermCount, elapsedSeconds);
           
   WriteToLog("Sieve %s at p=%" PRIu64".  Primes tested %" PRIu64".  Found %" PRIu64" factors.  %" PRIu64" terms remaining.  Time %.2f seconds\n",
           finishMethod, largestPrimeTested, primesTested, il_FactorCount, il_TermCount, elapsedSeconds);
}

void  FactorApp::GetReportStats(char *reportStats, double cpuUtilization)
{
   char     factoringRate[100];
   uint64_t checkpointPrime;
   uint32_t currnetStatusEntry;

   // Use this as our opportunity to checkpoint current progress
   if (time(NULL) > it_CheckpointTime)
   {
      checkpointPrime = GetLargestPrimeTested(false);
      
      WriteOutputTermsFile(checkpointPrime);
      
      it_CheckpointTime = time(NULL) + CHECKPOINT_SECONDS;
   }
   
   // Lock because workers can update il_FactorCount
   ip_FactorAppLock->Lock();
   
   currnetStatusEntry = ii_NextStatusEntry;
   
   if (ii_NextStatusEntry < MAX_STATUS_COUNT)
      ii_NextStatusEntry++;
   
   if (ii_NextStatusEntry == MAX_STATUS_COUNT) 
   {
      currnetStatusEntry = MAX_STATUS_COUNT - 1;
      
      // Since we have reached the end of the array, move everything one index up in the
      // array.  We don't care if we lose the oldest since we don't expect factor rates
      // to change much if we have been running for that long.
      // Eventually this will be changed to a vector.
      for (uint32_t i=0; i<currnetStatusEntry; i++)
      {   
         ir_ReportStatus[i].reportTimeUS = ir_ReportStatus[i+1].reportTimeUS;
         ir_ReportStatus[i].factorsFound = ir_ReportStatus[i+1].factorsFound;;
      }
   }
   
   ir_ReportStatus[currnetStatusEntry].reportTimeUS = Clock::GetCurrentMicrosecond();
   ir_ReportStatus[currnetStatusEntry].factorsFound = il_FactorCount;
   
   if (il_FactorCount > 0)
   {   
      if (!BuildFactorsPerSecondRateString(currnetStatusEntry, cpuUtilization, factoringRate))
         BuildSecondsPerFactorRateString(currnetStatusEntry, cpuUtilization, factoringRate);
   
      sprintf(reportStats, "%" PRIu64" factors found at %s", il_FactorCount, factoringRate);
   }
   else
      sprintf(reportStats, "no factors found");
      
   
   ip_FactorAppLock->Release();
}

bool  FactorApp::BuildFactorsPerSecondRateString(uint32_t currentStatusEntry, double cpuUtilization, char *factoringRate)
{
   uint32_t     previousStatusEntry;
   uint32_t     factorPrecision;
   const char  *factorRateUnit;
   uint64_t     factorsFound;
   uint64_t     factorTimeUS = 0;
   double       factorsPerSecond = 0;
   double       adjustedFactorTimeSeconds;
   uint64_t     currentReportTimeUS = ir_ReportStatus[currentStatusEntry].reportTimeUS;
   uint64_t     currentFactorsFound = ir_ReportStatus[currentStatusEntry].factorsFound;

   previousStatusEntry = currentStatusEntry;
   
   while (previousStatusEntry > 0)
   {
      previousStatusEntry--;
      
      factorTimeUS = currentReportTimeUS - ir_ReportStatus[previousStatusEntry].reportTimeUS;
      factorsFound = currentFactorsFound - ir_ReportStatus[previousStatusEntry].factorsFound;

      if (factorsFound == 0)
         continue;
      
      // Might add logic here in the future to allow for longer time periods for factor rate
      // calculation but the previous minute should be okay most of the time as most ranges
      // need to reach a removal rate of less than 1 per second which will trigger the 
      // seconds per factor logic.
      break;
   };
   
   // This will adjust based upon the CPU utilization, i.e. number of cores
   // and that the factorTimeUS is in microseconds.
   adjustedFactorTimeSeconds = ((double) factorTimeUS * cpuUtilization) / 1000000.0;
   
   // If finding at least one per second, then compute as factors per second
   if (factorsFound < (uint64_t) adjustedFactorTimeSeconds)
      return false;
   
   // Note that we are computing factors per second
   factorsPerSecond = ((double) factorsFound) / ((double) factorTimeUS);
   
   // Divide the CPU utilization to account for less or more than 1 core
   factorsPerSecond /= cpuUtilization;

   factorRateUnit = "M";
   if (factorsPerSecond < 1.0) factorsPerSecond *= 1000.0, factorRateUnit = "K";
   if (factorsPerSecond < 1.0) factorsPerSecond *= 1000.0, factorRateUnit = "";

   factorPrecision = 0;
   if (factorsPerSecond < 1000.0) factorPrecision = 1;
   if (factorsPerSecond < 100.0)  factorPrecision = 2;
   if (factorsPerSecond < 10.0)   factorPrecision = 3;
   
   sprintf(factoringRate, "%.*f%s f/sec (last %u min)", factorPrecision, factorsPerSecond, factorRateUnit, currentStatusEntry - previousStatusEntry);
   
   return true;
}

bool  FactorApp::BuildSecondsPerFactorRateString(uint32_t currentStatusEntry, double cpuUtilization, char *factoringRate)
{
   uint32_t previousStatusEntry;
   uint64_t factorsFound = 0;
   uint64_t factorTimeUS = 0;
   double   secondsPerFactor;
   uint64_t currentReportTimeUS = ir_ReportStatus[currentStatusEntry].reportTimeUS;
   uint64_t currentFactorsFound = ir_ReportStatus[currentStatusEntry].factorsFound;
   
   previousStatusEntry = currentStatusEntry;

   while (previousStatusEntry > 0)
   {
      previousStatusEntry--;
      
      factorTimeUS = currentReportTimeUS - ir_ReportStatus[previousStatusEntry].reportTimeUS;
      factorsFound = currentFactorsFound - ir_ReportStatus[previousStatusEntry].factorsFound;

      // If no factors found in this time slice, try the next biggest time slice
      if (factorsFound == 0)
         continue;

      // Note that we are computing seconds per factor
      secondsPerFactor = ((double) factorTimeUS) / ((double) factorsFound);
      
      // Convert from ms per factor to sec per factor.
      secondsPerFactor /= 1000000.0;

      // Multiply the CPU utilization to account for less or more than 1 core
      secondsPerFactor *= cpuUtilization;
      
      // The value of 5000 is arbitrary.  The idea is that as the factor rate decreases
      // we want a longer time slice so that natural fluctuations in factor density
      // do not skew the factoring rate too much with successive reports.  For example
      // once the rate is about 10 seconds per factor, it will compute the rate for the
      // past 50000 seconds.  If the rate is about 15 seconds per factor, it will
      // compute the rate based upon the lsat 75000 seconds.
      if ((double) factorsFound > 5000.0 * secondsPerFactor)
         break;
      
      // It is possible that we haven't removed that many factors in the window that
      // is being tracked, therefore this will compute based upon the oldest entry
      // in ir_ReportStatus.  If we want to track for longer then we can increase
      // MAX_STATUS_COUNT, which will be addressed when ir_RepportStatus is changed
      // to a vector.
   };

   if (factorsFound == 0)
   {
      sprintf(factoringRate, "no factors found (last %u min)", currentStatusEntry - previousStatusEntry);
      return true;
   }
   
   // Note that we are computing seconds per factor
   secondsPerFactor = ((double) factorTimeUS) / ((double) factorsFound);
   
   // Convert from ms per factor to sec per factor.
   secondsPerFactor /= 1000000.0;

   // Multiply the CPU utilization to account for less or more than 1 core
   secondsPerFactor *= cpuUtilization;
   
   if (secondsPerFactor > 100)
      sprintf(factoringRate, "%.0f sec per factor (last %u min)", secondsPerFactor, currentStatusEntry - previousStatusEntry);
   else
      sprintf(factoringRate, "%.2f sec per factor (last %u min)", secondsPerFactor, currentStatusEntry - previousStatusEntry);

   return true;      
}

void  FactorApp::LogFactor(uint64_t p, const char *fmt, ...)
{
   if (if_FactorFile == 0)
      return;
      
   fprintf(if_FactorFile, "%" PRIu64" | ", p);
   
   va_list args;

   va_start(args, fmt);
   vfprintf(if_FactorFile, fmt, args);
   va_end(args);
   
   fprintf(if_FactorFile, "\n");
   fflush(if_FactorFile);
}

void  FactorApp::LogFactor(char *factor, const char *fmt, ...)
{
   if (if_FactorFile == 0)
      return;
      
   fprintf(if_FactorFile, "(%s) | ", factor);
   
   va_list args;

   va_start(args, fmt);
   vfprintf(if_FactorFile, fmt, args);
   va_end(args);
   
   fprintf(if_FactorFile, "\n");
   fflush(if_FactorFile);
}
