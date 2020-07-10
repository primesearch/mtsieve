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

   // Use this as our opportunity to checkpoint current progress
   if (time(NULL) > it_CheckpointTime)
   {
      checkpointPrime = GetLargestPrimeTested(false);
      
      WriteOutputTermsFile(checkpointPrime);
      
      it_CheckpointTime = time(NULL) + CHECKPOINT_SECONDS;
   }
   
   // Lock because workers can update il_FactorCount
   ip_FactorAppLock->Lock();
   
   if (il_FactorCount > 0)
   {   
      GetFactoringRate(il_FactorCount, cpuUtilization, factoringRate);
   
      sprintf(reportStats, "%" PRIu64" factors found at %s", il_FactorCount, factoringRate);
   }
   else
      sprintf(reportStats, "no factors found");
      
   ip_FactorAppLock->Release();
}

void  FactorApp::GetFactoringRate(uint64_t factorCount, double cpuUtilization, char *factoringRate)
{
   int statusEntry = ii_NextStatusEntry;
   
   if (ii_NextStatusEntry == STATUS_COUNT) 
   {
      statusEntry = STATUS_COUNT - 1;
      
      // Since we have reached the end of the array, move everything one index up in the
      // array.  We don't care if we lose the oldest since we don't expect factor rates
      // to change much if we have been running for that long.
      for (int i=0; i<statusEntry; i++)
      {   
         ir_ReportStatus[i].reportTimeUS = ir_ReportStatus[i+1].reportTimeUS;
         ir_ReportStatus[i].factorsFound = ir_ReportStatus[i+1].factorsFound;;
      }
   }
   
   ir_ReportStatus[statusEntry].reportTimeUS = Clock::GetCurrentMicrosecond();
   ir_ReportStatus[statusEntry].factorsFound = factorCount;

   BuildFactoringRateString(cpuUtilization, factoringRate);
   
   if (ii_NextStatusEntry < STATUS_COUNT)
      ii_NextStatusEntry++;;
}

void  FactorApp::BuildFactoringRateString(double cpuUtilization, char *factoringRate)
{
   uint32_t previousStatusEntry;
   double   factorRate;
   uint32_t factorPrecision;
   uint32_t statusEntry = ii_NextStatusEntry;
   const char  *factorRateUnit;
   uint64_t factorsFound;
   uint64_t factorTimeUS, adjustedFactorTimeUS;
   uint64_t currentReportTimeUS = ir_ReportStatus[statusEntry].reportTimeUS;
   uint64_t currentFactorsFound = ir_ReportStatus[statusEntry].factorsFound;

   if (statusEntry == 1)
   {
      factorTimeUS = currentReportTimeUS - ir_ReportStatus[statusEntry-1].reportTimeUS;
      factorsFound = currentFactorsFound - ir_ReportStatus[statusEntry-1].factorsFound;
   }
   else
   {
      previousStatusEntry = ((statusEntry > 5) ? (statusEntry - 5) : statusEntry - 1);
      
      do {
         factorTimeUS = currentReportTimeUS - ir_ReportStatus[previousStatusEntry-1].reportTimeUS;
         factorsFound = currentFactorsFound - ir_ReportStatus[previousStatusEntry-1].factorsFound;

         if (previousStatusEntry == 0)
            break;
         
         // Add a minute until we have a range of time with at least one factor found
         previousStatusEntry--;
      } while (factorsFound == 0);
  
      // This will adjust based upon the CPU utilization, i.e. number of cores
      adjustedFactorTimeUS = (uint64_t) ((double) factorTimeUS * cpuUtilization);

      // If finding at least one per second in the previous 5 minutes, then compute as factors per second
      if (factorsFound > adjustedFactorTimeUS / 1000000)
      {
         // Note that we are computing factors per second
         factorRate = ((double) factorsFound) / ((double) factorTimeUS);
         
         // Divide the CPU utilization to account for less or more than 1 core
         factorRate /= cpuUtilization;

         factorRateUnit = "M";
         if (factorRate < 1.0) factorRate *= 1000.0, factorRateUnit = "K";
         if (factorRate < 1.0) factorRate *= 1000.0, factorRateUnit = "";

         factorPrecision = 0;
         if (factorRate < 1000.0) factorPrecision = 1;
         if (factorRate < 100.0)  factorPrecision = 2;
         if (factorRate < 10.0)   factorPrecision = 3;
         
         sprintf(factoringRate, "%.*f%s f/sec", factorPrecision, factorRate, factorRateUnit);
         return;
      }

      // Use the status entry from for the past hour
      previousStatusEntry = ((statusEntry > 60) ? (statusEntry - 60) : statusEntry - 1);

      do {
         factorTimeUS = currentReportTimeUS - ir_ReportStatus[previousStatusEntry].reportTimeUS;
         factorsFound = currentFactorsFound - ir_ReportStatus[previousStatusEntry].factorsFound;

         if (previousStatusEntry == 0)
            break;

         // Add a minute until we have a range of time with at least one factor found
         previousStatusEntry--;
      } while (factorsFound == 0);
         
      // This will adjust based upon the CPU utilization, i.e. number of cores
      adjustedFactorTimeUS = (uint64_t) ((double) factorTimeUS * cpuUtilization);

      // If finding less than 1 per minute for the past hour, use the past day
      if (adjustedFactorTimeUS < (statusEntry - previousStatusEntry))
      {
         // We will ignore the first hour as we want to avoid inflated numbers due to the
         // number of terms removed in the first hour.  This is arbitrary, but we're assuming
         // that most users will be running for at least 30 minutes before they stop sieving.
         // If we haven't been running for at least 30 minutes, compute the rate based upon
         // when we started sieving.
         previousStatusEntry = ((statusEntry > 30) ? 30 : 0);

         // If running for at least 24 hours, then use the last 23.5 hours to compute the rate.
         // If we use the last 24 hours, then we might get rates for the first hour of a new
         // sieve which will inflate the factoring rate for an hour.
         if (statusEntry > 1440)
            previousStatusEntry = statusEntry - 1410;
      
         factorTimeUS = currentReportTimeUS - ir_ReportStatus[previousStatusEntry].reportTimeUS;
         factorsFound = currentFactorsFound - ir_ReportStatus[previousStatusEntry].factorsFound;
         
         // If fewer than 24 factors found, then use the oldest status entry to compute the rate
         // as that is our best chance of getting a meaningful factor rate.
         if (factorsFound < 24)
         {
            factorTimeUS = currentReportTimeUS - ir_ReportStatus[0].reportTimeUS;
            factorsFound = currentFactorsFound - ir_ReportStatus[0].factorsFound;
         }
      }
   }

   if (factorsFound > 0)
   {
      // Note that we are computing seconds per factor
      factorRate = ((double) factorTimeUS) / ((double) factorsFound);
      
      // Convert from ms per factor to sec per factor.
      factorRate /= 1000000.0;

      // Multiply the CPU utilization to account for less or more than 1 core
      factorRate *= cpuUtilization;
      
      if (factorRate > 100)
         sprintf(factoringRate, "%.0f sec per factor", factorRate);
      else
         sprintf(factoringRate, "%.2f sec per factor", factorRate);
   }
   else
         sprintf(factoringRate, "no factors found");
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
