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
   
   is_InputTermsFileName = "";
   is_InputFactorsFileName = "";
   is_OutputTermsFileName = "";
   is_OutputFactorsFileName = "";
   
   it_CheckpointTime = time(NULL) + CHECKPOINT_SECONDS;
   il_FactorCount = 0;
   il_PreviousFactorCount = 0;
   il_TermCount = 0;
   if_FactorFile = 0;
   
   id_FPSTarget = 0.0;
   id_SPFTarget = 0.0;
   
   ii_MinutesForStatus = 0;
   
   ib_ApplyAndExit = false;
   
   ResetFactorStats();
}

FactorApp::~FactorApp(void)
{
   if (if_FactorFile)
      fclose(if_FactorFile);

   delete ip_FactorAppLock;
}

void FactorApp::ParentHelp(void)
{
   App::ParentHelp();

   printf("-A --applyandexit     apply factors and exit (used with -I)\n");
   printf("-i --inputterms=i     input file of remaining candidates\n");
   printf("-I --inputfactors=I   input file with factors (used with -A)\n");
   printf("-o --outputterms=o    output file of remaining candidates\n");
   printf("-O --outputfactors=O  output file with new factors\n");
   printf("-4 --fpstarget=4      stop sieving ASAP when factors per second falls below this value\n");
   printf("-5 --spftarget=5      stop sieving ASAP when seconds per factor goes above this value\n");
   printf("-6 --minutesforspf=6  maximum number of minutes to use for computed seconds per factor\n");
}

void  FactorApp::ParentAddCommandLineOptions(std::string &shortOpts, struct option *longOpts)
{
   App::ParentAddCommandLineOptions(shortOpts, longOpts);

   shortOpts += "Ai:o:I:O:4:5:6:";

   AppendLongOpt(longOpts, "applyandexit",   no_argument, 0, 'A');
   AppendLongOpt(longOpts, "inputterms",     required_argument, 0, 'i');
   AppendLongOpt(longOpts, "inputfactors",   required_argument, 0, 'I');
   AppendLongOpt(longOpts, "outputterms",    required_argument, 0, 'o');
   AppendLongOpt(longOpts, "outputfactors",  required_argument, 0, 'O');
   AppendLongOpt(longOpts, "fpstarget",      required_argument, 0, '4');
   AppendLongOpt(longOpts, "spftarget",      required_argument, 0, '5');
   AppendLongOpt(longOpts, "minutesforspf",  required_argument, 0, '6');
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
         
      case '4':
         status = Parser::Parse(arg, 0.1, 1000000000.0, id_FPSTarget);
         break;
         
      case '5':
         status = Parser::Parse(arg, 0.1, 1000000000.0, id_SPFTarget);
         break;

      case '6':
         status = Parser::Parse(arg, 1, MAX_FACTOR_REPORT_COUNT, ii_MinutesForStatus);
         break;
   }

   return status;
}

void  FactorApp::ParentValidateOptions(void)
{
   char     buffer[1000];
   char    *pos;
   uint32_t factors = 0, applied = 0;
   uint64_t thePrime;

   if (id_FPSTarget > 0.0 && id_SPFTarget > 0.0)
      FatalError("Cannot specify both -4 and -5");
   
   if (ii_MinutesForStatus > 0 && id_FPSTarget > 0.0)
      WriteToConsole(COT_OTHER, "Option -6 ignored as it only applies to seconds per factor calculations");
      
   if (ii_MinutesForStatus == 0)
      ii_MinutesForStatus = MAX_FACTOR_REPORT_COUNT;
      
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
   
      while (fgets(buffer, sizeof(buffer), factorFile) != NULL)
      {
         if (!StripCRLF(buffer))
            continue;

         if (sscanf(buffer, "%" SCNu64"", &thePrime) != 1)
            FatalError("Could not parse prime from string %s", buffer);
            
         // All factors are of the form "p | term"
         pos = strchr(buffer, '|');

         if (pos == 0)
            FatalError("Could not extract candidate from %s", buffer);

         *pos = 0;
         
         factors++;
         if (ApplyFactor(thePrime, pos + 2))
            applied++;
      }
      
      fclose(factorFile);
      
      snprintf(buffer, sizeof(buffer), "Read %u factors from %s which removed %u terms.", factors, is_InputFactorsFileName.c_str(), applied);
      
      WriteToConsole(COT_OTHER, "%s", buffer);

      WriteToLog("%s", buffer);
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

   il_PreviousFactorCount += il_FactorCount;
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
      
   snprintf(startOfMessage, sizeof(startOfMessage), "Sieve started: %s <= p <= %s with %" PRIu64" terms", minPrime, maxPrime, il_TermCount);

   GetExtraTextForSieveStartedMessage(extraText, 200);
   
   if (il_MaxPrime != il_AppMaxPrime)
   {
      uint64_t minPrime = ((il_MinPrime == 1) ? 2 : il_MinPrime);
      double expectedFactors = ((double) il_TermCount) * (1.0 - log(minPrime) / log(il_MaxPrime));
   
      snprintf(endOfMessage, sizeof(endOfMessage), "expecting %.f factors", expectedFactors);
   }

   if (*extraText != 0)
      if (*endOfMessage != 0)
         snprintf(fullMessage, sizeof(fullMessage), "%s (%s) (%s)", startOfMessage, extraText, endOfMessage);
      else
         snprintf(fullMessage, sizeof(fullMessage), "%s (%s)", startOfMessage, extraText);
   else
      if (*endOfMessage != 0)
         snprintf(fullMessage, sizeof(fullMessage), "%s (%s)", startOfMessage, endOfMessage);
      else
         snprintf(fullMessage, sizeof(fullMessage), "%s", startOfMessage);
      
   WriteToConsole(COT_OTHER, "%s", fullMessage);

   WriteToLog("%s", fullMessage);
}

void  FactorApp::Finish(const char *finishMethod, uint64_t elapsedTimeUS, uint64_t largestPrimeTested, uint64_t primesTested)
{
   double   elapsedSeconds = ((double) elapsedTimeUS) / 1000000.0;
   uint64_t factorCount = il_FactorCount + il_PreviousFactorCount;

   if (IsWritingOutputTermsFile())
   {
      WriteOutputTermsFile(largestPrimeTested);
   
      WriteToConsole(COT_OTHER, "%" PRIu64" terms written to %s", il_TermCount, is_OutputTermsFileName.c_str());
   }
   
   WriteToConsole(COT_OTHER, "Primes tested: %" PRIu64".  Factors found: %" PRIu64".  Remaining terms: %" PRIu64".  Time: %.2f seconds.",
           primesTested, factorCount, il_TermCount, elapsedSeconds);
           
   WriteToLog("Sieve %s at p=%" PRIu64".  Primes tested %" PRIu64".  Found %" PRIu64" factors.  %" PRIu64" terms remaining.  Time %.2f seconds\n",
           finishMethod, largestPrimeTested, primesTested, factorCount, il_TermCount, elapsedSeconds);
}

void  FactorApp::GetReportStats(char *reportStats, uint32_t maxStatsLength, double cpuUtilization)
{
   char     factoringRate[100];
   uint64_t checkpointPrime;
   uint32_t currentStatusEntry;

   // Use this as our opportunity to checkpoint current progress
   if (time(NULL) > it_CheckpointTime)
   {
      checkpointPrime = GetLargestPrimeTested(false);
      
      WriteOutputTermsFile(checkpointPrime);
      
      it_CheckpointTime = time(NULL) + CHECKPOINT_SECONDS;
   }
   
   // Lock because workers can update il_FactorCount
   ip_FactorAppLock->Lock();
   
   currentStatusEntry = ii_NextStatusEntry;
   
   if (ii_NextStatusEntry < MAX_FACTOR_REPORT_COUNT)
      ii_NextStatusEntry++;
   
   if (ii_NextStatusEntry == MAX_FACTOR_REPORT_COUNT) 
   {
      currentStatusEntry = MAX_FACTOR_REPORT_COUNT - 1;
      
      // Since we have reached the end of the array, move everything one index up in the
      // array.  We don't care if we lose the oldest since we don't expect factor rates
      // to change much if we have been running for that long.
      // Eventually this will be changed to a vector.
      for (uint32_t i=0; i<currentStatusEntry; i++)
      {   
         ir_ReportStatus[i].reportTimeUS = ir_ReportStatus[i+1].reportTimeUS;
         ir_ReportStatus[i].factorsFound = ir_ReportStatus[i+1].factorsFound;
      }
   }
   
   ir_ReportStatus[currentStatusEntry].reportTimeUS = Clock::GetCurrentMicrosecond();
   ir_ReportStatus[currentStatusEntry].factorsFound = il_FactorCount;
   
   if (il_FactorCount > 0)
   {   
      if (!BuildFactorsPerSecondRateString(currentStatusEntry, cpuUtilization, factoringRate))
         BuildSecondsPerFactorRateString(currentStatusEntry, cpuUtilization, factoringRate);
   
      snprintf(reportStats, maxStatsLength, "%" PRIu64" factors found at %s", il_FactorCount + il_PreviousFactorCount, factoringRate);
   }
   else
      snprintf(reportStats, maxStatsLength, "no factors found");
      
   
   ip_FactorAppLock->Release();
}

bool  FactorApp::BuildFactorsPerSecondRateString(uint32_t currentStatusEntry, double cpuUtilization, char *factoringRate)
{
   uint32_t     previousStatusEntry;
   uint32_t     factorPrecision;
   const char  *factorRateUnit;
   uint64_t     factorsFound = 0;
   uint64_t     factorTimeUS = 0;
   double       factorsPerUS = 0;
   double       adjustedFactorTimeSeconds;
   uint64_t     currentReportTimeUS = ir_ReportStatus[currentStatusEntry].reportTimeUS;
   uint64_t     currentFactorsFound = ir_ReportStatus[currentStatusEntry].factorsFound;

   previousStatusEntry = currentStatusEntry;

#if defined(USE_OPENCL) || defined(USE_METAL)
   // In the GPU we could have an uneven distribution of factors because factors are only
   // reported when the GPU is done with its chunk.  If each chunk requires more than a
   // few seconds in the GPU, the factor rate will bounce up and down making it harder to
   // compute the removal rate.  To reduce the size of the "bounce", we will use abort
   // larger time slice.  There will still be bouncing, but it will be less pronounced
   // as the runtime increases.
   if (previousStatusEntry > 60)
      previousStatusEntry -= 60;
   else
      previousStatusEntry = 1;
#endif

   while (previousStatusEntry > 0)
   {
      previousStatusEntry--;
      
      factorTimeUS = currentReportTimeUS - ir_ReportStatus[previousStatusEntry].reportTimeUS;
      factorsFound = currentFactorsFound - ir_ReportStatus[previousStatusEntry].factorsFound;

      // If no factors found in this time slice, use a larger time slice
      if (factorsFound == 0) {
         if (previousStatusEntry > 60)
            previousStatusEntry -= 60;
         else
            previousStatusEntry = 1;
         continue;
      }
      
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
   factorsPerUS = ((double) factorsFound) / ((double) factorTimeUS);
   
   // Divide the CPU utilization to account for less or more than 1 core
   factorsPerUS /= cpuUtilization;

   if (id_FPSTarget > factorsPerUS * 1000000.0)
      Interrupt("f/sec is slower than target");
   
   factorRateUnit = "M";
   if (factorsPerUS < 1.0) factorsPerUS *= 1000.0, factorRateUnit = "K";
   if (factorsPerUS < 1.0) factorsPerUS *= 1000.0, factorRateUnit = "";

   factorPrecision = 0;
   if (factorsPerUS < 1000.0) factorPrecision = 1;
   if (factorsPerUS < 100.0)  factorPrecision = 2;
   if (factorsPerUS < 10.0)   factorPrecision = 3;
   
   snprintf(factoringRate, 100, "%.*f%s f/sec (last %u min)", factorPrecision, factorsPerUS, factorRateUnit, currentStatusEntry - previousStatusEntry);
   
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

#if defined(USE_OPENCL) || defined(USE_METAL)
   // In the GPU we could have an uneven distribution of factors because factors are only
   // reported when the GPU is done with its chunk.  If each chunk requires more than a
   // few seconds in the GPU, the factor rate will bounce up and down making it harder to
   // compute the removal rate.  To reduce the size of the "bounce", we will use abort
   // larger time slice.  There will still be bouncing, but it will be less pronounced
   // as the runtime increases.
   if (previousStatusEntry > 60)
      previousStatusEntry -= 60;
   else
      previousStatusEntry = 1;
#endif

   while (previousStatusEntry > 0)
   {
      previousStatusEntry--;
      
      factorTimeUS = currentReportTimeUS - ir_ReportStatus[previousStatusEntry].reportTimeUS;
      factorsFound = currentFactorsFound - ir_ReportStatus[previousStatusEntry].factorsFound;

      // If no factors found in this time slice, go back another minute
      if (factorsFound == 0) {
         previousStatusEntry = 1;
         continue;
      }

      // Note that we are computing seconds per factor
      secondsPerFactor = ((double) factorTimeUS) / ((double) factorsFound);
      
      // Convert from ms per factor to sec per factor.
      secondsPerFactor /= 1000000.0;

      // Multiply the CPU utilization to account for less or more than 1 core
      secondsPerFactor *= cpuUtilization;
      
      if (currentStatusEntry - previousStatusEntry >= ii_MinutesForStatus)
         break;
      
      // It is possible that we haven't removed that many factors in the window that
      // is being tracked, therefore this will compute based upon the oldest entry
      // in ir_ReportStatus.  If we want to track for longer then we can increase
      // MAX_FACTOR_REPORT_COUNT, which will be addressed when ir_RepportStatus is changed
      // to a vector.
   };
   
   if (factorsFound == 0)
   {
      snprintf(factoringRate, 100, "no factors found (last %u min)", currentStatusEntry - previousStatusEntry);
      return true;
   }
   
   // Note that we are computing seconds per factor
   secondsPerFactor = ((double) factorTimeUS) / ((double) factorsFound);
   
   // Convert from ms per factor to sec per factor.
   secondsPerFactor /= 1000000.0;

   // Multiply the CPU utilization to account for less or more than 1 core
   secondsPerFactor *= cpuUtilization;
   
   if (id_FPSTarget > 0.0)
      Interrupt("sec/f is slower than target");
   
   if (id_SPFTarget > 0.0 && id_SPFTarget < secondsPerFactor)
      Interrupt("sec/f is slower than target");
   
   if (secondsPerFactor > 100)
      snprintf(factoringRate, 100, "%.0f sec per factor (last %u min)", secondsPerFactor, currentStatusEntry - previousStatusEntry);
   else
      snprintf(factoringRate, 100, "%.2f sec per factor (last %u min)", secondsPerFactor, currentStatusEntry - previousStatusEntry);
   
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
