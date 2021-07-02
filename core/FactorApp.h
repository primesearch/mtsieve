/* FactorApp.h -- (C) Mark Rodenkirch, January 2018

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _FactorApp_H
#define _FactorApp_H

#include <stdio.h>

#include "App.h"
#include "SharedMemoryItem.h"

// As long as we don't expect the factor rate to fall below 1 per day
// then this should be sufficient to capture the rate.
#define MAX_STATUS_COUNT 60 * 5 * 24

typedef struct {
   uint64_t reportTimeUS;
   uint64_t factorsFound;
} report_t;

class FactorApp : public App
{  
public:
   FactorApp(void);
   ~FactorApp(void);
   
protected:
   virtual void      ProcessInputTermsFile(bool haveBitMap) = 0;
   virtual bool      IsWritingOutputTermsFile(void) = 0;
   virtual void      WriteOutputTermsFile(uint64_t largestPrime) = 0;
   virtual bool      ApplyFactor(uint64_t thePrime, const char *term) = 0;
   virtual void      GetExtraTextForSieveStartedMessage(char *extraText) = 0;
   
   void              ParentHelp(void);
   void              ParentAddCommandLineOptions(string &shortOpts, struct option *longOpts);
   parse_t           ParentParseOption(int opt, char *arg, const char *source);
   void              ParentValidateOptions(void);
   
   void              LogStartSievingMessage(void);
   void              Finish(const char *finishMethod, uint64_t elapsedTimeUS, uint64_t largestPrimeTested, uint64_t primesTested);
   void              GetReportStats(char *reportStats, double cpuUtilization);
   bool              StripCRLF(char *line);
   
   void              ResetFactorStats(void);
   
   // Only call this if ip_FactorAppLock has been locked, then release upon return
   void              LogFactor(uint64_t p, const char *fmt, ...) __attribute__ ((format (printf, 3, 4)));
   void              LogFactor(char *factor, const char *fmt, ...) __attribute__ ((format (printf, 3, 4)));
   
   bool              ib_ApplyAndExit;
   
   SharedMemoryItem *ip_FactorAppLock;
   
   // These are only updated by the child class, but any reads/writes of these
   // variables must use ip_FactorAppLock to lock them. 
   uint64_t          il_FactorCount;
   uint64_t          il_PreviousFactorCount;
   uint64_t          il_TermCount;
   
   string            is_InputTermsFileName;
   string            is_InputFactorsFileName;
   string            is_OutputTermsFileName;
   string            is_OutputFactorsFileName;
   
private:
   bool              BuildFactorsPerSecondRateString(uint32_t currentStatusEntry, double cpuUtilization, char *factoringRate);
   bool              BuildSecondsPerFactorRateString(uint32_t currentStatusEntry, double cpuUtilization, char *factoringRate);
   
   FILE             *if_FactorFile;
   time_t            it_CheckpointTime;
   
   // I could use a vector, but I'm lazy
   report_t          ir_ReportStatus[MAX_STATUS_COUNT];
   uint32_t          ii_NextStatusEntry;
};

#endif
