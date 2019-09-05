/* MultiFactorialApp.h -- (C) Mark Rodenkirch, October 2012

   This class inherits from App.h and has the implementation for this project

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _MultiFactorialApp_H
#define _MultiFactorialApp_H

#include <vector>
#include "../core/FactorApp.h"
#include "../core/SharedMemoryItem.h"

class MultiFactorialApp : public FactorApp
{
public:
   MultiFactorialApp();

   ~MultiFactorialApp() {};

   void              Help(void);
   void              AddCommandLineOptions(string &shortOpts, struct option *longOpts);
   parse_t           ParseOption(int opt, char *arg, const char *source);
   void              ValidateOptions(void);
   bool              ApplyFactor(const char *term);
   void              GetExtraTextForSieveStartedMessage(char *extraText);
   
   bool              IsMultiFactorial(void) { return (ii_MultiFactorial > 1); };
   uint32_t          GetMultiFactorial(void) { return ii_MultiFactorial; };
   uint32_t          GetMinN(void) { return ii_MinN; };
   uint32_t          GetMaxN(void) { return ii_MaxN; };

#ifdef HAVE_GPU_WORKERS
   uint32_t          GetStepN(void) { return ii_StepN; };
#endif

   bool              ReportFactor(uint64_t p, uint32_t n, int32_t c);
   void              ReportPrime(uint64_t p, uint32_t n, int32_t c);

protected:
   void              PreSieveHook(void) {};
   bool              PostSieveHook(void) { return true; };
   
   void              NotifyAppToRebuild(void) {};
   
   void              ProcessInputTermsFile(bool haveBitMap);
   bool              IsWritingOutputTermsFile(void){ return true; };
   void              WriteOutputTermsFile(uint64_t largestPrime);
   
   Worker           *CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested);

private:
   vector<bool>      iv_PlusTerms;
   vector<bool>      iv_MinusTerms;

   uint32_t          ii_MultiFactorial;
   uint32_t          ii_MinN;
   uint32_t          ii_MaxN;
   
#ifdef HAVE_GPU_WORKERS
   uint32_t          ii_StepN;
#endif
};

#endif

