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

typedef struct {
   uint32_t  mf;
   uint64_t  maxNForTerm;
   uint64_t  termCount;
   uint64_t *termList;
} terms_t;

class MultiFactorialApp : public FactorApp
{
public:
   MultiFactorialApp();

   ~MultiFactorialApp() {};

   void              Help(void);
   void              AddCommandLineOptions(string &shortOpts, struct option *longOpts);
   parse_t           ParseOption(int opt, char *arg, const char *source);
   void              ValidateOptions(void);
   bool              ApplyFactor(uint64_t thePrime, const char *term);
   void              GetExtraTextForSieveStartedMessage(char *extraText);
   
   bool              IsMultiFactorial(void) { return (ii_MultiFactorial > 1); };
   uint32_t          GetMultiFactorial(void) { return ii_MultiFactorial; };
   uint32_t          GetMinN(void) { return ii_MinN; };
   uint32_t          GetMaxN(void) { return ii_MaxN; };

#ifdef HAVE_GPU_WORKERS
   uint32_t          GetMaxGpuSteps(void) { return ii_MaxGpuSteps; };
   uint32_t          GetMaxGpuFactors(void) { return ii_MaxGpuFactors; };
#endif

   bool              ReportFactor(uint64_t p, uint32_t n, int32_t c);

   terms_t          *GetTerms(void);
   
protected:
   void              PreSieveHook(void) {};
   bool              PostSieveHook(void) { return true; };
   
   void              NotifyAppToRebuild(uint64_t largestPrimeTested) {};
   
   void              ProcessInputTermsFile(bool haveBitMap);
   bool              IsWritingOutputTermsFile(void){ return true; };
   void              WriteOutputTermsFile(uint64_t largestPrime);
   
   Worker           *CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested);
   

private:
   bool              VerifyFactor(bool badFactorIsFatal, uint64_t p, uint32_t n, int32_t c);
   
   vector<bool>      iv_PlusTerms;
   vector<bool>      iv_MinusTerms;

   uint32_t          ii_MultiFactorial;
   uint32_t          ii_MinN;
   uint32_t          ii_MaxN;
   
#ifdef HAVE_GPU_WORKERS
   uint32_t          ii_MaxGpuSteps;
   uint32_t          ii_MaxGpuFactors;
#endif
};

#endif

