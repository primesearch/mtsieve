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
   uint64_t *base;
   uint32_t *power;
   uint32_t  count;
} terms_t;

class MultiFactorialApp : public FactorApp
{
public:
   MultiFactorialApp();

   ~MultiFactorialApp();

   void              Help(void);
   void              AddCommandLineOptions(std::string &shortOpts, struct option *longOpts);
   parse_t           ParseOption(int opt, char *arg, const char *source);
   void              ValidateOptions(void);
   bool              ApplyFactor(uint64_t theFactor, const char *term);
   void              GetExtraTextForSieveStartedMessage(char *extraText, uint32_t maxTextLength);
   
   bool              IsMultiFactorial(void) { return (ii_MultiFactorial > 1); };
   uint32_t          GetMultiFactorial(void) { return ii_MultiFactorial; };
   uint32_t          GetMinN(void) { return ii_MinN; };
   uint32_t          GetMaxN(void) { return ii_MaxN; };

#if defined(USE_OPENCL) || defined(USE_METAL)
   uint32_t          GetMaxGpuSteps(void) { return ii_MaxGpuSteps; };
   uint32_t          GetMaxGpuFactors(void) { return ii_MaxGpuFactors; };
#endif

   bool              ReportFactor(uint64_t theFactor, uint32_t n, int32_t c);

   void              BuildTerms(void);
   terms_t          *GetTerms(void) { return ip_Terms; };
   
protected:
   void              PreSieveHook(void) {};
   bool              PostSieveHook(void) { return true; };
   
   void              NotifyAppToRebuild(uint64_t largestPrimeTested) {};
   
   void              ProcessInputTermsFile(bool haveBitMap);
   bool              IsWritingOutputTermsFile(void){ return true; };
   void              WriteOutputTermsFile(uint64_t largestPrime);
   void              OuptutAdditionalConsoleMessagesUponFinish(void) {};
   
   Worker           *CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested);
   

private:
   void              VerifyFactor(uint64_t theFactor, uint32_t n, int32_t c);
   
   std::vector<bool> iv_PlusTerms;
   std::vector<bool> iv_MinusTerms;

   uint32_t          ii_MultiFactorial;
   uint32_t          ii_MinN;
   uint32_t          ii_MaxN;
   
   terms_t          *ip_Terms;
   
#if defined(USE_OPENCL) || defined(USE_METAL)
   uint32_t          ii_MaxGpuSteps;
   uint32_t          ii_MaxGpuFactors;
#endif
};

#endif

