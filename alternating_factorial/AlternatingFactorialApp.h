/* AlternatingFactorialApp.h -- (C) Mark Rodenkirch, July 2017

   This class inherits from App.h and has the implementation for this project

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _AlternatingFactorialApp_H
#define _AlternatingFactorialApp_H

#include <vector>
#include "../core/FactorApp.h"
#include "../core/SharedMemoryItem.h"

class AlternatingFactorialApp : public FactorApp
{
public:
   AlternatingFactorialApp(void);

   ~AlternatingFactorialApp(void) {};

   void              Help(void);
   void              AddCommandLineOptions(std::string &shortOpts, struct option *longOpts);
   parse_t           ParseOption(int opt, char *arg, const char *source);
   void              ValidateOptions(void);
   bool              ApplyFactor(uint64_t theFactor, const char *term);
   void              GetExtraTextForSieveStartedMessage(char *extraText);
   
   uint32_t          GetMinN(void) { return ii_MinN; };
   uint32_t          GetMaxN(void) { return ii_MaxN; };

#if defined(USE_OPENCL) || defined(USE_METAL)
   uint32_t          GetMaxGpuSteps(void) { return ii_MaxGpuSteps; };
   uint32_t          GetMaxGpuFactors(void) { return ii_MaxGpuFactors; };
#endif

   bool              ReportFactor(uint64_t theFactor, uint32_t term);

protected:
   void              PreSieveHook(void) {};
   bool              PostSieveHook(void) { return true; };
   
   void              NotifyAppToRebuild(uint64_t largestPrimeTested) {};
   
   void              ProcessInputTermsFile(bool haveBitMap);
   bool              IsWritingOutputTermsFile(void){ return true; };
   void              WriteOutputTermsFile(uint64_t largestPrime);
   
   Worker           *CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested);
   
private:    
   void              VerifyFactor(uint64_t theFactor, uint32_t term);
   
   std::vector<bool> iv_Terms;
 
   uint32_t          ii_MinN;
   uint32_t          ii_MaxN;

#if defined(USE_OPENCL) || defined(USE_METAL)
   uint32_t          ii_MaxGpuSteps;
   uint32_t          ii_MaxGpuFactors;
#endif
};

#endif

