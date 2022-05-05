/* PrimorialApp.h -- (C) Mark Rodenkirch, July 2018

   This class inherits from App.h and has the implementation for this project

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _PrimorialApp_H
#define _PrimorialApp_H

#include <vector>
#include "../core/FactorApp.h"
#include "../core/SharedMemoryItem.h"

// This has to be less than minPrime
#define FIRST_PRIMORIAL       2*3*5
#define FIRST_PRIMORIAL_PRIME 5

class PrimorialApp : public FactorApp
{
public:
   PrimorialApp();

   ~PrimorialApp() {};

   void              Help(void);
   void              AddCommandLineOptions(std::string &shortOpts, struct option *longOpts);
   parse_t           ParseOption(int opt, char *arg, const char *source);
   void              ValidateOptions(void);
   bool              ApplyFactor(uint64_t thePrime, const char *term);
   void              GetExtraTextForSieveStartedMessage(char *extraText);
   
   uint32_t          GetMinPrimorial(void) { return ii_MinPrimorial; };
   uint32_t          GetMaxPrimorial(void) { return ii_MaxPrimorial; };
   
   uint32_t         *GetPrimorialPrimes(uint32_t &numberOfPrimorialPrimes) { numberOfPrimorialPrimes = ii_NumberOfPrimorialPrimes; return ip_PrimorialPrimes; };
   uint16_t         *GetPrimorialPrimeGaps(uint16_t &biggestGap) { biggestGap = ii_BiggestGap; return ip_PrimorialPrimeGaps; };   
   
#if defined(USE_OPENCL) || defined(USE_METAL)
   uint32_t          GetMaxGpuSteps(void) { return ii_MaxGpuSteps; };
   uint32_t          GetMaxGpuFactors(void) { return ii_MaxGpuFactors; };
#endif

   bool              ReportFactor(uint64_t theFactor, uint32_t primorial, int32_t c);

protected:
   void              PreSieveHook(void) {};
   bool              PostSieveHook(void) { return true; };
   
   void              NotifyAppToRebuild(uint64_t largestPrimeTested) {};
   
   void              ProcessInputTermsFile(bool haveBitMap);
   bool              IsWritingOutputTermsFile(void){ return true; };
   void              WriteOutputTermsFile(uint64_t largestPrime);
   
   Worker           *CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested);

   void              VerifyFactor(uint64_t theFactor, uint32_t primorial, int32_t c);

private:
   std::vector<bool> iv_PlusTerms;
   std::vector<bool> iv_MinusTerms;

   uint32_t         *ip_PrimorialPrimes;
   uint32_t          ii_NumberOfPrimorialPrimes;
   
   uint16_t         *ip_PrimorialPrimeGaps;
   uint16_t          ii_BiggestGap;
   
   uint32_t          ii_MinPrimorial;
   uint32_t          ii_MaxPrimorial;
   
#if defined(USE_OPENCL) || defined(USE_METAL)
   uint32_t          ii_MaxGpuSteps;
   uint32_t          ii_MaxGpuFactors;
#endif
};

#endif

