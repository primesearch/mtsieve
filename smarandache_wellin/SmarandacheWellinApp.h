/* SmarandacheWellinApp.h -- (C) Mark Rodenkirch, March 2025

   This class inherits from App.h and has the implementation for this project

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _SmarandacheWellinApp_H
#define _SmarandacheWellinApp_H

#include <vector>
#include "../core/FactorApp.h"
#include "../core/SharedMemoryItem.h"

typedef struct {
   uint32_t  termCount;
   uint32_t *termList;
} terms_t;

class SmarandacheWellinApp : public FactorApp
{
public:
   SmarandacheWellinApp();

   ~SmarandacheWellinApp() {};

   void              Help(void);
   void              AddCommandLineOptions(std::string &shortOpts, struct option *longOpts);
   parse_t           ParseOption(int opt, char *arg, const char *source);
   void              ValidateOptions(void);
   bool              ApplyFactor(uint64_t theFactor, const char *term);
   void              GetExtraTextForSieveStartedMessage(char *extraText, uint32_t maxTextLength);

   uint32_t          GetMinN(void) { return ii_MinN; };
   uint32_t          GetMaxN(void) { return ii_MaxN; };

   uint32_t         *GetPrimes(uint32_t &numberOfPrimes) { numberOfPrimes = ii_NumberOfPrimes; return ip_Primes; };
   uint16_t         *GetPrimeGaps(uint16_t &biggestGap) { biggestGap = ii_BiggestGap; return ip_PrimeGaps; };   
   
#if defined(USE_OPENCL) || defined(USE_METAL)
   uint32_t          GetMaxGpuSteps(void) { return ii_MaxGpuSteps; };
   uint32_t          GetMaxGpuFactors(void) { return ii_MaxGpuFactors; };
#endif

   bool              ReportFactor(uint64_t theFactor, uint32_t n);

   void              FillTerms(uint8_t *terms);
   
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
   void              VerifyFactor(uint64_t theFactor, uint32_t n);
   
   std::vector<bool> iv_Terms;

   uint32_t          ii_MinN;
   uint32_t          ii_MaxN;

   uint32_t         *ip_Primes;
   uint32_t          ii_NumberOfPrimes;

   uint16_t         *ip_PrimeGaps;
   uint16_t          ii_BiggestGap;
   
#if defined(USE_OPENCL) || defined(USE_METAL)
   uint32_t          ii_MaxGpuSteps;
   uint32_t          ii_MaxGpuFactors;
#endif
};

#endif

