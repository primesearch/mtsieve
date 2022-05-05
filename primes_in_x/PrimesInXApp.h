/* PrimesInXApp.h -- (C) Mark Rodenkirch, October 2012

   This class inherits from App.h and has the implementation for this project

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _PrimesInXApp_H
#define _PrimesInXApp_H

#include "../core/FactorApp.h"
#include "../core/SharedMemoryItem.h"

class PrimesInXApp : public FactorApp
{
public:
   PrimesInXApp(void);

   ~PrimesInXApp(void);

   void              Help(void);
   void              AddCommandLineOptions(std::string &shortOpts, struct option *longOpts);
   parse_t           ParseOption(int opt, char *arg, const char *source);
   void              ValidateOptions(void);
   bool              ApplyFactor(uint64_t theFactor, const char *term);
   void              GetExtraTextForSieveStartedMessage(char *extraText);
   
#if defined(USE_OPENCL) || defined(USE_METAL)
   uint32_t          GetMaxGpuSteps(void) { return ii_MaxGpuSteps; };
   uint32_t          GetMaxGpuFactors(void) { return ii_MaxGpuFactors; };
#endif

   uint32_t          GetMinLength(void) { return ii_MinLength; };
   uint32_t          GetMinLengthRemaining(void) { return ii_MinLengthRemaining; };
   uint32_t          GetMaxLength(void) { return ii_MaxLength; };

   uint32_t         *Get1DigitTerms(void) { return ii_e1Terms; };
   uint32_t         *Get3DigitTermsCopy(void);
   uint32_t         *Get6DigitTermsCopy(void);
   uint32_t         *Get9DigitTermsCopy(void);
   
   bool              ReportFactor(uint64_t theFactor, uint32_t n);
   void              ReportPrime(uint64_t thePrime, uint32_t n);

protected:
   void              PreSieveHook(void) {};
   bool              PostSieveHook(void) { return true; };
   
   void              NotifyAppToRebuild(uint64_t largestPrimeTested) {};
   
   void              ProcessInputTermsFile(bool haveBitMap);
   bool              IsWritingOutputTermsFile(void){ return true; };
   void              WriteOutputTermsFile(uint64_t largestPrime);
   
   Worker           *CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested);

private:
   void              ProcessInputStringFile(void);
   void              BuildTerms(char *inputTerm);
   void              VerifyFactor(uint64_t theFactor, uint32_t termLength);
   
   std::vector<bool> iv_Terms;
   
   std::string       is_SearchString;
   std::string       is_FullTerm;
  std:: string       is_StringFileName;

   uint32_t          ii_MaxGpuSteps;
   uint32_t          ii_MaxGpuFactors;
   uint32_t          ii_MinLength;
   uint32_t          ii_MinLengthRemaining;
   uint32_t          ii_MaxLength;
   
   uint32_t         *ii_e1Terms;                        // The string as an array of terms < 10
   uint32_t         *ii_e3Terms;                        // The string as an array of terms < 1000
   uint32_t         *ii_e6Terms;                        // The string as an array of terms < 1000000 
   uint32_t         *ii_e9Terms;                        // The string as an array of trems < 1000000000
   uint32_t          ii_e1TermCount;
   uint32_t          ii_e3TermCount;
   uint32_t          ii_e6TermCount;
   uint32_t          ii_e9TermCount;
};

#endif

