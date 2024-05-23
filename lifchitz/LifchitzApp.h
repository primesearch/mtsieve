/* LifchitzApp.h -- (C) Mark Rodenkirch, April 2024
 
   This class inherits from App.h and has the implementation for this project
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _LifchitzApp_H
#define _LifchitzApp_H

#include <vector>
#include <set>

#include "../core/FactorApp.h"

#define P_ONE     0x01
#define M_ONE     0x02

typedef struct {
   uint32_t    x;
   uint32_t    y;
   uint8_t     signs;
} term_t;

class LifchitzApp : public FactorApp
{
public:
   LifchitzApp(void);

   ~LifchitzApp(void) {};

   void              Help(void);
   void              AddCommandLineOptions(std::string &shortOpts, struct option *longOpts);
   parse_t           ParseOption(int opt, char *arg, const char *source);
   void              ValidateOptions(void);
   bool              ApplyFactor(uint64_t theFactor, const char *term);
   void              GetExtraTextForSieveStartedMessage(char *extraText, uint32_t maxTextLength);
   
   bool              ReportFactor(uint64_t theFactor, uint32_t x, uint32_t y, int32_t sign);

   uint32_t          GetMinX(void) { return ii_MinX; };
   uint32_t          GetMaxX(void) { return ii_MaxX; };
   uint32_t          GetMinY(void) { return ii_MinY; };
   uint32_t          GetMaxY(void) { return ii_MaxY; };
      
#if defined(USE_OPENCL) || defined(USE_METAL)
   uint32_t          GetXPerChunk(void) { return ii_XPerChunk; };
   uint32_t          GetMaxGpuFactors(void) { return ii_MaxGpuFactors; };
#endif

protected:
   void              PreSieveHook(void) {};
   bool              PostSieveHook(void) { return true; };
   
   void              NotifyAppToRebuild(uint64_t largestPrimeTested) {};
   
   Worker           *CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested);
   
   void              ProcessInputTermsFile(bool haveBitMap);
   bool              IsWritingOutputTermsFile(void){ return true; };
   void              WriteOutputTermsFile(uint64_t largestPrime);
   void              OuptutAdditionalConsoleMessagesUponFinish(void) {};

private:
   void              SetInitialTerms(void);
   void              CollapseTerms(void);
   bool              RemoveTerm(uint64_t theFactor, term_t *term, int32_t sign, bool logFactor);
   void              VerifyFactor(uint64_t theFactor, uint32_t x, uint32_t y, int32_t sign);
      
   bool              ib_IsPlus;
   bool              ib_IsMinus;
   uint32_t          ii_MinX;
   uint32_t          ii_MaxX;
   uint32_t          ii_MinY;
   uint32_t          ii_MaxY;
   uint64_t          il_TermArraySize;
   
   term_t           *ip_Terms;

#if defined(USE_OPENCL) || defined(USE_METAL)
   uint32_t          ii_XPerChunk;
   uint32_t          ii_MaxGpuFactors;
#endif
};

#endif

