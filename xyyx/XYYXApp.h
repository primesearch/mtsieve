/* XYYXApp.h -- (C) Mark Rodenkirch, May 2014

   This class inherits from App.h and has the implementation for this project
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _XYYXApp_H
#define _XYYXApp_H

#include <vector>
#include <set>

#include "../core/FactorApp.h"

typedef struct {
   uint32_t    y;
   uint64_t   *fpuRemainders;
   double     *avxRemainders;
} powerofx_t;

typedef struct {
   uint32_t    x;
   powerofx_t *powerOfX;
} powerofy_t;

typedef struct {
   uint32_t    base;
   uint32_t    powerCount;
   uint32_t    maxPowerDiff;
   uint32_t   *powerIndices;
   uint64_t   *fpuRemaindersPtr;
   double     *avxRemaindersPtr;
   powerofx_t *powersOfX;
   powerofy_t *powersOfY;
} base_t;

typedef struct {
   base_t     *xPowY;
   base_t     *yPowX;
} bases_t;

typedef struct {
   uint32_t    x;
   uint32_t    y;
   bool        haveFactor;
} term_t;

typedef struct {
   uint32_t    x;
   uint32_t    y;
} gputerm_t;

class XYYXApp : public FactorApp
{
public:
   XYYXApp(void);

   ~XYYXApp(void) {};

   void              Help(void);
   void              AddCommandLineOptions(std::string &shortOpts, struct option *longOpts);
   parse_t           ParseOption(int opt, char *arg, const char *source);
   void              ValidateOptions(void);
   bool              ApplyFactor(uint64_t theFactor, const char *term);
   void              GetExtraTextForSieveStartedMessage(char *extraText, uint32_t maxTextLength);
   
   bool              ReportFactor(uint64_t theFactor, uint32_t x, uint32_t y);
   
   uint32_t          GetMinX(void) { return ii_MinX; };
   uint32_t          GetMaxX(void) { return ii_MaxX; };
   uint32_t          GetMinY(void) { return ii_MinY; };
   uint32_t          GetMaxY(void) { return ii_MaxY; };
   uint64_t          GetTermCount(void) { return il_TermCount; }
   uint32_t          GetXCount(void) { return (ii_MaxX - ii_MinX + 1); };
   uint32_t          GetYCount(void) { return (ii_MaxY - ii_MinY + 1); };

   bool              UseAvxIfAvailable(void) { return ib_UseAvx; };
   bool              IsPlus(void) { return ib_IsPlus; };
   bool              IsMinus(void) { return ib_IsMinus; };

   void              GetTerms(uint32_t fpuRemaindersCount, uint32_t avxRemaindersCount, bases_t *bases);
   term_t           *GetSparseTerms(void);
   
#if defined(USE_OPENCL) || defined(USE_METAL)
   uint32_t          GetNumberOfGroups(void);
   uint32_t          GetGroupedTerms(uint32_t *terms);
   gputerm_t        *GetSparseGroupedTerms(void);
   
   uint32_t          GetMaxGpuSteps(void) { return ii_MaxGpuSteps; };
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
   void              VerifyFactor(uint64_t theFactor, uint32_t x, uint32_t y);
   
   std::vector<bool>  iv_Terms;
   term_t            *ip_Terms;
   
   bool              ib_Sparse;
   bool              ib_UseAvx;
   bool              ib_IsPlus;
   bool              ib_IsMinus;
   uint32_t          ii_MinX;
   uint32_t          ii_MaxX;
   uint32_t          ii_MinY;
   uint32_t          ii_MaxY;
   uint64_t          il_InitialSparseTermCount;

#if defined(USE_OPENCL) || defined(USE_METAL)
   uint32_t          ii_MaxGpuSteps;
   uint32_t          ii_MaxGpuFactors;
#endif
};

#endif

