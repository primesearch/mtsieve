/* GFNDivisorTester.h -- (C) Mark Rodenkirch, February 2021
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _GFNDivisorTester_H
#define _GFNDivisorTester_H

#include <gmp.h>
#include "../core/App.h"

class GFNDivisorTester 
{
public:
   GFNDivisorTester(App *theApp);

   ~GFNDivisorTester(void) {};

   void              StartedSieving(void);
   void              TestRemainingTerms(uint64_t totalTerms, uint64_t termsInChunk, uint64_t termCount);
   
protected:

private:
   bool              IsFermatDivisor(uint64_t k, uint32_t n);
   void              CheckRedc(mp_limb_t *xp, uint32_t xn, uint32_t b, uint32_t m, uint64_t k, uint32_t n);
   void              VerifyFactor(uint64_t thePrime, uint64_t k, uint32_t n);
   
   App              *ip_App;
   
   std::vector<std::vector<bool>>  iv_Terms;
   uint64_t          il_MinK;
   uint64_t          il_MaxK;
   uint32_t          ii_MinN;
   uint32_t          ii_MaxN;
   
   uint64_t          il_StartSievingUS;
   uint64_t          il_TotalTermsEvaluated;
};

#endif

