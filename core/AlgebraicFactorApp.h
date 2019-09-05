/* AlgebraicFactorApp.h -- (C) Mark Rodenkirch, January 2018

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _AlgebraicFactorApp_H
#define _AlgebraicFactorApp_H

#include <stdio.h>

#include "FactorApp.h"
#include "SharedMemoryItem.h"


class AlgebraicFactorApp : public FactorApp
{  
public:
   AlgebraicFactorApp(void);
   ~AlgebraicFactorApp(void);

protected:
   virtual void      EliminateTermsWithAlgebraicFactors(void) = 0;
   
   bool              IsGfnOrMersenneForm(uint64_t k, uint32_t base, int32_t c);
   
   void              GetRoot(uint64_t theNumber, uint64_t *root, uint32_t *power);
   
private:   
   // These variables are used by this class when looking for algebraic factors.
   vector<uint64_t>  iv_SmallPrimes;
   
   void              GetSmallPrimes(void);
   uint32_t          GetFactorList(uint64_t the_number, uint64_t *factor_list, uint32_t *power_list);
};

#endif
