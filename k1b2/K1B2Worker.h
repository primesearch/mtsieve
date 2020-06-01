/* K1B2Worker.h -- (C) Mark Rodenkirch, June 2020

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _K1B2Worker_H
#define _K1B2Worker_H

#include "K1B2App.h"
#include "../core/Worker.h"

using namespace std;

class K1B2Worker : public Worker
{
public:
   K1B2Worker(uint32_t myId, App *theApp);

   ~K1B2Worker(void) {};

   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   void              CleanUp(void);

protected:

private:
   void              RemoveTermsSmallP(uint64_t prime, uint32_t n, uint64_t twoExpN);
   void              RemoveTermsLargeP(uint64_t prime, uint32_t n, uint64_t twoExpN);
   void              VerifyFactor(uint64_t prime, uint32_t n, int64_t c);

   K1B2App          *ip_K1B2App;

   uint32_t          ii_MinN;
   uint32_t          ii_MaxN;
   int64_t           il_MinC;
   int64_t           il_MaxC;
};

#endif

