/* FixedKBNWorker.h -- (C) Mark Rodenkirch, July 2017

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _FixedKBNWorker_H
#define _FixedKBNWorker_H

#include "FixedKBNApp.h"
#include "../core/Worker.h"

using namespace std;

class FixedKBNWorker : public Worker
{
public:
   FixedKBNWorker(uint32_t myId, App *theApp);

   ~FixedKBNWorker(void) {};

   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   void              CleanUp(void);

protected:

private:
   void              RemoveTerms(uint64_t prime, uint64_t b);
   void              VerifyFactor(uint64_t prime, int64_t c);

   FixedKBNApp      *ip_FixedKBNApp;

   uint64_t          il_KBpowN;
   uint64_t          il_K;
   uint32_t          ii_Base;
   uint32_t          ii_N;
   int64_t           il_MinC;
   int64_t           il_MaxC;
};

#endif

