/* SmarandacheWorker.h -- (C) Mark Rodenkirch, January 2022

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _SmarandacheWorker_H
#define _SmarandacheWorker_H

#include "SmarandacheApp.h"
#include "../core/Worker.h"

using namespace std;

class SmarandacheWorker : public Worker
{
public:
   SmarandacheWorker(uint32_t myId, App *theApp);

   ~SmarandacheWorker() {};

   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   void              CleanUp(void);

protected:
   bool              TestSixDigitN(void);
   bool              TestSevenDigitN(void);
   
   SmarandacheApp   *ip_SmarandacheApp;
   
   uint32_t          ii_MinN;
   uint32_t          ii_MaxN;
   
   terms_t          *ip_Terms;
};

#endif

