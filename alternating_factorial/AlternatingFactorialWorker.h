/* AlternatingFactorialWorker.h -- (C) Mark Rodenkirch, July 2017

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _AlternatingFactorialWorker_H
#define _AlternatingFactorialWorker_H

#include "AlternatingFactorialApp.h"
#include "../core/Worker.h"

using namespace std;

class AlternatingFactorialWorker : public Worker
{
public:
   AlternatingFactorialWorker(uint32_t myId, App *theApp);

   ~AlternatingFactorialWorker(void) {};

   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   void              CleanUp(void);

protected:

private:
   AlternatingFactorialApp      *ip_AlternatingFactorialApp;
   
   void              ExtractFactors(uint64_t p);
};

#endif

