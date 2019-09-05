/* AbstractWorker.cpp -- (C) Mark Rodenkirch, October 2018

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <stdint.h>

#include "AbstractWorker.h"
#include "../x86_asm/sse-asm-x86.h"

AbstractWorker::AbstractWorker(uint32_t myId, App *theApp, AbstractSubsequenceHelper *appHelper) : Worker(myId, theApp)
{
   ip_SierpinskiRieselApp = (SierpinskiRieselApp *) theApp;
   ip_AppHelper = appHelper;
   
   ii_Base = ip_SierpinskiRieselApp->GetBase();
   ii_MinN = ip_SierpinskiRieselApp->GetMinN();
   ii_MaxN = ip_SierpinskiRieselApp->GetMaxN();
   
   ip_Sequences = 0;
   ii_SequenceCount = 0;
   
   ip_Subsequences = 0;
   ii_SubsequenceCount = 0;
}
