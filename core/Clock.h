/* Clock.h -- (c) Mark Rodenkirch, December 2013
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _CLOCK_H
#define _CLOCK_H

#include <stdio.h>

class Clock
{  
public:
   Clock();
   ~Clock(void);

   static uint64_t       GetCurrentMicrosecond(void);
   static uint64_t       GetThreadMicroseconds(void);
   static uint64_t       GetProcessMicroseconds(void);
};

#endif
