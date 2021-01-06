/* Parser.h -- (C) Mark Rodenkirch, February 2012

   This class provides the interface functions to parse strings as
   integers and verify ensure the value is in a range.
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _PARSER_H
#define _PARSER_H

#include <inttypes.h>

typedef enum { P_SUCCESS, P_FAILURE, P_OUT_OF_RANGE, P_UNSUPPORTED } parse_t;

class Parser
{  
public:
   Parser();
   ~Parser(void);

   static parse_t  Parse(const char *str, const char *valid, char &value);
   static parse_t  Parse(const char *str,  int32_t lo,  int32_t hi,  int32_t &value);
   static parse_t  Parse(const char *str, uint32_t lo, uint32_t hi, uint32_t &value);
   static parse_t  Parse(const char *str,  int64_t lo,  int64_t hi,  int64_t &value);
   static parse_t  Parse(const char *str, uint64_t lo, uint64_t hi, uint64_t &value);
   static parse_t  Parse(const char *str,   double lo,   double hi,   double &value);
};

#endif
