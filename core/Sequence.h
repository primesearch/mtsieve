

#ifndef _SubSequence_H
#define _SubSequence_H

typedef struct
{
  uint32_t seq;
  uint32_t d;
  uint32_t filter;
  uint_fast32_t *M;
  uint32_t mcount;
  uint32_t a,b;
} subseq_t;

#endif