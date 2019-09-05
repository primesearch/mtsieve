#ifdef __cplusplus
#ifndef _MAGIC_CL
#define _MAGIC_CL
const char *magic= \
"__kernel void magic_kernel(__global const long *primes,\n" \
"__global      ulong *magicNumbers,\n" \
"__global      ulong *magicShifts)\n" \
"{\n" \
"int   gid = get_global_id(0);\n" \
"long   thePrime = primes[gid];\n" \
"ulong  magicNumber, magicShift;\n" \
"ulong  mN, mS;\n" \
"ulong two63 = 0x8000000000000000;\n" \
"ulong d = thePrime;\n" \
"ulong t = two63;\n" \
"ulong anc = t - 1 - t%d;\n" \
"ulong p = 63;\n" \
"ulong q1 = two63/anc;\n" \
"ulong r1 = two63 - q1*anc;\n" \
"ulong q2 = two63/d;\n" \
"ulong r2 = two63- q2*d;\n" \
"ulong delta;\n" \
"do {\n" \
"p = p + 1;\n" \
"q1 = 2*q1;\n" \
"r1 = 2*r1;\n" \
"if (r1 >= anc) {\n" \
"q1 = q1 + 1;\n" \
"r1 = r1 - anc;\n" \
"}\n" \
"q2 = 2*q2;\n" \
"r2 = 2*r2;\n" \
"if (r2 >= d) {\n" \
"q2 = q2 + 1;\n" \
"r2 = r2 - d;\n" \
"}\n" \
"delta = d - r2;\n" \
"} while (q1 < delta || (q1 == delta && r1 == 0));\n" \
"magicNumbers[gid] = (q2 + 1);\n" \
"magicShifts[gid] = (p - 64);\n" \
"}\n" \
;
#endif
#endif
