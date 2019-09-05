#ifdef __cplusplus
#ifndef _MAGIC_CL
#define _MAGIC_CL
const char* magic= \
"void compute_magic(long theDivisor,\n" \
"ulong *magicNumber,\n" \
"ulong *magicShift);\n" \
"__kernel void magic_kernel(__global const long2 *primes,\n" \
"__global      ulong2 *magicNumbers,\n" \
"__global      ulong2 *magicShifts)\n" \
"{\n" \
"int gid = get_global_id(0);\n" \
"long2 pl = primes[gid];\n" \
"ulong2 magicNumber, magicShift;\n" \
"ulong  mN, mS;\n" \
"if (pl.x == 0) return;\n" \
"if (pl.y == 0) pl.y = pl.x;\n" \
"compute_magic(pl.x, &mN, &mS);\n" \
"magicNumber.x = mN;\n" \
"magicShift.x = mS;\n" \
"compute_magic(pl.y, &mN, &mS);\n" \
"magicNumber.y = mN;\n" \
"magicShift.y = mS;\n" \
"magicNumbers[gid] = magicNumber;\n" \
"magicShifts[gid] = magicShift;\n" \
"}\n" \
"void compute_magic(long theDivisor,\n" \
"ulong *magicNumber,\n" \
"ulong *magicShift)\n" \
"{\n" \
"ulong two63 = 0x8000000000000000;\n" \
"ulong d = theDivisor;\n" \
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
"*magicNumber = (q2 + 1);\n" \
"*magicShift = (p - 64);\n" \
"}\n" \
;
#endif
#endif
