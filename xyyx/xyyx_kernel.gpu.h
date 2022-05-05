#ifdef __cplusplus
#ifndef _XYYX_KERNEL_CL
#define _XYYX_KERNEL_CL
const char *xyyx_kernel= \
"#pragma OPENCL EXTENSION cl_khr_int64_extended_atomics: enable\n" \
"#define MAX_POWERS 50\n" \
"#if defined(USE_OPENCL)\n" \
"void collectFactor(uint    x,\n" \
"uint    y,\n" \
"ulong   p,\n" \
"volatile __global uint   *factorCount,\n" \
"__global ulong4 *factors);\n" \
"#else\n" \
"void collectFactor(uint        x,\n" \
"uint        y,\n" \
"ulong       p,\n" \
"volatile device atomic_int *factorCount,\n" \
"device ulong4     *factors);\n" \
"#endif\n" \
"void computeMagicValues(ulong thePrime,\n" \
"ulong *magicNumber,\n" \
"ulong *magicShift);\n" \
"ulong powmod(ulong base,\n" \
"ulong exp,\n" \
"ulong thePrime,\n" \
"ulong magicNumber,\n" \
"ulong magicShift);\n" \
"ulong mulmod(ulong a,\n" \
"ulong b,\n" \
"ulong thePrime,\n" \
"ulong magicNumber,\n" \
"ulong magicShift);\n" \
"__kernel void xyyx_kernel(__global const ulong  *primes,\n" \
"__global const uint   *terms,\n" \
"volatile __global       uint   *factorCount,\n" \
"__global       ulong4 *factors)\n" \
"{\n" \
"uint   gid = get_global_id(0);\n" \
"uint   termIndex;\n" \
"uint   currentX;\n" \
"uint   currentY, previousY;\n" \
"uint   pIndex;\n" \
"ulong  thePrime = primes[gid];\n" \
"ulong  magicNumber;\n" \
"ulong  magicShift;\n" \
"ulong  xPowY, xPowN;\n" \
"ulong  yPowX;\n" \
"ulong powers[MAX_POWERS];\n" \
"if (thePrime == 0)\n" \
"return;\n" \
"computeMagicValues(thePrime, &magicNumber, &magicShift);\n" \
"termIndex = 0;\n" \
"while (terms[termIndex] > 0)\n" \
"{\n" \
"currentX = terms[termIndex];\n" \
"termIndex++;\n" \
"previousY = 0;\n" \
"if (terms[termIndex] > 0)\n" \
"{\n" \
"currentY = terms[termIndex];\n" \
"xPowY = mulmod(currentX, currentX, thePrime, magicNumber, magicShift);\n" \
"powers[0] = 1;\n" \
"powers[1] = xPowY;\n" \
"for (pIndex=2; pIndex<=MAX_POWERS; pIndex++)\n" \
"powers[pIndex] = mulmod(powers[pIndex-1], xPowY, thePrime, magicNumber, magicShift);\n" \
"xPowY = powmod(currentX, currentY, thePrime, magicNumber, magicShift);\n" \
"previousY = currentY;\n" \
"while (terms[termIndex] > 0)\n" \
"{\n" \
"currentY = terms[termIndex];\n" \
"pIndex = (currentY - previousY) >> 1;\n" \
"if (pIndex >= MAX_POWERS)\n" \
"xPowN = powmod(currentX, currentY - previousY, thePrime, magicNumber, magicShift);\n" \
"else\n" \
"xPowN = powers[pIndex];\n" \
"xPowY = mulmod(xPowY, xPowN, thePrime, magicNumber, magicShift);\n" \
"yPowX = powmod(currentY, currentX, thePrime, magicNumber, magicShift);\n" \
"#ifdef IS_MINUS\n" \
"if (xPowY == yPowX)\n" \
"collectFactor(currentX, currentY, thePrime, factorCount, factors);\n" \
"#else\n" \
"if (xPowY + yPowX == thePrime)\n" \
"collectFactor(currentX, currentY, thePrime, factorCount, factors);\n" \
"#endif\n" \
"termIndex++;\n" \
"previousY = currentY;\n" \
"}\n" \
"}\n" \
"termIndex++;\n" \
"};\n" \
"}\n" \
"void computeMagicValues(ulong  thePrime,\n" \
"ulong *magicNumber,\n" \
"ulong *magicShift)\n" \
"{\n" \
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
"*magicNumber = (q2 + 1);\n" \
"*magicShift = (p - 64);\n" \
"}\n" \
"ulong powmod(ulong base,\n" \
"ulong exp,\n" \
"ulong thePrime,\n" \
"ulong magicNumber,\n" \
"ulong magicShift)\n" \
"{\n" \
"long x = base, y = 1;\n" \
"while (true)\n" \
"{\n" \
"if (exp & 1)\n" \
"y = mulmod(x, y, thePrime, magicNumber, magicShift);\n" \
"exp >>= 1;\n" \
"if (!exp)\n" \
"return y;\n" \
"x = mulmod(x, x, thePrime, magicNumber, magicShift);\n" \
"}\n" \
"}\n" \
"ulong mulmod(ulong a,\n" \
"ulong b,\n" \
"ulong thePrime,\n" \
"ulong magicNumber,\n" \
"ulong magicShift)\n" \
"{\n" \
"ulong xa_low, xa_high;\n" \
"long xa_rem, xa_quot;\n" \
"ulong x196m, x196m1, x196m2, x196h;\n" \
"xa_low = a * b;\n" \
"xa_high = mul_hi(a, b);\n" \
"x196m1 = mul_hi(xa_low, magicNumber);\n" \
"x196m2 = xa_high * magicNumber;\n" \
"x196h = mul_hi(xa_high, magicNumber);\n" \
"x196m = x196m1 + x196m2;\n" \
"if (x196m < x196m1) x196h++;\n" \
"xa_quot  = (x196m >> magicShift);\n" \
"xa_quot |= (x196h << (64 - magicShift));\n" \
"xa_rem = xa_low - (xa_quot * thePrime);\n" \
"if (xa_rem < 0) { xa_rem += thePrime; xa_quot -= 1; }\n" \
"return xa_rem;\n" \
"}\n" \
"#if defined(USE_OPENCL)\n" \
"void collectFactor(uint    x,\n" \
"uint    y,\n" \
"ulong   p,\n" \
"volatile __global uint   *factorCount,\n" \
"__global ulong4 *factors)\n" \
"#else\n" \
"void collectFactor(uint        x,\n" \
"uint        y,\n" \
"ulong       p,\n" \
"volatile device atomic_int *factorCount,\n" \
"device long4      *factors)\n" \
"#endif\n" \
"{\n" \
"#if defined(USE_OPENCL)\n" \
"int old = atomic_inc(factorCount);\n" \
"#else\n" \
"int old = atomic_fetch_add_explicit(factorCount, 1, memory_order_relaxed);\n" \
"#endif\n" \
"if (old >= D_MAX_FACTORS)\n" \
"return;\n" \
"factors[old].x = x;\n" \
"factors[old].y = y;\n" \
"factors[old].z = p;\n" \
"}\n" \
;
#endif
#endif
