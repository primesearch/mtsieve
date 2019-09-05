#ifdef __cplusplus
#ifndef _XYYX_KERNEL_CL
#define _XYYX_KERNEL_CL
const char *xyyx_kernel= \
"#pragma OPENCL EXTENSION cl_khr_int64_extended_atomics: enable\n" \
"#define MAX_POWERS 50\n" \
"void collect_factor(int  x,\n" \
"int  y,\n" \
"int  c,\n" \
"long p,\n" \
"volatile __global int *counter,\n" \
"__global long4 *factors);\n" \
"long powmod(long base,\n" \
"long exp,\n" \
"ulong thePrime,\n" \
"ulong magicNumber,\n" \
"ulong magicShift);\n" \
"long mulmod(long a,\n" \
"long b,\n" \
"ulong thePrime,\n" \
"ulong magicNumber,\n" \
"ulong magicShift);\n" \
"__kernel void xyyx_kernel(__global const  long  *primes,\n" \
"__global const ulong  *magicNumbers,\n" \
"__global const ulong  *magicShifts,\n" \
"__global const  int   *terms,\n" \
"volatile  __global        long  *minusFactors,\n" \
"volatile  __global        long  *plusFactors)\n" \
"{\n" \
"int gid = get_global_id(0);\n" \
"int   termIndex;\n" \
"int   currentX;\n" \
"int   currentY, previousY;\n" \
"int   pIndex;\n" \
"long  thePrime = primes[gid];\n" \
"ulong magicNumber = magicNumbers[gid];\n" \
"ulong magicShift = magicShifts[gid];\n" \
"long  xPowY, xPowN;\n" \
"long  yPowX;\n" \
"ulong powers[MAX_POWERS];\n" \
"if (thePrime == 0)\n" \
"return;\n" \
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
"if (pIndex > MAX_POWERS)\n" \
"xPowN = powmod(currentX, currentY - previousY, thePrime, magicNumber, magicShift);\n" \
"else\n" \
"xPowN = powers[pIndex];\n" \
"xPowY = mulmod(xPowY, xPowN, thePrime, magicNumber, magicShift);\n" \
"yPowX = powmod(currentY, currentX, thePrime, magicNumber, magicShift);\n" \
"if (xPowY == yPowX) atom_min(&minusFactors[termIndex], thePrime);\n" \
"if (xPowY + yPowX == thePrime) atom_min(&plusFactors[termIndex], thePrime);\n" \
"termIndex++;\n" \
"previousY = currentY;\n" \
"}\n" \
"}\n" \
"termIndex++;\n" \
"};\n" \
"}\n" \
"long powmod(long base,\n" \
"long exp,\n" \
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
"long mulmod(long a,\n" \
"long b,\n" \
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
;
#endif
#endif
