#ifdef __cplusplus
#ifndef _CW_KERNEL_CL
#define _CW_KERNEL_CL
const char *cw_kernel= \
"#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable\n" \
"void collect_factor(int  n,\n" \
"int  c,\n" \
"long p,\n" \
"volatile __global int *factorsCount,\n" \
"__global long4 *factors);\n" \
"long expmod(long base,\n" \
"long exp,\n" \
"ulong thePrime,\n" \
"ulong magicNumber,\n" \
"ulong magicShift);\n" \
"long mulmod(long a,\n" \
"long b,\n" \
"ulong thePrime,\n" \
"ulong magicNumber,\n" \
"ulong magicShift);\n" \
"long compute_inverse(long a, long modulus);\n" \
"#define MAX_POWERS 50\n" \
"__kernel void cw_kernel(__global const  long  *primes,\n" \
"__global const ulong  *magicNumbers,\n" \
"__global const ulong  *magicShifts,\n" \
"__global const  int   *terms,\n" \
"volatile __global        int   *factorsCount,\n" \
"__global        long4 *factors)\n" \
"{\n" \
"int gid = get_global_id(0);\n" \
"int   termIndex;\n" \
"int   theN, prevN;\n" \
"int   power, idx;\n" \
"long  thePrime = primes[gid];\n" \
"ulong magicNumber = magicNumbers[gid];\n" \
"ulong magicShift = magicShifts[gid];\n" \
"long  rem = 0;\n" \
"ulong inverse;\n" \
"ulong powers[MAX_POWERS+1];\n" \
"inverse = compute_inverse(BASE, thePrime);\n" \
"theN = terms[0];\n" \
"rem = expmod(inverse, theN, thePrime, magicNumber, magicShift);\n" \
"#ifdef CHECK_WOODALL\n" \
"if (rem == theN)\n" \
"collect_factor(theN, -1, thePrime, factorsCount, factors);\n" \
"#endif\n" \
"#ifdef CHECK_CULLEN\n" \
"if (rem == thePrime - theN)\n" \
"collect_factor(theN, +1, thePrime, factorsCount, factors);\n" \
"#endif\n" \
"powers[1] = BASE;\n" \
"for (idx=2; idx<MAX_POWERS+1; idx++)\n" \
"{\n" \
"if ((idx > 2) && (BASE & 1))\n" \
"{\n" \
"if (idx & 1)\n" \
"continue;\n" \
"powers[idx] = mulmod(powers[idx-2], powers[2], thePrime, magicNumber, magicShift);\n" \
"}\n" \
"else\n" \
"{\n" \
"powers[idx] = mulmod(powers[idx-1], powers[1], thePrime, magicNumber, magicShift);\n" \
"}\n" \
"}\n" \
"prevN = terms[0];\n" \
"termIndex = 1;\n" \
"while (terms[termIndex] > 0)\n" \
"{\n" \
"theN = terms[termIndex];\n" \
"power = prevN - theN;\n" \
"while (power > MAX_POWERS)\n" \
"{\n" \
"rem = mulmod(rem, powers[MAX_POWERS], thePrime, magicNumber, magicShift);\n" \
"power -= MAX_POWERS;\n" \
"}\n" \
"rem = mulmod(rem, powers[power], thePrime, magicNumber, magicShift);\n" \
"#ifdef CHECK_WOODALL\n" \
"if (rem == theN)\n" \
"collect_factor(theN, -1, thePrime, factorsCount, factors);\n" \
"#endif\n" \
"#ifdef CHECK_CULLEN\n" \
"if (rem == thePrime - theN)\n" \
"collect_factor(theN, +1, thePrime, factorsCount, factors);\n" \
"#endif\n" \
"prevN = theN;\n" \
"termIndex++;\n" \
"};\n" \
"}\n" \
"void collect_factor(int  n,\n" \
"int  c,\n" \
"long p,\n" \
"volatile __global int *factorsCount,\n" \
"__global long4 *factors)\n" \
"{\n" \
"int old = atomic_inc(factorsCount);\n" \
"if (old >= MAX_FACTORS)\n" \
"return;\n" \
"factors[old].x = n;\n" \
"factors[old].y = c;\n" \
"factors[old].z = p;\n" \
"}\n" \
"long expmod(long base,\n" \
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
"long compute_inverse(long a, long p)\n" \
"{\n" \
"long     m = p;\n" \
"long     y = 0, x = 1;\n" \
"long     q, t;\n" \
"while (a > 1)\n" \
"{\n" \
"q = a / m;\n" \
"t = m;\n" \
"m = a % m;\n" \
"a = t;\n" \
"t = y;\n" \
"y = x - q * y;\n" \
"x = t;\n" \
"}\n" \
"if (x < 0)\n" \
"x += p;\n" \
"return x;\n" \
"}\n" \
;
#endif
#endif
