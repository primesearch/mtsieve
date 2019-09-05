#ifdef __cplusplus
#ifndef _GFN_KERNEL_CL
#define _GFN_KERNEL_CL
const char *gfn_kernel= \
"void compute_magic(long theDivisor,\n" \
"ulong *magicNumber,\n" \
"ulong *magicShift);\n" \
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
"__kernel void gfn_kernel(__global const long *primes,\n" \
"__global      ulong *ks)\n" \
"{\n" \
"int    gid = get_global_id(0);\n" \
"long   p = primes[gid];\n" \
"long   k, rem;\n" \
"ulong  magicNumber, magicShift;\n" \
"ulong  mN, mS;\n" \
"int    n, index;\n" \
"compute_magic(p, &magicNumber, &magicShift);\n" \
"k = (1 + p) >> 1;\n" \
"rem = powmod(k, N_MIN, p, magicNumber, magicShift);\n" \
"index = gid * N_COUNT;\n" \
"k = p - rem;\n" \
"for (n=0; n<N_COUNT; n++)\n" \
"{\n" \
"if (k & 1)\n" \
"{\n" \
"if (k <= K_MAX && k >= K_MIN)\n" \
"ks[index + n] = k;\n" \
"k += p;\n" \
"}\n" \
"else\n" \
"ks[index + n] = 0;\n" \
"k >>= 1;\n" \
"}\n" \
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
