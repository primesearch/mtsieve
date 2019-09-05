#ifdef __cplusplus
#ifndef _MF_KERNEL_CL
#define _MF_KERNEL_CL
const char *mf_kernel= \
"#pragma OPENCL EXTENSION cl_khr_int64_extended_atomics: enable\n" \
"__kernel void mf_kernel(__global const  long2 *primes,\n" \
"__global const ulong2 *magicNumbers,\n" \
"__global const ulong2 *magicShifts,\n" \
"__global        long2 *rems,\n" \
"__global        int   *nRange,\n" \
"volatile __global        long  *minusFactors,\n" \
"volatile __global        long  *plusFactors)\n" \
"{\n" \
"int gid = get_global_id(0);\n" \
"long n_min = nRange[0];\n" \
"long n_max = nRange[1];\n" \
"long n_start = n_min;\n" \
"long n_idx = 0;\n" \
"long maxMul;\n" \
"ulong2 n_mul;\n" \
"ulong2 quotient, xa_temp, xa_low, xa_high;\n" \
"ulong2 x196m1, x196m2, x196h, x196m;\n" \
"long2  pl = primes[gid];\n" \
"long2  xa_rem = rems[gid];\n" \
"ulong2 magicNumber = magicNumbers[gid];\n" \
"ulong2 magicShift = magicShifts[gid];\n" \
"if (pl.x == 0) return;\n" \
"if (pl.y == 0) pl.y = pl.x;\n" \
"maxMul = pl.y / n_max;\n" \
"n_mul.x = n_min;\n" \
"n_mul.y = n_min;\n" \
"do\n" \
"{\n" \
"if (xa_rem.x < maxMul && xa_rem.y < maxMul)\n" \
"{\n" \
"xa_rem.x = xa_rem.x * n_mul.x;\n" \
"xa_rem.y = xa_rem.y * n_mul.y;\n" \
"}\n" \
"else\n" \
"{\n" \
"xa_temp.x = xa_rem.x;\n" \
"xa_temp.y = xa_rem.y;\n" \
"xa_low = xa_temp * n_mul;\n" \
"xa_high = mul_hi(xa_temp, n_mul);\n" \
"x196m1 = mul_hi(xa_low, magicNumber);\n" \
"x196m2 = xa_high * magicNumber;\n" \
"x196h = mul_hi(xa_high, magicNumber);\n" \
"x196m = x196m1 + x196m2;\n" \
"if (x196m.x < x196m1.x) x196h.x++;\n" \
"if (x196m.y < x196m1.y) x196h.y++;\n" \
"quotient = (x196m >> magicShift);\n" \
"quotient |= (x196h << (64 - magicShift));\n" \
"xa_rem.x = (long) (xa_low.x) - ((long) quotient.x * pl.x);\n" \
"xa_rem.y = (long) (xa_low.y) - ((long) quotient.y * pl.y);\n" \
"if (xa_rem.x < 0) xa_rem.x += pl.x;\n" \
"if (xa_rem.y < 0) xa_rem.y += pl.y;\n" \
"if (xa_rem.x ==      1) atom_min(&minusFactors[n_idx], pl.x);\n" \
"if (xa_rem.y ==      1) atom_min(&minusFactors[n_idx], pl.y);\n" \
"if (xa_rem.x == pl.x-1) atom_min(&plusFactors[n_idx], pl.x);\n" \
"if (xa_rem.y == pl.y-1) atom_min(&plusFactors[n_idx], pl.y);\n" \
"}\n" \
"n_mul.x += D_ADDER;\n" \
"n_mul.y += D_ADDER;\n" \
"n_min += D_ADDER;\n" \
"n_idx++;\n" \
"#ifdef D_INVERSE\n" \
"if (n_min == n_skip)\n" \
"{\n" \
"n_mul.x += D_ADDER;\n" \
"n_mul.y += D_ADDER;\n" \
"n_min += D_ADDER;\n" \
"n_skip += D_ADDER;\n" \
"}\n" \
"#endif\n" \
"} while (n_min < n_max);\n" \
"rems[gid] = xa_rem;\n" \
"}\n" \
;
#endif
#endif
