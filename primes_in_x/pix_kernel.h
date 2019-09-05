#ifdef __cplusplus
#ifndef _PIX_KERNEL_CL
#define _PIX_KERNEL_CL
const char *pix_kernel= \
"#pragma OPENCL EXTENSION cl_khr_int64_extended_atomics: enable\n" \
"__kernel void pix_kernel(__global const  long2 *primes,\n" \
"__global const ulong2 *magicNumbers,\n" \
"__global const ulong2 *magicShifts,\n" \
"__global        long2 *rems,\n" \
"__global       uint   *digits,\n" \
"volatile __global        long  *factors)\n" \
"{\n" \
"int gid = get_global_id(0);\n" \
"long maxMul;\n" \
"ulong2 multiplier;\n" \
"ulong2 quotient, xa_temp, xa_low, xa_high;\n" \
"ulong2 x196m1, x196m2, x196h, x196m;\n" \
"int  index = 0;\n" \
"long2  pl = primes[gid];\n" \
"long2  xa_rem = rems[gid];\n" \
"ulong2 magicNumber = magicNumbers[gid];\n" \
"ulong2 magicShift = magicShifts[gid];\n" \
"if (pl.x == 0) return;\n" \
"if (pl.y == 0) pl.y = pl.x;\n" \
"multiplier.x = digits[1];\n" \
"multiplier.y = digits[1];\n" \
"maxMul = pl.y / multiplier.y;\n" \
"index = 2;\n" \
"while (digits[index] < multiplier.x)\n" \
"{\n" \
"if (xa_rem.x < maxMul && xa_rem.y < maxMul)\n" \
"{\n" \
"xa_rem.x = xa_rem.x * multiplier.x;\n" \
"xa_rem.y = xa_rem.y * multiplier.y;\n" \
"}\n" \
"else\n" \
"{\n" \
"xa_temp.x = xa_rem.x;\n" \
"xa_temp.y = xa_rem.y;\n" \
"xa_low = xa_temp * multiplier;\n" \
"xa_high = mul_hi(xa_temp, multiplier);\n" \
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
"}\n" \
"xa_rem.x += digits[index];\n" \
"xa_rem.y += digits[index];\n" \
"if (xa_rem.x >= pl.x) xa_rem.x -= pl.x;\n" \
"if (xa_rem.y >= pl.y) xa_rem.y -= pl.y;\n" \
"if (xa_rem.x == 0) atom_min(&factors[index], pl.x);\n" \
"if (xa_rem.y == 0) atom_min(&factors[index], pl.y);\n" \
"index++;\n" \
"}\n" \
"rems[gid] = xa_rem;\n" \
"}\n" \
;
#endif
#endif
