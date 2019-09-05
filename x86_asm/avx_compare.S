/* avx_compare.S -- (C) Mark Rodenkirch, May 2018

   All AVX functions require parameters to be aligned to a 32-byte
   boundary such as shown here:
      double __attribute__((aligned(32))) p[1];

   uint16_t  _avx_pos_compare_1v(double *comparator)
   uint16_t  _avx_pos_compare_16v(double *comparator)
      Return 16 bits where each bit indicates if the corresponding double
      in ymm12-ymm15 is equal to the comparator.
      
   uint16_t  _avx_neg_compare_1v(double *comparator, double *p)
   uint16_t  _avx_neg_compare_16v(double *comparator, double *p)
      Return 16 bits where each bit indicates if the corresponding double
      in ymm12-ymm15 is equal to p - comparator.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

  .text
  
#ifdef _WIN64
#define ARG1      %rcx
#define ARG2      %rdx
#else
#define ARG1      %rdi
#define ARG2      %rsi
#endif

  .p2align 4,,15

  .globl _avx_neg_compare_1v
  .globl avx_neg_compare_1v

_avx_neg_compare_1v:
avx_neg_compare_1v:
   vbroadcastsd   0(ARG1), %ymm8
   vbroadcastsd   0(ARG1), %ymm9
   vbroadcastsd   0(ARG1), %ymm10
   vbroadcastsd   0(ARG1), %ymm11
   jmp            avx_neg_compare
   
  .globl _avx_neg_compare_16v
  .globl avx_neg_compare_16v

_avx_neg_compare_16v:
avx_neg_compare_16v:
   vmovapd	      0x00(ARG1), %ymm8   
   vmovapd	      0x20(ARG1), %ymm9   
   vmovapd	      0x40(ARG1), %ymm10   
   vmovapd	      0x60(ARG1), %ymm11
   
avx_neg_compare:
   vmovapd        0x00(ARG2), %ymm4
   vmovapd        0x20(ARG2), %ymm5
   vmovapd        0x40(ARG2), %ymm6
   vmovapd        0x60(ARG2), %ymm7
   
   vsubpd         %ymm8 , %ymm4, %ymm0
   vsubpd         %ymm9 , %ymm5, %ymm1
   vsubpd         %ymm10, %ymm6, %ymm2
   vsubpd         %ymm11, %ymm7, %ymm3
   jmp            avx_compare
  
  .globl _avx_pos_compare_1v
  .globl avx_pos_compare_1v

_avx_pos_compare_1v:
avx_pos_compare_1v: 
   vbroadcastsd   0(ARG1), %ymm0
   vbroadcastsd   0(ARG1), %ymm1
   vbroadcastsd   0(ARG1), %ymm2
   vbroadcastsd   0(ARG1), %ymm3
   jmp            avx_pos_compare
   
  .globl _avx_pos_compare_16v
  .globl avx_pos_compare_16v

_avx_pos_compare_16v:
avx_pos_compare_16v:
   vmovapd	      0x00(ARG1), %ymm0   
   vmovapd	      0x20(ARG1), %ymm1   
   vmovapd	      0x40(ARG1), %ymm2   
   vmovapd	      0x60(ARG1), %ymm3
   
avx_pos_compare:
avx_compare:
   vcmppd         $0, %ymm12, %ymm0, %ymm4
   vcmppd         $0, %ymm13, %ymm1, %ymm5
   vcmppd         $0, %ymm14, %ymm2, %ymm6
   vcmppd         $0, %ymm15, %ymm3, %ymm7

   // Note that the bits are in reverse sequence
   vmovmskpd      %ymm4, %rax
   vmovmskpd      %ymm5, %r8
   vmovmskpd      %ymm6, %r9
   vmovmskpd      %ymm7, %rdx

   shl            $12, %rax
   shl            $8, %r8
   shl            $4, %r9
      
   or             %r8, %rax
   or             %r9, %rdx
   or             %rdx, %rax
      
   ret

