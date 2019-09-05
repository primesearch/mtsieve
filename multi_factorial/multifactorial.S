/* multifactorial.S -- (C) Mark Rodenkirch, January 2018

   int multifactorial(uint32_t seed, uint32_t mf, const uint32_t max, uint64_t *P);
      Returns if any value n for n!mf = -1 or +1 mod p
         p is a array of 4 primes.

      Assumes FPU is set to double extended precision and round to zero.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

  .text
  
#ifdef _WIN64
#define ARG1      %rcx
#define ARG2      %rdx
#define ARG3      %r8
#define ARG4      %r9
#define mf      0(%rsp)
#define term    8(%rsp)
#define max    16(%rsp)
#define r1     24(%rsp)
#define r2     32(%rsp)
#define r3     40(%rsp)
#define r4     48(%rsp)
#define t1     56(%rsp)
#define t2     64(%rsp)
#define t3     72(%rsp)
#define t4     80(%rsp)
#define p1m1   88(%rsp)
#define p2m1   96(%rsp)
#define p3m1  104(%rsp)
#define p4m1  112(%rsp)
#else
#define ARG1      %rdi
#define ARG2      %rsi
#define ARG3      %rdx
#define ARG4      %rcx
#define mf     -0(%rsp)
#define term   -8(%rsp)
#define max   -16(%rsp)
#define r1    -24(%rsp)
#define r2    -32(%rsp)
#define r3    -40(%rsp)
#define r4    -48(%rsp)
#define t1    -56(%rsp)
#define t2    -64(%rsp)
#define t3    -72(%rsp)
#define t4    -80(%rsp)
#define p1m1  -88(%rsp)
#define p2m1  -96(%rsp)
#define p3m1 -104(%rsp)
#define p4m1 -112(%rsp)
#endif

#define p1      0(%rsi)
#define p2      8(%rsi)
#define p3     16(%rsi)
#define p4     24(%rsi)

  .globl _multifactorial
  .globl multifactorial

  .p2align 4,,15

_multifactorial:
multifactorial:
    
   push    %rbp
   push    %rbx
   push    %r12
   push    %r13
   push    %r14
   push    %r15
  
#ifdef _WIN64
   push    %rsi
   push    %rdi
   
   sub     $160, %rsp
#endif
   mov     ARG1, r1
   mov     ARG1, r2
   mov     ARG1, r3
   mov     ARG1, r4
   mov     ARG1, %rbp
   
   mov     ARG2, mf
   mov     ARG3, max
   mov     ARG4, %rsi
   
// Compute 1.0 / p for each p
   fildll  p4
   fld1                       // put 1.0 into %st(0)
   fdivp   %st(0), %st(1)     // compute 1.0 / p4
   
   fildll  p3
   fld1                       // put 1.0 into %st(0)
   fdivp   %st(0), %st(1)     // compute 1.0 / p3
   
   fildll  p2
   fld1                       // put 1.0 into %st(0)
   fdivp   %st(0), %st(1)     // compute 1.0 / p2
   
   fildll  p1
   fld1                       // put 1.0 into %st(0)
   fdivp   %st(0), %st(1)     // compute 1.0 / p1

   mov     p1, %r8
   mov     p2, %r9
   mov     p3, %r10
   mov     p4, %r11
           
   sub     $1, %r8
   sub     $1, %r9
   sub     $1, %r10
   sub     $1, %r11
   
   // Now we have p-1 for the 4 input values
   mov     %r8, p1m1
   mov     %r9, p2m1
   mov     %r10, p3m1
   mov     %r11, p4m1
   
   // Set the return code here.  Assume we found a factor
   mov     $1, %rax
   jmp     next_term

mf_loop:
   // add term off of the FP stack
   mov     %rbp, term
   fildll  term
   
   fildll  r1
   fmul    %st(1), %st(0)
   fmul    %st(2), %st(0)
   fistpll t1                 // t1 = trunc(r1*term / p1)

   fildll  r2
   fmul    %st(1), %st(0)
   fmul    %st(3), %st(0)
   fistpll t2                 // t2 = trunc(r2*term / p2)

   fildll  r3
   fmul    %st(1), %st(0)
   fmul    %st(4), %st(0)
   fistpll t3                 // t3 = trunc(r3*term / p3)
   
   fildll  r4
   fmul    %st(1), %st(0)
   fmul    %st(5), %st(0)
   fistpll t4                 // t4 = trunc(r4*term / p4)

   // remove term from the FP stack
   fstp    %st(0)
   
   mov     r1, %r8
   mov     r2, %r9
   mov     r3, %r10
   mov     r4, %r11
   
   // Yes, this can be a 124-bit product, but we only care
   // about the lower 62 bits because p < 2^62.
   imul    term, %r8          // r8  = r1*term
   imul    term, %r9          // r9  = r2*term
   imul    term, %r10         // r10 = r3*term
   imul    term, %r11         // r11 = r4*term

   mov     t1, %r12
   mov     t2, %r13
   mov     t3, %r14
   mov     t4, %r15

   imul    p1, %r12           // r12 = p1 * trunc(r1*term / p1) = r1*term - (r1*term % p1)
   imul    p2, %r13           // r13 = p2 * trunc(r2*term / p2) = r2*term - (r2*term % p2)
   imul    p3, %r14           // r14 = p3 * trunc(r3*term / p3) = r3*term - (r3*term % p3)
   imul    p4, %r15           // r15 = p4 * trunc(r4*term / p4) = r4*term - (r4*term % p4)
   
   sub     %r12, %r8          // r8  = r1*term - (r1*term - (r1*term % p1)) = r1*term % p1
   sub     %r13, %r9          // r9  = r2*term - (r2*term - (r2*term % p2)) = r2*term % p2
   sub     %r14, %r10         // r10 = r3*term - (r3*term - (r3*term % p3)) = r3*term % p3
   sub     %r15, %r11         // r11 = r4*term - (r4*term - (r4*term % p4)) = r4*term % p4
   
   mov     %r8, r1            // r1 = r1*term % p1
   mov     %r9, r2            // r2 = r2*term % p2
   mov     %r10, r3           // r3 = r3*term % p3
   mov     %r11, r4           // r4 = r4*term % p4

   // In some cases the remainder is p, so subtract it out
   sub     p1, %r8           // if r1 > p1, then r1 -= p1
   jl      nomove0            
   mov     %r8, r1          
                              
nomove0:                      
   sub     p2, %r9           // if r2 > p2, then r2 -= p2
   jl      nomove1            
   mov     %r9, r2         
                              
nomove1:                      
   sub     p3, %r10           // if r3 > p3, then r3 -= p3
   jl      nomove2            
   mov     %r10, r3      
                              
nomove2:                      
   sub     p4, %r11           // if r4 > p4, then r4 -= p4
   jl      nomove3
   mov     %r11, r3
           
nomove3:   
   // If any remainder is 1, then we will exit early
   cmp     $1, %r8
   je      all_done
   cmp     $1, %r9
   je      all_done
   cmp     $1, %r10
   je      all_done
   cmp     $1, %r11
   je      all_done

   mov     r1, %r8
   mov     r2, %r9
   mov     r3, %r10
   mov     r4, %r11
   
   // If any remainder is p-1, then we will exit early
   cmp     p1m1, %r8
   je      all_done
   cmp     p2m1, %r9
   je      all_done
   cmp     p3m1, %r10
   je      all_done
   cmp     p4m1, %r11
   je      all_done
           
next_term:  
   // term = term + mf
   add     mf, %rbp
   
   // If we have not reached the max term, continue
   cmp     max, %rbp
   jle     mf_loop

   // no factors found
   mov     $0, %rax

all_done:
   fstp    %st(0)
   fstp    %st(0)
   fstp    %st(0)
   fstp    %st(0)

#ifdef _WIN64
   add     $160, %rsp
   
   pop     %rdi
   pop     %rsi
#endif

   pop     %r15
   pop     %r14
   pop     %r13
   pop     %r12
   pop     %rbx
   pop     %rbp
   
   ret
