/* mfsieve.S -- Mark Rodenkirch

   int mfsieve(const uint32_t start, const uint32_t mf, const uint32_t *minmax, const uint64_t *P) __attribute__ ((pure));

   Find a factor for the smallest multifactorial in the series.
   P is an array of 4 primes.
   The value of start must be less than the smallest P in the array.

   Assumes p < 2^51 for all p in P.
   (2^25 < n < 2^31 allowed if p > n).

   Here is how it works:
      term = start
      while (term < max) {
         term += mf
         term *= (term % p)
         if (term = 1 or term = -1)
            return 1
      }
      
      return 0
*/

   .text

#if defined(__APPLE__)
# define FUN_NAME _mfsieve
#else
# define FUN_NAME mfsieve
#endif

#if _WIN64
   #define old_mxcsr 80(%rsp)
   #define new_mxcsr 84(%rsp)
   #define mf   0(%rsp)
   #define min  8(%rsp)
   #define max  16(%rsp)
   #define p1m1 24(%rsp)
   #define p2m1 32(%rsp)
   #define p3m1 40(%rsp)
   #define p4m1 48(%rsp)
#else
   #define old_mxcsr -4(%rsp)
   #define new_mxcsr -8(%rsp)
   #define mf   -24(%rsp)
   #define min  -32(%rsp)
   #define max  -40(%rsp)
   #define p1m1 -48(%rsp)
   #define p2m1 -56(%rsp)
   #define p3m1 -64(%rsp)
   #define p4m1 -72(%rsp)
#endif

   .p2align 4,,15
   .globl   FUN_NAME

FUN_NAME:
   push   %rbp
   push   %rbx
   push   %r12
   push   %r13
   push   %r14
   push   %r15

#if _WIN64
   push   %rsi
   push   %rdi
   
   sub    $120, %rsp

   mov    %rcx, %rdi
   mov    %rdx, mf
   mov    %r8, max
   mov    %r9, %rsi
   
   movdqa %xmm6, 32(%rsp)
   movdqa %xmm7, 48(%rsp)
   movdqa %xmm8, 64(%rsp)
#else
   mov    %rsi, mf
   mov    %rdx, max
   mov    %rcx, %rsi
#endif

   // %rdi = start
   // %rsi = P

   stmxcsr old_mxcsr
   mov     old_mxcsr, %eax
   or      $0x6000, %eax      /* Round to zero */
   mov     %eax, new_mxcsr
   ldmxcsr new_mxcsr
   
   // Now grab the primes for the mod operation
   mov   0(%rsi), %r8
   mov   8(%rsi), %r9
   mov   16(%rsi), %r10
   mov   24(%rsi), %r11
   
   mov   $1, %ebp

   cvtsi2sd %ebp, %xmm5
   cvtsi2sd %ebp, %xmm6
   cvtsi2sd %ebp, %xmm7
   cvtsi2sd %ebp, %xmm8

   cvtsi2sd %r8, %xmm1
   cvtsi2sd %r9, %xmm2
   cvtsi2sd %r10, %xmm3
   cvtsi2sd %r11, %xmm4

   // Now we have 1/p for the 4 input values
   divsd   %xmm1, %xmm5
   divsd   %xmm2, %xmm6
   divsd   %xmm3, %xmm7
   divsd   %xmm4, %xmm8
   
   sub   $1, %r8
   sub   $1, %r9
   sub   $1, %r10
   sub   $1, %r11
   
   // Now we have p-1 for the 4 input values
   mov   %r8, p1m1
   mov   %r9, p2m1
   mov   %r10, p3m1
   mov   %r11, p4m1
   
   add   $1, %r8
   add   $1, %r9
   add   $1, %r10
   add   $1, %r11
   
   // Starting term
   mov   %rdi, %rax
   mov   %rdi, %rbx
   mov   %rdi, %rcx
   mov   %rdi, %rdx
   
   mov   %rdi, %rbp
   jmp   nextterm

   .p2align 4,,15
loop4:
   // Convert rem to double
   cvtsi2sd %rax, %xmm1
   cvtsi2sd %rbx, %xmm2
   cvtsi2sd %rcx, %xmm3
   cvtsi2sd %rdx, %xmm4
   
   // Convert the term to double
   cvtsi2sd %rbp, %xmm11

   // Compute (as integer) rem = rem * term
   imul    %rbp, %rax
   imul    %rbp, %rbx
   imul    %rbp, %rcx
   imul    %rbp, %rdx

   // Compute (as double) rem = rem * term
   mulsd   %xmm11, %xmm1
   mulsd   %xmm11, %xmm2
   mulsd   %xmm11, %xmm3
   mulsd   %xmm11, %xmm4
   
   // Compute rem * 1/p
   mulsd   %xmm5, %xmm1
   mulsd   %xmm6, %xmm2
   mulsd   %xmm7, %xmm3
   mulsd   %xmm8, %xmm4

   // Convert back to integer which will effectively 
   // be a floor of (rem * term) / p
   cvtsd2si %xmm1, %r12
   cvtsd2si %xmm2, %r13
   cvtsd2si %xmm3, %r14
   cvtsd2si %xmm4, %r15

   // Now we can get the largest value of x where
   // x <= (rem * term) and x % p = 0
   imul  %r8, %r12
   imul  %r9, %r13
   imul  %r10, %r14
   imul  %r11, %r15

   // rem * term - x is our remainder
   sub   %r12, %rax
   sub   %r13, %rbx
   sub   %r14, %rcx
   sub   %r15, %rdx

   // Now we have modded by p and have our remainder
   mov   %rax, %r12
   mov   %rbx, %r13
   mov   %rcx, %r14
   mov   %rdx, %r15

   // But in some cases the remainder is p, so subtract it out
   sub   %r8, %r12
   jl    nomove0
   mov   %r12, %rax
   
nomove0:
   sub   %r9, %r13
   jl    nomove1
   mov   %r13, %rbx
   
nomove1:
   sub   %r10, %r14
   jl    nomove2
   mov   %r14, %rcx
   
nomove2:
   sub   %r11, %r15
   jl    nomove3
   mov   %r15, %rdx

nomove3: 
   // If any remainder is 1, then we will exit early
   cmp   $1, %rax
   je    earlyout
   cmp   $1, %rbx
   je    earlyout
   cmp   $1, %rcx
   je    earlyout
   cmp   $1, %rdx
   je    earlyout

   // If any remainder is p-1, then we will exit early
   cmp   p1m1, %rax
   je    earlyout
   cmp   p2m1, %rbx
   je    earlyout
   cmp   p3m1, %rcx
   je    earlyout
   cmp   p4m1, %rdx
   je    earlyout
   
nextterm:
   // term = term + mf
   add   mf, %rbp

   // If we have not reached the max term, continue
   cmp   max, %rbp
   jle   loop4

   mov   $0, %rax
   jmp   alldone
   
earlyout:
   mov   $1, %rax

alldone:
   ldmxcsr   old_mxcsr

#if _WIN64
   movdqa   32(%rsp), %xmm6
   movdqa   48(%rsp), %xmm7
   movdqa   64(%rsp), %xmm8
   add   $120, %rsp
   pop   %rdi
   pop   %rsi
#endif

   pop   %r15
   pop   %r14
   pop   %r13
   pop   %r12
   pop   %rbx
   pop   %rbp

   ret
