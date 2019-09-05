/* pixsieve.S -- Mark Rodenkirch

   int pixsieve(const uint32_t *N1, uint32_t *N2, const uint64_t *P, uint64_t mult);

   N1:N2 are two arrays of integers where end-of-list is equal to mult.
   P is an array of 4 primes.

   Assumes p < 2^51 for all p in P.
   (2^25 < n < 2^31 allowed if p > n).

   Returns 1 if we find a term that is divisible by some p in P.
   Returns zero if there is no such n.
*/

   .text

#if defined(__APPLE__)
# define FUN_NAME _pixsieve
#else
# define FUN_NAME pixsieve
#endif

#if _WIN64
   #define old_mxcsr 80(%rsp)
   #define new_mxcsr 84(%rsp)
   #define mult  0(%rsp)
   #define n2    8(%rsp)
#else
   #define old_mxcsr -4(%rsp)
   #define new_mxcsr -8(%rsp)
   #define mult -24(%rsp)
   #define n1   -32(%rsp)
   #define n2   -40(%rsp)
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
   mov    %rdx, n2
   mov    %r8,  %rsi
   mov    %r9,  mult
   
   movdqa %xmm6, 32(%rsp)
   movdqa %xmm7, 48(%rsp)
   movdqa %xmm8, 64(%rsp)
#else
   mov    %rdi, n1
   mov    %rsi, n2
   mov    %rdx, %rsi
   mov    %rcx, mult
#endif

   // %rdi = N
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
   
   mov   $0, %rbp
   mov   $0, %rbx
   mov   $0, %rcx
   mov   $0, %rdx
   
   mov   $0, %rax
   jmp   nextterm

   .p2align 4,,15
loop4:
#if PREFETCH
   prefetchnta 4*(PREFETCH-1)(%rdi)
#endif

   cvtsi2sd %rbp, %xmm1
   cvtsi2sd %rbx, %xmm2
   cvtsi2sd %rcx, %xmm3
   cvtsi2sd %rdx, %xmm4

   mov     mult, %r12

   // Mulitple the remainder by the multiplier, which is a power of 10
   imul    %r12, %rbp
   imul    %r12, %rbx
   imul    %r12, %rcx
   imul    %r12, %rdx

   // Add our next term
   add     %rax, %rbp
   add     %rax, %rbx
   add     %rax, %rcx
   add     %rax, %rdx

   cvtsi2sd %r12, %xmm10
   cvtsi2sd %eax, %xmm11

   mulsd   %xmm10, %xmm1
   mulsd   %xmm10, %xmm2
   mulsd   %xmm10, %xmm3
   mulsd   %xmm10, %xmm4

   addsd   %xmm11, %xmm1
   addsd   %xmm11, %xmm2
   addsd   %xmm11, %xmm3
   addsd   %xmm11, %xmm4

   mulsd   %xmm5, %xmm1
   mulsd   %xmm6, %xmm2
   mulsd   %xmm7, %xmm3
   mulsd   %xmm8, %xmm4

   cvtsd2si %xmm1, %r12
   cvtsd2si %xmm2, %r13
   cvtsd2si %xmm3, %r14
   cvtsd2si %xmm4, %r15

   imul   %r8, %r12
   imul   %r9, %r13
   imul   %r10, %r14
   imul   %r11, %r15

   sub   %r12, %rbp
   sub   %r13, %rbx
   sub   %r14, %rcx
   sub   %r15, %rdx

   // Now we have modded by p and have our remainder
   mov   %rbp, %r12
   mov   %rbx, %r13
   mov   %rcx, %r14
   mov   %rdx, %r15

   // But in some cases the remainder is p, so subtract it out
   sub   %r8, %r12
   jl    nomove0
   mov   %r12, %rbp
   
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
   // Only when our multiplier is 10, do we want to compare the remainder to 0
   cmpl  $10, mult
   jne   nextterm

   // If any remainder is 0, then we will exit early
   cmp   $0, %rbp
   je    earlyout
   cmp   $0, %rbx
   je    earlyout
   cmp   $0, %rcx
   je    earlyout
   cmp   $0, %rdx
   je    earlyout

nextterm:
   mov   0(%rdi), %eax
   
   // If we have not exhausted the list, then multiply by the next term
   cmp   mult, %rax
   lea   4(%rdi), %rdi
   jl    loop4
  
   xor   %eax, %eax

   // If mult = 10, then we have exhausted the second list
   cmpl  $10, mult
   je    alldone

   // Otherwise set it to 10 and point to the second list
   mov   $10, %rdi
   mov   %rdi, mult
   mov   n2, %rdi
   jmp   nextterm

earlyout:
   mov   $1, %eax

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
