/* primorial.S -- Geoffrey Reynolds, October 2008.

   int primorial(const uint32_t *N, const uint64_t *P);

   N is an array of consecutive odd primes 3,..,max,0
   P is an array of 4 primes.

   Assumes n < 2^25 for all n in N.
   Assumes p < 2^51 for all p in P.
   (2^25 < n < 2^31 allowed if p > n).

   Returns the first n in N such that n#+/-1 is divisible by some p in P.
   Returns zero if there is no such n.
*/

   .text

   .p2align 4,,15
   .globl   _primorial
   .globl   primorial

_primorial:
primorial:

#if _WIN64
#define ARG1     %rcx
#define ARG2     %rdx
#define ARG3     %r8
#define nmax  88(%rsp)
#define M0     0(%rsp)
#define M1     8(%rsp)
#define M2    16(%rsp)
#define M3    24(%rsp)
#else
#define ARG1     %rdi
#define ARG2     %rsi
#define ARG3     %rdx
#define M0   -16(%rsp)
#define M1   -24(%rsp)
#define M2   -32(%rsp)
#define M3   -40(%rsp)
#define nmax -44(%rsp)
#endif

   push     %rbp
   push     %rbx
   push     %r12
   push     %r13
   push     %r14
   push     %r15

#if _WIN64
   push     %rsi
   push     %rdi
   sub      $104, %rsp
   
   mov      %rcx, %rdi
   
   movdqa   %xmm6, 32(%rsp)
   movdqa   %xmm7, 48(%rsp)
   movdqa   %xmm8, 64(%rsp)
#endif

   mov      ARG3, nmax

   mov      0(ARG2), %r8
   mov      8(ARG2), %r9
   mov      16(ARG2), %r10
   mov      24(ARG2), %r11

   mov      $1, %ebp

   cvtsi2sd %ebp, %xmm5
   cvtsi2sd %ebp, %xmm6
   cvtsi2sd %ebp, %xmm7
   cvtsi2sd %ebp, %xmm8

   cvtsi2sd %r8, %xmm1
   cvtsi2sd %r9, %xmm2
   cvtsi2sd %r10, %xmm3
   cvtsi2sd %r11, %xmm4

   divsd    %xmm1, %xmm5
   divsd    %xmm2, %xmm6
   divsd    %xmm3, %xmm7
   divsd    %xmm4, %xmm8

   lea      -1(%r8), %r12
   lea      -1(%r9), %r13
   lea      -1(%r10), %r14
   lea      -1(%r11), %r15

   mov      %r12, M0
   mov      %r13, M1
   mov      %r14, M2
   mov      %r15, M3

   mov      $1, %ebp
   mov      %ebp, %ebx
   mov      %ebp, %ecx
   mov      %ebp, %edx

   jmp      test4

loop4:
   cvtsi2sd %eax, %xmm0

   cvtsi2sd %rbp, %xmm1
   cvtsi2sd %rbx, %xmm2
   cvtsi2sd %rcx, %xmm3
   cvtsi2sd %rdx, %xmm4

   imul     %rax, %rbp
   imul     %rax, %rbx
   imul     %rax, %rcx
   imul     %rax, %rdx

   mulsd    %xmm0, %xmm1
   mulsd    %xmm0, %xmm2
   mulsd    %xmm0, %xmm3
   mulsd    %xmm0, %xmm4
   
   mulsd    %xmm5, %xmm1
   mulsd    %xmm6, %xmm2
   mulsd    %xmm7, %xmm3
   mulsd    %xmm8, %xmm4

   cvtsd2si %xmm1, %r12
   cvtsd2si %xmm2, %r13
   cvtsd2si %xmm3, %r14
   cvtsd2si %xmm4, %r15

   imul     %r8, %r12
   imul     %r9, %r13
   imul     %r10, %r14
   imul     %r11, %r15
   
   sub      %r12, %rbp
   sub      %r13, %rbx
   sub      %r14, %rcx
   sub      %r15, %rdx
   
   mov      %rbp, %r12
   mov      %rbx, %r13
   mov      %rcx, %r14
   mov      %rdx, %r15

   sub      %r8, %r12
   jl       p1
   mov      %r12, %rbp
   
p1:   
   sub      %r9, %r13
   jl       p2
   mov      %r13, %rbx
   
p2:
   sub      %r10, %r14
   jl       p3
   mov      %r14, %rcx
   
p3:
   sub      %r11, %r15
   jl       p4
   mov      %r15, %rdx

p4:
   cmp      $1, %rbp
   je       done4
   cmp      M0, %rbp
   je       done4
   cmp      $1, %rbx
   je       done4
   cmp      M1, %rbx
   je       done4
   cmp      $1, %rcx
   je       done4
   cmp      M2, %rcx
   je       done4
   cmp      $1, %rdx
   je       done4
   cmp      M3, %rdx
   je       done4
   
test4:   
   mov      0(%rdi), %eax
   cmp      $0, %eax
   lea      4(%rdi), %rdi
   jg       loop4
   
   xor      %eax, %eax

done4:

#if _WIN64
   movdqa   32(%rsp), %xmm6
   movdqa   48(%rsp), %xmm7
   movdqa   64(%rsp), %xmm8
   add      $104, %rsp
   pop      %rdi
   pop      %rsi
#endif

   pop      %r15
   pop      %r14
   pop      %r13
   pop      %r12
   pop      %rbx
   pop      %rbp

   ret
