/* mulmod192.S -- (C) Geoffrey Reynolds, March 2008.

   These routines are based on Pierrick Gaudrys REDC code in GMP-ECM.


   void mulmod192(uint64_t *R, const uint64_t *A, const uint64_t *B,
                  const uint64_t *N, uint64_t inv);

   Assign R <-- A*B (mod N) in Montgomery form, where 0 <= A,B < N < 2^192.
   Assumes that A,B are in Montgomery form, N is odd, and N[0]*inv = -1.


   void mulmod192_proth0(uint64_t *R, const uint64_t *A, const uint64_t *B,
                        const uint64_t *N);

   As above, but assumes that N[2] < 2^63, N[0] = 1.


   void mulmod192_proth1(uint64_t *R, const uint64_t *A, const uint64_t *B,
                         const uint64_t *N);

   As above, but assumes that N[1] = 0, N[0] = 1.


   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

	.globl	_mulmod192
	.globl	mulmod192
	.p2align 4,,15

_mulmod192:
mulmod192:
#ifdef _WIN64
	push	%rdi
	push	%rsi
	mov	%rcx, %rdi
	mov	%rdx, %rsi
	mov	%r8, %rdx
	mov	%r9, %rcx
	mov	64(%rsp), %r8
#endif
	mov	(%rdx), %rax
	mov	8(%rdx), %r10
	mov	16(%rdx), %rdx
	mov	%rcx, %r9

	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14

	mov	%rax, (%rdi)
	mov	%r10, 8(%rdi)
	mov	%rdx, 16(%rdi)

	/* %r8  = inv
	   %r9  = N
	   %rsi = A
	   %rdi = B,R */

	mov	(%rsi), %rax
	mulq	(%rdi)
	mul	%r8
	mov	%rax, %r10

	mov	(%r9), %rax
	mul	%r10
	mov	%rax, %r12
	mov	%rdx, %r14

	mov	8(%r9), %rax
	mul	%r10
	add	%rax, %r14
	mov	%rdx, %rbp
	adc	$0, %rbp

	mov	16(%r9), %rax
	mul	%r10
	add	%rax, %rbp
	mov	%rdx, %r11
	adc	$0, %r11

	mov	(%rsi), %r10
	mov	(%rdi), %rax
	mul	%r10
	mov	%rax, %rbx
	mov	%rdx, %rcx

	mov	8(%rdi), %rax
	mul	%r10
	add	%rbx, %r12
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx

	mov	16(%rdi), %rax
	mul	%r10
	add	%rcx, %r14
	adc	%rax, %rbx
	adc	$0, %rdx
	xor	%r12, %r12
	add	%rbx, %rbp
	adc	%rdx, %r11
	adc	$0, %r12


	mov	8(%rsi), %rax
	mulq	(%rdi)
	add	%r14, %rax
	mul	%r8
	mov	%rax, %r10

	mov	(%r9), %rax
	mul	%r10
	mov	%rax, %rbx
	mov	%rdx, %rcx

	mov	8(%r9), %rax
	mul	%r10
	add	%rbx, %r14
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx

	mov	16(%r9), %rax
	mul	%r10
	add	%rcx, %rbp
	adc	%rax, %rbx
	adc	$0, %rdx
	xor	%r13, %r13
	add	%rbx, %r11
	adc	%rdx, %r12
	adc	$0, %r13

	mov	8(%rsi), %r10
	mov	(%rdi), %rax
	mul	%r10
	mov	%rax, %rbx
	mov	%rdx, %rcx

	mov	8(%rdi), %rax
	mul	%r10
	add	%rbx, %r14
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx

	mov	16(%rdi), %rax
	mul	%r10
	add	%rcx, %rbp
	adc	%rax, %rbx
	adc	$0, %rdx
	add	%rbx, %r11
	adc	%rdx, %r12
	adc	$0, %r13


	mov	16(%rsi), %rax
	mulq	(%rdi)
	add	%rbp, %rax
	mul	%r8
	mov	%rax, %r10

	mov	(%r9), %rax
	mul	%r10
	mov	%rax, %rbx
	mov	%rdx, %rcx

	mov	8(%r9), %rax
	mul	%r10
	add	%rbx, %rbp
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx

	mov	16(%r9), %rax
	mul	%r10
	add	%rcx, %r11
	adc	%rax, %rbx
	adc	$0, %rdx
	xor	%r8, %r8
	add	%rbx, %r12
	adc	%rdx, %r13
	adc	$0, %r8

	mov	16(%rsi), %r10
	mov	(%rdi), %rax
	mul	%r10
	mov	%rax, %rbx
	mov	%rdx, %rcx

	mov	8(%rdi), %rax
	mul	%r10
	add	%rbx, %rbp
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx

	mov	16(%rdi), %rax
	mul	%r10
	add	%rcx, %r11
	adc	%rax, %rbx
	adc	$0, %rdx
	add	%rbx, %r12
	adc	%rdx, %r13
	adc	$0, %r8

0:	mov	%r11, (%rdi)
	mov	%r12, 8(%rdi)
	mov	%r13, 16(%rdi)
	sub	(%r9), %r11
	sbb	8(%r9), %r12
	sbb	16(%r9), %r13
	sbb	$0, %r8
	jnc	0b

	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx
#ifdef _WIN64
	pop	%rsi
	pop	%rdi
#endif
	ret
   


	.text
	.globl	_mulmod192_proth0
	.globl	mulmod192_proth0
	.p2align 4,,15

_mulmod192_proth0:
mulmod192_proth0:
#ifdef _WIN64
	push	%rdi
	push	%rsi
	mov	%rcx, %rdi
	mov	%rdx, %rsi
	mov	%r8, %rdx
	mov	%r9, %rcx
#endif
	mov	(%rdx), %rax
	mov	8(%rdx), %r10
	mov	16(%rdx), %rdx

	mov	%rcx, %r9
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13

	mov	%rax, (%rdi)
	mov	%r10, 8(%rdi)
	mov	%rdx, 16(%rdi)

	/* %r9  = N
	   %rsi = A
	   %rdi = B,R  */

	mov	(%rsi), %rax
	mulq	(%rdi)
	neg	%rax
	mov	%rax, %r12

	mov	8(%r9), %rax
	mul	%r12
	mov	%rax, %r8
	mov	%rdx, %rbp

	mov	16(%r9), %rax	/* < 2^63 */
	mul	%r12
	xor	%r11, %r11
	add	%rax, %rbp
	adc	%rdx, %r11	/* < 2^63 */

	mov	(%rsi), %r10
	mov	(%rdi), %rax
	mul	%r10
	mov	%rax, %rbx
	mov	%rdx, %rcx

	mov	8(%rdi), %rax
	mul	%r10
	add	%rbx, %r12
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx

	mov	16(%rdi), %rax	/* < 2^63-1 */
	mul	%r10
	add	%rcx, %r8
	adc	%rax, %rbx
	adc	$0, %rdx	/* < 2^63-1 */
	add	%rbx, %rbp
	adc	%rdx, %r11	/* < 2^63, so CF=0 */


	mov	8(%rsi), %rax
	mulq	(%rdi)
	add	%r8, %rax
	neg	%rax
	mov	%rax, %r10

	mov	8(%r9), %rax
	mul	%r10
	xor	%ecx, %ecx
	xor	%ebx, %ebx
	add	%r10, %r8
	adc	%rax, %rcx
	adc	%rdx, %rbx

	mov	16(%r9), %rax	/* < 2^63 */
	mul	%r10
	add	%rcx, %rbp
	adc	%rax, %rbx
	adc	$0, %rdx	/* < 2^63 */
	xor	%r12, %r12
	add	%rbx, %r11
	adc	%rdx, %r12	/* <= 2^63, so CF=0 */

	mov	8(%rsi), %r10
	mov	(%rdi), %rax
	mul	%r10
	mov	%rax, %rbx
	mov	%rdx, %rcx

	mov	8(%rdi), %rax
	mul	%r10
	add	%rbx, %r8
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx

	mov	16(%rdi), %rax	/* < 2^63-1 */
	mul	%r10
	add	%rcx, %rbp
	adc	%rax, %rbx
	adc	$0, %rdx	/* < 2^63-1 */
	add	%rbx, %r11
	adc	%rdx, %r12	/* CF=0 */


	mov	16(%rsi), %rax	/* < 2^63-1 */
	mulq	(%rdi)
	add	%rbp, %rax
	neg	%rax
	mov	%rax, %r10

	mov	8(%r9), %rax
	mul	%r10
	xor	%ecx, %ecx
	xor	%ebx, %ebx
	add	%r10, %rbp
	adc	%rax, %rcx
	adc	%rdx, %rbx

	mov	16(%r9), %rax	/* < 2^63 */
	mul	%r10
	add	%rcx, %r11
	adc	%rax, %rbx
	adc	$0, %rdx	/* < 2^63 */
	xor	%r13, %r13
	add	%rbx, %r12
	adc	%rdx, %r13	/* <= 2^63, so CF=0 */

	mov	16(%rsi), %r10	/* < 2^63-1 */
	mov	(%rdi), %rax
	mul	%r10
	mov	%rax, %rbx
	mov	%rdx, %rcx	/* < 2^63-2 */

	mov	8(%rdi), %rax
	mul	%r10
	add	%rbx, %rbp
	adc	%rax, %rcx
	mov	%rdx, %rbx	/* < 2^63-2 */
	adc	$0, %rbx

	mov	16(%rdi), %rax	/* < 2^63-1 */
	mul	%r10
	add	%r11, %rcx
	adc	%rbx, %rax
	adc	$0, %rdx	/* < 2^62-1 */
	add	%r12, %rax
	adc	%r13, %rdx	/* CF=0 */

0:	mov	%rcx, (%rdi)
	mov	%rax, 8(%rdi)
	mov	%rdx, 16(%rdi)
	sub	$1, %rcx
	sbb	8(%r9), %rax
	sbb	16(%r9), %rdx
	jnc	0b

	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx
#ifdef _WIN64
	pop	%rsi
	pop	%rdi
#endif
	ret




	.text
	.globl	_mulmod192_proth1
	.globl	mulmod192_proth1
	.p2align 4,,15

_mulmod192_proth1:
mulmod192_proth1:
#ifdef _WIN64
	push	%rdi
	push	%rsi
	mov	%rcx, %rdi
	mov	%rdx, %rsi
	mov	%r8, %rdx
	mov	%r9, %rcx
#endif
	mov	(%rdx), %rax
	mov	8(%rdx), %r11
	mov	16(%rdx), %rdx
	mov	16(%rcx), %r9

	mov	(%rsi), %rcx
	mov	8(%rsi), %r10
	mov	16(%rsi), %r8

	push	%rbx
	push	%rbp
	push	%r12
	push	%r13

	mov	%rax, (%rdi)
	mov	%r11, 8(%rdi)
	mov	%rdx, 16(%rdi)

	/* %r9  = N[2]
	   %rcx = A[0]
	   %r10 = A[1]
	   %r8  = A[2]
	   %rdi = B,R */


	mul	%rcx
	mov	%rax, %rbx
	mov	%rdx, %rsi
	neg	%rax
	mov	%rax, %r12

	mul	%r9
	mov	%rax, %rbp
	mov	%rdx, %r11

	mov	8(%rdi), %rax
	mul	%rcx
	xor	%r13, %r13
	add	%rbx, %r12	/* set carry if nonzero */
	adc	%rax, %rsi
	adc	%rdx, %r13

	mov	16(%rdi), %rax
	mul	%rcx
	add	%rax, %r13
	adc	$0, %rdx
	xor	%r12, %r12
	add	%r13, %rbp
	adc	%rdx, %r11
	adc	$0, %r12


	mov	(%rdi), %rax
	mul	%r10
	add	%rsi, %rax
	neg	%rax

	xor	%ecx, %ecx
	add	%rax, %rsi
	adc	$0, %rcx

	mul	%r9
	add	%rcx, %rbp
	adc	%rax, %r11
	adc	%rdx, %r12	/* CF=0 */

	mov	(%rdi), %rax
	mul	%r10
	mov	%rax, %rbx
	mov	%rdx, %rcx

	mov	8(%rdi), %rax
	mul	%r10
	add	%rbx, %rsi
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx

	mov	16(%rdi), %rax
	mul	%r10
	add	%rcx, %rbp
	adc	%rax, %rbx
	adc	$0, %rdx
	xor	%r13, %r13
	add	%rbx, %r11
	adc	%rdx, %r12
	adc	$0, %r13


	mov	(%rdi), %rax
	mul	%r8
	add	%rbp, %rax
	neg	%rax

	xor	%ecx, %ecx
	add	%rax, %rbp
	adc	$0, %rcx

	mul	%r9
	xor	%ebx, %ebx
	add	%rcx, %r11
	adc	%rax, %rbx
	adc	$0, %rdx
	add	%rbx, %r12
	adc	%rdx, %r13	/* CF=0 */

	mov	(%rdi), %rax
	mul	%r8
	mov	%rax, %rbx
	mov	%rdx, %rcx

	mov	8(%rdi), %rax
	mul	%r8
	add	%rbx, %rbp
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx

	mov	16(%rdi), %rax
	mul	%r8
	add	%r11, %rcx
	adc	%rbx, %rax
	adc	$0, %rdx
	xor	%esi, %esi
	add	%r12, %rax
	adc	%r13, %rdx
	adc	$0, %esi

0:	mov	%rcx, (%rdi)
	mov	%rax, 8(%rdi)
	mov	%rdx, 16(%rdi)
	sub	$1, %rcx
	sbb	$0, %rax
	sbb	%r9, %rdx
	sbb	$0, %esi
	jnc	0b

	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx
#ifdef _WIN64
	pop	%rsi
	pop	%rdi
#endif
	ret
