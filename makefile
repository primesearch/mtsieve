# Makefile for mtsieve applications. (C) Mark Rodenkirch
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.

# Set DEBUG=yes to compile with debugging information and internal checks.

# On Windows, the OpenCL headers are cloned from:
#    git clone --recursive https://github.com/KhronosGroup/OpenCL-SDK.git

# On Windows, I am using the llvm-mingw compiler.  I have had problems trying
# to use other compilers on Windows.  This one seems to work.

# Metal is only available on Mac OS X

DEBUG=no
CC=g++
PERL=perl

PERL_FLAGS=cltoh.pl

CPP_FLAGS=-Isieve -m64 -Wall
CPP_FLAGS_OPENCL=-DUSE_OPENCL
CPP_FLAGS_METAL=-DUSE_METAL
CPP_FLAGS_GMP=
CPP_FLAGS_SIEVE=-DNDEBUG

LD_FLAGS=
LD_FLAGS_STDC=-lstdc++
LD_FLAGS_OPENCL=-lstdc++
LD_FLAGS_METAL=
LD_FLAGS_GMP=-lgmp

GPULIBS=

ifeq ($(strip $(DEBUG)),yes)
   OPT_CPP_FLAGS_02=-g
   OPT_CPP_FLAGS=-g
else
   OPT_CPP_FLAGS_02=-O2
   OPT_CPP_FLAGS=-O3
endif

ifeq ($(OS),Windows_NT)
   HAS_X86=yes
   HAS_ARM=no
   LD_FLAGS+=-static
   CPP_FLAGS+=-I"E:\stuff\gmp-6.1.2" -DUSE_X86
   CPP_FLAGS_OPENCL+=-I"E:\stuff\OpenCL-SDK\external\OpenCL-Headers" -DCL_TARGET_OPENCL_VERSION=300
   LD_FLAGS_GMP+=-L"E:\stuff\gmp-6.1.2\.libs"
   LD_FLAGS_OPENCL="C:\Windows\System32\OpenCl.dll"
else
   UNAME_S := $(shell uname -s)
   UNAME_A := $(shell uname -a)
   
   ifeq ($(findstring ARM,$(UNAME_A)),ARM)
      HAS_X86=no
      HAS_ARM=yes
      CPP_FLAGS+=-DUSE_ARM
   else
      HAS_X86=yes
      HAS_ARM=no
      CPP_FLAGS+=-DUSE_X86
   endif
   
   ifeq ($(UNAME_S),Darwin)
      XC=xcrun
      XC_FLAGS=-sdk macosx metal
      HAS_METAL=yes
      CPP_FLAGS+=-std=c++17
      CPP_FLAGS_METAL=-DUSE_METAL -Igpu_metal
      LD_FLAGS_OPENCL+=-framework OpenCL
      LD_FLAGS_METAL+=-framework Metal -framework Foundation -framework CoreGraphics
   else
      CPP_FLAGS+=-std=c++11
      LD_FLAGS+=-lpthread
      LD_FLAGS_GMP=-lgmp
      LD_FLAGS_OPENCL+=-lOpenCL
   endif
endif
   
# No changes should be needed below here.

ifeq ($(strip $(HAS_X86)),yes)
CPU_PROGS=afsieve ccsieve cksieve dmdsieve fbncsieve fkbnsieve gcwsieve gfndsieve hcwsieve k1b2sieve kbbsieve lifsieve \
   mfsieve pixsieve psieve sgsieve smsieve smwsieve srsieve2 twinsieve xyyxsieve
else
# no non-x86 builds for afsieve, gfndsieve, pixsieve, and xyyxsieve
CPU_PROGS=ccsieve cksieve dmdsieve fbncsieve fbncsieve fkbnsieve gcwsieve hcwsievek1b2sieve kbbsieve lifsieve \
   mfsieve psieve sgsieve smsieve smwsieve srsieve2 twinsieve
endif

ifeq ($(strip $(HAS_X86)),yes)
OPENCL_PROGS=afsievecl cksievecl gcwsievecl gfndsievecl hcwsievecl lifsievecl mfsievecl pixsievecl smsievecl smwsievecl srsieve2cl xyyxsievecl psievecl
else
OPENCL_PROGS=cksievecl cwsievecl gcwsievecl gfndsievecl hcwsievecl lifsievecl mfsievecl smsievecl smwsievecl srsieve2cl psievecl
endif

METAL_PROGS=cksievemtl cwsievemtl gcwsievemtl gfndsievemtl hcwsievemtl lifsievemtl mfsievemtl psievemtl smsievemtl srsieve2mtl

CPU_CORE_OBJS=core/App_cpu.o core/FactorApp_cpu.o core/AlgebraicFactorApp_cpu.o \
   core/Clock_cpu.o core/Parser_cpu.o core/Worker_cpu.o core/main_cpu.o core/SharedMemoryItem_cpu.o \
   core/HashTable_cpu.o core/BigHashTable_cpu.o core/SmallHashTable_cpu.o core/TinyHashTable_cpu.o 
   
OPENCL_CORE_OBJS=core/App_opencl.o core/FactorApp_opencl.o core/AlgebraicFactorApp_opencl.o core/GpuDevice_opencl.o core/GpuKernel_opencl.o \
   core/Clock_opencl.o core/Parser_opencl.o core/Worker_opencl.o core/main_opencl.o core/SharedMemoryItem_opencl.o \
   core/HashTable_opencl.o core/BigHashTable_opencl.o core/SmallHashTable_opencl.o core/TinyHashTable_opencl.o \
   gpu_opencl/OpenCLDevice_opencl.o gpu_opencl/OpenCLKernel_opencl.o gpu_opencl/OpenCLErrorChecker_opencl.o

METAL_CORE_OBJS=core/App_metal.o core/FactorApp_metal.o core/AlgebraicFactorApp_metal.o core/GpuDevice_metal.o core/GpuKernel_metal.o \
   core/Clock_metal.o core/Parser_metal.o core/Worker_metal.o core/main_metal.o core/SharedMemoryItem_metal.o \
   core/HashTable_metal.o core/BigHashTable_metal.o core/SmallHashTable_metal.o core/TinyHashTable_metal.o \
   gpu_metal/MetalDevice_metal.o gpu_metal/MetalKernel_metal.o

ifeq ($(strip $(HAS_X86)),yes)
   ASM_OBJS=x86_asm/fpu_mod_init_fini.o x86_asm/fpu_push_pop.o \
      x86_asm/fpu_mulmod.o x86_asm/fpu_powmod.o x86_asm/fpu_powmod_4b_1n_4p.o \
      x86_asm/fpu_mulmod_iter.o x86_asm/fpu_mulmod_iter_4a.o x86_asm/fpu_mulmod_4a_4b_4p.o \
      x86_asm/avx_set_a.o x86_asm/avx_set_b.o x86_asm/avx_get.o \
      x86_asm/avx_compute_reciprocal.o x86_asm/avx_compare.o \
      x86_asm/avx_mulmod.o x86_asm/avx_powmod.o

   ASM_EXT_OBJS=x86_asm_ext/m320.o x86_asm_ext/m384.o x86_asm_ext/m448.o x86_asm_ext/m512.o \
      x86_asm_ext/m576.o x86_asm_ext/m640.o x86_asm_ext/m704.o x86_asm_ext/m768.o \
      x86_asm_ext/mulmod128.o x86_asm_ext/mulmod192.o x86_asm_ext/mulmod256.o \
      x86_asm_ext/sqrmod128.o x86_asm_ext/sqrmod192.o x86_asm_ext/sqrmod256.o \
      x86_asm_ext/redc.o
endif

PRIMESIEVE_OBJS=sieve/Erat.o sieve/EratBig.o sieve/EratMedium.o sieve/EratSmall.o sieve/PreSieve.o \
   sieve/CpuInfo.o sieve/MemoryPool.o sieve/PrimeGenerator.o sieve/PrimeSieve.o \
   sieve/IteratorHelper.o sieve/LookupTables.o sieve/popcount.o sieve/nthPrime.o sieve/CountPrintPrimes.o \
   sieve/ParallelSieve.o sieve/iterator.o sieve/api.o sieve/SievingPrimes.o

AF_OBJS=alternating_factorial/AlternatingFactorialApp_cpu.o alternating_factorial/AlternatingFactorialWorker_cpu.o alternating_factorial/afsieve.o
CC_OBJS=cunningham_chain/CunninghamChainApp.o cunningham_chain/CunninghamChainWorker.o
CK_OBJS=carol_kynea/CarolKyneaApp.o carol_kynea/CarolKyneaWorker.o
DMD_OBJS=dm_divisor/DMDivisorApp.o dm_divisor/DMDivisorWorker.o
FBNC_OBJS=fixed_bnc/FixedBNCApp.o fixed_bnc/FixedBNCWorker.o
FKBN_OBJS=fixed_kbn/FixedKBNApp.o fixed_kbn/FixedKBNWorker.o
GCW_OBJS=cullen_woodall/CullenWoodallApp_cpu.o cullen_woodall/CullenWoodallWorker_cpu.o
GFND_OBJS=gfn_divisor/GFNDivisorApp_cpu.o gfn_divisor/GFNDivisorTester_cpu.o gfn_divisor/GFNDivisorWorker_cpu.o
HCW_OBJS=hyper_cullen_woodall/HyperCullenWoodallApp_cpu.o hyper_cullen_woodall/HyperCullenWoodallWorker_cpu.o
K1B2_OBJS=k1b2/K1B2App.o k1b2/K1B2Worker.o
KBB_OBJS=kbb/KBBApp.o kbb/KBBWorker.o
LIF_OBJS=lifchitz/LifchitzApp_cpu.o lifchitz/LifchitzWorker_cpu.o
MF_OBJS=multi_factorial/MultiFactorialApp_cpu.o multi_factorial/MultiFactorialWorker_cpu.o
PIX_OBJS=primes_in_x/PrimesInXApp_cpu.o primes_in_x/PrimesInXWorker_cpu.o primes_in_x/pixsieve.o
PRIM_OBJS=primorial/PrimorialApp_cpu.o primorial/PrimorialWorker_cpu.o
TWIN_OBJS=twin/TwinApp.o twin/TwinWorker.o
SG_OBJS=sophie_germain/SophieGermainApp.o sophie_germain/SophieGermainWorker.o
SM_OBJS=smarandache/SmarandacheApp_cpu.o smarandache/SmarandacheWorker_cpu.o
SMW_OBJS=smarandache_wellin/SmarandacheWellinApp_cpu.o smarandache_wellin/SmarandacheWellinWorker_cpu.o
SR2_OBJS=sierpinski_riesel/SierpinskiRieselApp_cpu.o sierpinski_riesel/AlgebraicFactorHelper_cpu.o \
   sierpinski_riesel/AbstractSequenceHelper_cpu.o sierpinski_riesel/AbstractWorker_cpu.o \
   sierpinski_riesel/GenericSequenceHelper_cpu.o sierpinski_riesel/GenericWorker_cpu.o \
   sierpinski_riesel/CisOneSequenceHelper_cpu.o sierpinski_riesel/CisOneWithOneSequenceHelper_cpu.o \
   sierpinski_riesel/CisOneWithOneSequenceWorker_cpu.o \
   sierpinski_riesel/CisOneWithMultipleSequencesHelper_cpu.o sierpinski_riesel/CisOneWithMultipleSequencesWorker_cpu.o
XYYX_OBJS=xyyx/XYYXApp_cpu.o xyyx/XYYXWorker_cpu.o xyyx/XYYXSparseWorker_cpu.o

AF_OPENCL_OBJS=alternating_factorial/AlternatingFactorialApp_opencl.o alternating_factorial/AlternatingFactorialWorker_opencl.o alternating_factorial/afsieve.o alternating_factorial/AlternatingFactorialGpuWorker_opencl.o
CK_OPENCL_OBJS=carol_kynea/CarolKyneaApp_opencl.o carol_kynea/CarolKyneaWorker_opencl.o carol_kynea/CarolKyneaGpuWorker_opencl.o
GCW_OPENCL_OBJS=cullen_woodall/CullenWoodallApp_opencl.o cullen_woodall/CullenWoodallWorker_opencl.o cullen_woodall/CullenWoodallGpuWorker_opencl.o
GFND_OPENCL_OBJS=gfn_divisor/GFNDivisorApp_opencl.o gfn_divisor/GFNDivisorTester_opencl.o gfn_divisor/GFNDivisorWorker_opencl.o gfn_divisor/GFNDivisorGpuWorker_opencl.o
HCW_OPENCL_OBJS=hyper_cullen_woodall/HyperCullenWoodallApp_opencl.o hyper_cullen_woodall/HyperCullenWoodallWorker_opencl.o hyper_cullen_woodall/HyperCullenWoodallGpuWorker_opencl.o
LIF_OPENCL_OBJS=lifchitz/LifchitzApp_opencl.o lifchitz/LifchitzWorker_opencl.o lifchitz/LifchitzGpuWorker_opencl.o
MF_OPENCL_OBJS=multi_factorial/MultiFactorialApp_opencl.o multi_factorial/MultiFactorialWorker_opencl.o multi_factorial/MultiFactorialGpuWorker_opencl.o
PIX_OPENCL_OBJS=primes_in_x/PrimesInXApp_opencl.o primes_in_x/PrimesInXWorker_opencl.o primes_in_x/pixsieve.o primes_in_x/PrimesInXGpuWorker_opencl.o
PRIM_OPENCL_OBJS=primorial/PrimorialApp_opencl.o primorial/PrimorialWorker_opencl.o primorial/PrimorialGpuWorker_opencl.o
SM_OPENCL_OBJS=smarandache/SmarandacheApp_opencl.o smarandache/SmarandacheWorker_opencl.o smarandache/SmarandacheGpuWorker_opencl.o
SMW_OPENCL_OBJS=smarandache_wellin/SmarandacheWellinApp_opencl.o smarandache_wellin/SmarandacheWellinWorker_opencl.o smarandache_wellin/SmarandacheWellinGpuWorker_opencl.o
SR2_OPENCL_OBJS=sierpinski_riesel/SierpinskiRieselApp_opencl.o sierpinski_riesel/AlgebraicFactorHelper_opencl.o \
   sierpinski_riesel/AbstractSequenceHelper_opencl.o sierpinski_riesel/AbstractWorker_opencl.o \
   sierpinski_riesel/GenericSequenceHelper_opencl.o sierpinski_riesel/GenericWorker_opencl.o sierpinski_riesel/GenericGpuWorker_opencl.o \
   sierpinski_riesel/CisOneSequenceHelper_opencl.o sierpinski_riesel/CisOneWithOneSequenceHelper_opencl.o \
   sierpinski_riesel/CisOneWithOneSequenceWorker_opencl.o sierpinski_riesel/CisOneWithOneSequenceGpuWorker_opencl.o \
   sierpinski_riesel/CisOneWithMultipleSequencesHelper_opencl.o sierpinski_riesel/CisOneWithMultipleSequencesWorker_opencl.o \
   sierpinski_riesel/CisOneWithMultipleSequencesGpuWorker_opencl.o
XYYX_OPENCL_OBJS=xyyx/XYYXApp_opencl.o xyyx/XYYXWorker_opencl.o xyyx/XYYXGpuWorker_opencl.o xyyx/XYYXSparseWorker_opencl.o xyyx/XYYXSparseGpuWorker_opencl.o

CK_METAL_OBJS=carol_kynea/CarolKyneaApp_metal.o carol_kynea/CarolKyneaWorker_metal.o carol_kynea/CarolKyneaGpuWorker_metal.o
GCW_METAL_OBJS=cullen_woodall/CullenWoodallApp_metal.o cullen_woodall/CullenWoodallWorker_metal.o cullen_woodall/CullenWoodallGpuWorker_metal.o
HCW_METAL_OBJS=hyper_cullen_woodall/HyperCullenWoodallApp_gpu.o hyper_cullen_woodall/HyperCullenWoodallWorker_metal.o \
   hyper_cullen_woodall/HyperCullenWoodallSparseWorker_metal.o hyper_cullen_woodall/HyperCullenWoodallGpuWorker_metal.o hyper_cullen_woodall/HyperCullenWoodallSparseGpuWorker_metal.o
LIF_METAL_OBJS=lifchitz/LifchitzApp_metal.o lifchitz/LifchitzWorker_metal.o
MF_METAL_OBJS=multi_factorial/MultiFactorialApp_metal.o multi_factorial/MultiFactorialWorker_metal.o multi_factorial/MultiFactorialGpuWorker_metal.o
PRIM_METAL_OBJS=primorial/PrimorialApp_metal.o primorial/PrimorialWorker_metal.o primorial/PrimorialGpuWorker_metal.o
SM_METAL_OBJS=smarandache/SmarandacheApp_metal.o smarandache/SmarandacheWorker_metal.o smarandache/SmarandacheGpuWorker_metal.o
SMW_METAL_OBJS=smarandache_wellin/SmarandacheWellinApp_metal.o smarandache_wellin/SmarandacheWellinWorker_metal.o smarandache_wellin/SmarandacheWellinGpuWorker_metal.o
SR2_METAL_OBJS=sierpinski_riesel/SierpinskiRieselApp_metal.o sierpinski_riesel/AlgebraicFactorHelper_metal.o \
   sierpinski_riesel/AbstractSequenceHelper_metal.o sierpinski_riesel/AbstractWorker_metal.o \
   sierpinski_riesel/GenericSequenceHelper_metal.o sierpinski_riesel/GenericWorker_metal.o sierpinski_riesel/GenericGpuWorker_metal.o \
   sierpinski_riesel/CisOneSequenceHelper_metal.o sierpinski_riesel/CisOneWithOneSequenceHelper_metal.o \
   sierpinski_riesel/CisOneWithOneSequenceWorker_metal.o sierpinski_riesel/CisOneWithOneSequenceGpuWorker_metal.o \
   sierpinski_riesel/CisOneWithMultipleSequencesHelper_metal.o sierpinski_riesel/CisOneWithMultipleSequencesWorker_metal.o \
   sierpinski_riesel/CisOneWithMultipleSequencesGpuWorker_metal.o
XYYX_METAL_OBJS=xyyx/XYYXApp_metal.o xyyx/XYYXWorker_metal.o xyyx/XYYXGpuWorker_metal.o xyyx/XYYXSparseWorker_metal.o xyyx/XYYXSparseGpuWorker_metal.o

AIR_LIBS=alternating_factorial/af_kernel.air cullen_woodall/cw_kernel.air gfn_divisor/gfn_kernel.air \
   multi_factorial/mf_kernel.air primes_in_x/pix_kernel.air primorial/primorial_kernel.air \
   sierpinski_riesel/cisonesingle_kernel.air sierpinski_riesel/generic_kernel.air \
   smarandache/sm_kernel.air smarandache_wellin/smw_kernel.air xyyx/xyyx_kernel.air
METAL_LIBS=alternating_factorial/af_kernel.metallib cullen_woodall/cw_kernel.metallib gfn_divisor/gfn_kernel.metallib \
   multi_factorial/mf_kernel.metallib primes_in_x/pix_kernel.metallib primorial/primorial_kernel.metallib \
   sierpinski_riesel/cisonesingle_kernel.metallib sierpinski_riesel/generic_kernel.metallib \
   smarandache/sm_kernel.metallib smarandache_wellin/smw_kernel.metallib xyyx/xyyx_kernel.metallib lifchitz/lifchitz_kernel.metallib
GPU_HEADERS=alternating_factorial/af_kernel.gpu.h cullen_woodall/cw_kernel.gpu.h gfn_divisor/gfn_kernel.gpu.h \
   multi_factorial/mf_kernel.gpu.h primes_in_x/pix_kernel.gpu.h primorial/primorial_kernel.gpu.h \
   sierpinski_riesel/cisonesingle_kernel.gpu.h sierpinski_riesel/generic_kernel.gpu.h \
   smarandache/sm_kernel.gpu.h smarandache_wellin/smw_kernel.gpu.h xyyx/xyyx_kernel.gpu.h lifchitz/lifchitz_kernel.gpu.h


ALL_OBJS=$(PRIMESIEVE_OBJS) $(ASM_OBJS) $(ASM_EXT_OBJS) $(CPU_CORE_OBJS) $(OPENCL_CORE_OBJS) $(METAL_CORE_OBJS) \
   $(AF_OBJS) $(MF_OBJS) $(FBNC_OBJS) $(FKBN_OBJS) $(GFND_OBJS) $(CC_OBJS) $(CK_OBJS) \
   $(PIX_OBJS) $(XYYX_OBJS) $(KBB_OBJS) $(GCW_OBJS) $(PRIM_OBJS) $(TWIN_OBJS) \
   $(DMD_OBJS) $(SR2_OBJS) $(K1B2_OBJS) $(SG_OBJS) $(PRIM_OPENCL_OBJS) $(LIF_OBJS) \
   $(LIF_OPENCL_OBJS) $(LIF_METAL_OBJS) $(SM_OBJS) $(SMW_OBJS) \
   $(HCW_OBJS) $(HCW_OPENCL_OBJS) $(HCW_METAL_OBJS) \
   $(AF_OPENCL_OBJS) $(GCW_OPENCL_OBJS) $(GFND_OPENCL_OBJS) $(MF_OPENCL_OBJS) $(PIX_OPENCL_OBJS) \
   $(XYYX_OPENCL_OBJS) $(SR2_OPENCL_OBJS) $(PRIM_METAL_OBJS) $(SM_OPENCL_OBJS) \
   $(CW_METAL_OBJS) $(MF_METAL_OBJS) $(SR2_METAL_OBJS) $(CK_OPENCL_OBJS)

ifeq ($(strip $(HAS_METAL)),yes)
all: $(CPU_PROGS) $(OPENCL_PROGS) $(METAL_PROGS)
cpu_all: $(CPU_PROGS)
gpu_all: $(OPENCL_PROGS) $(METAL_PROGS) 
else
all: $(CPU_PROGS) $(OPENCL_PROGS)
cpu_all: $(CPU_PROGS)
gpu_all: $(OPENCL_PROGS)
endif

%.o: %.S
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) -c -o $@ $< 
   
%.o: %.cpp
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) -c -o $@ $< 

sieve/%.o: sieve/%.cpp
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(CPP_FLAGS_SIEVE) -c -o $@ $< 
   
ifeq ($(OS),Windows_NT)
xyyx/%_cpu.o: xyyx/%.cpp
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS_02) -c -o $@ $< 
    
xyyx/%_opencl.o: xyyx/%.cpp
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS_02) $(CPP_FLAGS_OPENCL) -c -o $@ $< 

primes_in_x/%_cpu.o: primes_in_x/%.cpp
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS_02) -c -o $@ $< 
    
primes_in_x/%_opencl.o: primes_in_x/%.cpp
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS_02) $(CPP_FLAGS_OPENCL) -c -o $@ $< 
   
alternating_factorial/%_cpu.o: alternating_factorial/%.cpp
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS_02) -c -o $@ $< 
    
alternating_factorial/%_opencl.o: alternating_factorial/%.cpp
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS_02) $(CPP_FLAGS_OPENCL) -c -o $@ $< 
endif

%_cpu.o: %.cpp
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) -c -o $@ $< 
    
%_opencl.o: %.cpp
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(CPP_FLAGS_OPENCL) -c -o $@ $<

%_metal.o: %.cpp
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(CPP_FLAGS_METAL) -c -o $@ $<

# This is here to allow us to verify that the .mh file will not give errors when
# initializing the device and the function in MetalKernel.InitializeFunction()
%.air: %.gpu
	$(XC) $(XC_FLAGS) -c $< -o $@ 

# This is here to allow us to verify that the .mh file will not give errors when
# initializing the device and the function in MetalKernel.InitializeFunction()
%.metallib: %.air
	$(XC) $(XC_FLAGS) $< -o $@ 

%.gpu.h: %.gpu
	$(PERL) $(PERL_FLAGS) $< > $@

# Add the targets here to build the .h files needed for the OpenCL and Metal GPU workers.
# This means that it is no longer need to explicitly use perl to create these files.  Yay!
# Only some are here right now, but it shouldn't take long to add all of them.
ifeq ($(strip $(HAS_METAL)),yes)
alternating_factorial/AlternatingFactorialGpuWorker_metal.o: alternating_factorial/af_kernel.gpu.h alternating_factorial/af_kernel.metallib
carol_kynea/CarolKynealGpuWorker_metal.o: carol_kynea/ck_kernel.gpu.h carol_kynea/ck_kernel.metallib
cullen_woodall/CullenWoodallGpuWorker_metal.o: cullen_woodall/cw_kernel.gpu.h cullen_woodall/cw_kernel.metallib
gfn_divisor/GFNDivisorGpuWorker_metal.o: gfn_divisor/gfn_kernel.gpu.h gfn_divisor/gfn_kernel.metallib
hyper_cullen_woodall/HyperCullenWoodallGpuWorker_metal.o: hyper_cullen_woodall/hcw_kernel.gpu.h hyper_cullen_woodall/hcw_kernel.metallib
lifchitz/LitchitzGpuWorker_metal.o: lifchitz/lifchitz_kernel.gpu.h lifchitz/lifchitz_kernel.metallib
multi_factorial/MultiFactorialGpuWorker_metal.o: multi_factorial/mf_kernel.gpu.h multi_factorial/mf_kernel.metallib
primes_in_x/PrimesInXGpuWorker_metal.o: primes_in_x/pix_kernel.gpu.h primes_in_x/pix_kernel.metallib
primorial/PrimorialGpuWorker_metal.o: primorial/primorial_kernel.gpu.h primorial/primorial_kernel.metallib
sierpinski_riesel/CisOneWithMultipleSequencesGpuWorker_metal.o: sierpinski_riesel/cisonemultiple_kernel.gpu.h sierpinski_riesel/cisonemultiple_kernel.metallib
sierpinski_riesel/CisOneWithOneSequenceGpuWorker_metal.o: sierpinski_riesel/cisonesingle_kernel.gpu.h sierpinski_riesel/cisonesingle_kernel.metallib
sierpinski_riesel/GenericGpuWorker_metal.o: sierpinski_riesel/generic_kernel.gpu.h sierpinski_riesel/generic_kernel.metallib
smarandache/SmarandacheGpuWorker_metal.o: smarandache/sm_kernel.gpu.h smarandache/sm_kernel.metallib
smarandache_wellin/SmarandacheWellinGpuWorker_metal.o: smarandache_wellin/smw_kernel.gpu.h smarandache_wellin/smw_kernel.metallib
xyyx/XYYXGpuWorker_metal.o: xyyx/xyyx_kernel.gpu.h xyyx/xyyx_kernel.metallib xyyx/xyyx_sparse_kernel.gpu.h xyyx/xyyx_sparse_kernel.metallib
endif

alternating_factorial/AlternatingFactorialGpuWorker_opencl.o: alternating_factorial/af_kernel.gpu.h
carol_kynea/CarolKyneaGpuWorker_opencl.o: carol_kynea/ck_kernel.gpu.h
cullen_woodall/CullenWoodallGpuWorker_opencl.o: cullen_woodall/cw_kernel.gpu.h
gfn_divisor/GFNDivisorGpuWorker_opencl.o: gfn_divisor/gfn_kernel.gpu.h
hyper_cullen_woodall/HyperCullenWoodallGpuWorker_opencl.o: hyper_cullen_woodall/hcw_kernel.gpu.h
lifchitz/LifchitzGpuWorker_opencl.o: lifchitz/lifchitz_kernel.gpu.h
multi_factorial/MultiFactorialGpuWorker_opencl.o: multi_factorial/mf_kernel.gpu.h
primes_in_x/PrimesInXGpuWorker_opencl.o: primes_in_x/pix_kernel.gpu.h
primorial/PrimorialGpuWorker_opencl.o: primorial/primorial_kernel.gpu.h
sierpinski_riesel/CisOneWithMultipleSequencesGpuWorker_opencl.o: sierpinski_riesel/cisonemultiple_kernel.gpu.h
sierpinski_riesel/CisOneWithOneSequenceGpuWorker_opencl.o: sierpinski_riesel/cisonesingle_kernel.gpu.h
sierpinski_riesel/GenericGpuWorker_opencl.o: sierpinski_riesel/generic_kernel.gpu.h
smarandache/SmarandacheGpuWorker_opencl.o: smarandache/sm_kernel.gpu.h
smarandache_wellin/SmarandacheWellinGpuWorker_opencl.o: smarandache_wellin/smw_kernel.gpu.h
xyyx/XYYXGpuWorker_opencl.o: xyyx/xyyx_kernel.gpu.h xyyx/xyyx_sparse_kernel.gpu.h

afsieve: $(CPU_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(AF_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(LD_FLAGS_STDC) -o $@ $^ $(LD_FLAGS)

afsievecl: $(OPENCL_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(AF_OPENCL_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(CPP_FLAGS_OPENCL) $(LD_FLAGS_STDC) -o $@ $^ $(LD_FLAGS_OPENCL) $(LD_FLAGS)

ccsieve: $(CPU_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(CC_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(LD_FLAGS_STDC) -o $@ $^ $(LD_FLAGS)

cksieve: $(CPU_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(CK_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(LD_FLAGS_STDC) -o $@ $^ $(LD_FLAGS)

cksievecl: $(OPENCL_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(CK_OPENCL_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(CPP_FLAGS_OPENCL) $(LD_FLAGS_STDC) -o $@ $^ $(LD_FLAGS_OPENCL) $(LD_FLAGS)

cksievemtl: $(OPENCL_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(CK_METAL_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(CPP_FLAGS_OPENMTL) $(LD_FLAGS_OPENMTL) -o $@ $^ $(LD_FLAGS_OPENMTL) $(LD_FLAGS)
   
dmdsieve: $(CPU_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(ASM_EXT_OBJS) $(DMD_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(LD_FLAGS_STDC) -o $@ $^ $(LD_FLAGS_GMP) $(LD_FLAGS)
   
fbncsieve: $(CPU_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(FBNC_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(LD_FLAGS_STDC) -o $@ $^ $(LD_FLAGS)

fkbnsieve: $(CPU_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(FKBN_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(LD_FLAGS_STDC) -o $@ $^ $(LD_FLAGS)
   
gcwsieve: $(CPU_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(GCW_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(LD_FLAGS_STDC) -o $@ $^ $(LD_FLAGS)

gcwsievecl: $(OPENCL_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(GCW_OPENCL_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(CPP_FLAGS_OPENCL) $(LD_FLAGS_STDC) -o $@ $^ $(LD_FLAGS_OPENCL) $(LD_FLAGS)

gcwsievemtl: $(OPENCL_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(GCW_METAL_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(CPP_FLAGS_OPENMTL) $(LD_FLAGS_OPENMTL) -o $@ $^ $(LD_FLAGS_OPENMTL) $(LD_FLAGS)
   
gfndsieve: $(CPU_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(ASM_EXT_OBJS) $(GFND_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(LD_FLAGS_STDC) -o $@ $^  $(LD_FLAGS_GMP) $(LD_FLAGS)

gfndsievecl: $(OPENCL_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(ASM_EXT_OBJS) $(GFND_OPENCL_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(CPP_FLAGS_OPENCL) $(LD_FLAGS_STDC) -o $@ $^  $(LD_FLAGS_GMP) $(LD_FLAGS_OPENCL) $(LD_FLAGS)

hcwsieve: $(CPU_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(HCW_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(LD_FLAGS_STDC) -o $@ $^ $(LD_FLAGS)

hcwsievecl: $(OPENCL_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(HCW_OPENCL_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(CPP_FLAGS_OPENCL) $(LD_FLAGS_STDC) -o $@ $^ $(LD_FLAGS_OPENCL) $(LD_FLAGS)

hcwsievemtl: $(METAL_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(HCW_METAL_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(CPP_FLAGS_METAL) $(LD_FLAGS_METAL) -o $@ $^ $(LD_FLAGS_METAL) $(LD_FLAGS)
   
k1b2sieve: $(CPU_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(K1B2_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(LD_FLAGS_STDC) -o $@ $^ $(LD_FLAGS)

kbbsieve: $(CPU_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(KBB_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(LD_FLAGS_STDC) -o $@ $^ $(LD_FLAGS)

lifsieve: $(CPU_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(LIF_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(LD_FLAGS_STDC) -o $@ $^ $(LD_FLAGS)

lifsievecl: $(OPENCL_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(LIF_OPENCL_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(CPP_FLAGS_OPENCL) $(LD_FLAGS_STDC) -o $@ $^ $(LD_FLAGS_OPENCL) $(LD_FLAGS)
  
lifsievemtl: $(METAL_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(LIF_METAL_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(CPP_FLAGS_METAL) $(LD_FLAGS_STDC) -o $@ $^ $(LD_FLAGS_METAL) $(LD_FLAGS)
   
mfsieve: $(CPU_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(MF_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(LD_FLAGS_STDC) -o $@ $^ $(LD_FLAGS)

mfsievecl: $(OPENCL_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(MF_OPENCL_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(CPP_FLAGS_OPENCL) $(LD_FLAGS_STDC) -o $@ $^ $(LD_FLAGS_OPENCL) $(LD_FLAGS)
  
mfsievemtl: $(METAL_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(MF_METAL_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(CPP_FLAGS_METAL) $(LD_FLAGS_STDC) -o $@ $^ $(LD_FLAGS_METAL) $(LD_FLAGS)
   
pixsieve: $(CPU_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(PIX_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(LD_FLAGS_STDC) -o $@ $^ $(LD_FLAGS)

pixsievecl: $(OPENCL_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(PIX_OPENCL_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(CPP_FLAGS_OPENCL) $(LD_FLAGS_STDC) -o $@ $^ $(LD_FLAGS_OPENCL) $(LD_FLAGS)

psieve: $(CPU_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(PRIM_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(LD_FLAGS_STDC) -o $@ $^ $(LD_FLAGS)

psievecl: $(OPENCL_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(PRIM_OPENCL_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(CPP_FLAGS_OPENCL) $(LD_FLAGS_STDC) -o $@ $^ $(LD_FLAGS_OPENCL) $(LD_FLAGS)
   
psievemtl: $(METAL_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(PRIM_METAL_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(CPP_FLAGS_METAL) $(LD_FLAGS_METAL) -o $@ $^ $(LD_FLAGS_METAL) $(LD_FLAGS)
   
sgsieve: $(CPU_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(SG_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(LD_FLAGS_STDC) -o $@ $^ $(LD_FLAGS)
   
smsieve: $(CPU_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(SM_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(LD_FLAGS_STDC) -o $@ $^ $(LD_FLAGS)
   
smsievecl: $(OPENCL_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(SM_OPENCL_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(CPP_FLAGS_OPENCL) $(LD_FLAGS_STDC) -o $@ $^ $(LD_FLAGS_OPENCL) $(LD_FLAGS)
   
smsievemtl: $(METAL_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(SM_METAL_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(CPP_FLAGS_METAL) $(LD_FLAGS_METAL) -o $@ $^ $(LD_FLAGS_METAL) $(LD_FLAGS)
   
smwsieve: $(CPU_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(SMW_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(LD_FLAGS_STDC) -o $@ $^ $(LD_FLAGS)
   
smwsievecl: $(OPENCL_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(SMW_OPENCL_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(CPP_FLAGS_OPENCL) $(LD_FLAGS_STDC) -o $@ $^ $(LD_FLAGS_OPENCL) $(LD_FLAGS)
   
smwsievemtl: $(METAL_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(SMW_METAL_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(CPP_FLAGS_METAL) $(LD_FLAGS_METAL) -o $@ $^ $(LD_FLAGS_METAL) $(LD_FLAGS)
   
srsieve2: $(CPU_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(SR2_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(LD_FLAGS_STDC) -o $@ $^ $(LD_FLAGS)
   
srsieve2cl: $(OPENCL_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(SR2_OPENCL_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(CPP_FLAGS_OPENCL) $(LD_FLAGS_STDC) -o $@ $^ $(LD_FLAGS_OPENCL) $(LD_FLAGS)

srsieve2mtl: $(METAL_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(SR2_METAL_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(CPP_FLAGS_METAL) $(LD_FLAGS_METAL) -o $@ $^ $(LD_FLAGS_METAL) $(LD_FLAGS)

twinsieve: $(CPU_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(TWIN_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(LD_FLAGS_STDC) -o $@ $^ $(LD_FLAGS)

xyyxsieve: $(CPU_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(XYYX_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(LD_FLAGS_STDC) -o $@ $^ $(LD_FLAGS)

xyyxsievecl: $(OPENCL_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(XYYX_OPENCL_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(CPP_FLAGS_OPENCL) $(LD_FLAGS_STDC) -o $@ $^ $(LD_FLAGS_OPENCL) $(LD_FLAGS)

xyyxsievemtl: $(METAL_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(XYYX_METAL_OBJS)
	$(CC) $(CPP_FLAGS) $(OPT_CPP_FLAGS) $(CPP_FLAGS_METAL) $(LD_FLAGS_METAL) -o $@ $^ $(LD_FLAGS_METAL) $(LD_FLAGS)
   
clean_objs:
	rm -f $(AIR_LIBS) $(METAL_LIBS) $(GPU_HEADERS) *.log
	rm -f $(ALL_OBJS)
   
clean:
	rm -f $(CPU_PROGS) $(OPENCL_PROGS) $(AIR_LIBS) $(METAL_LIBS) $(GPU_HEADERS) *.log
	rm -f $(ALL_OBJS)
