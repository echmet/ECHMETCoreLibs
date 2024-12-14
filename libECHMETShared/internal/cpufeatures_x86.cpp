#include "cpufeatures.h"

#include <cassert>
#include <cstring>

#ifdef ECHMET_COMPILER_MSVC
	#include <intrin.h>
	#include <immintrin.h>
#endif // ECHMET_COMPILER_MSVC

#ifdef ECHMET_DEBUG_OUTPUT
	#include <cstdio>
#endif // ECHMET_DEBUG_OUTPUT

#define ECHMET_IMPORT_INTERNAL
#include <echmetelems.h>

namespace ECHMET {

static const uint8_t SSE2_FEATURE_BIT_EDX{26};
static const uint8_t SSE3_FEATURE_BIT_ECX{0};
static const uint8_t SSSE3_FEATURE_BIT_ECX{9};
static const uint8_t SSE41_FEATURE_BIT_ECX{19};
static const uint8_t SSE42_FEATURE_BIT_ECX{20};
static const uint8_t AVX_FEATURE_BIT_ECX{28};
static const uint8_t FMA3_FEATURE_BIT_ECX{12};
static const uint8_t XSAVE_FEATURE_BIT_ECX{27};
static const uint8_t AVX2_FEATURE_BIT_EBX{5};

// Queried by cpuid EAX=07h, ECX=0
static const uint8_t AVX512F_FEATURE_BIT_EBX{16};
static const uint8_t AVX512CD_FEATURE_BIT_EBX{28};
static const uint8_t AVX512PF_FEATURE_BIT_EBX{26};
static const uint8_t AVX512ER_FEATURE_BIT_EBX{27};
static const uint8_t AVX512VL_FEATURE_BIT_EBX{31};
static const uint8_t AVX512BW_FEATURE_BIT_EBX{30};
static const uint8_t AVX512DQ_FEATURE_BIT_EBX{17};
static const uint8_t AVX512IFMA_FEATURE_BIT_EBX{21};
static const uint8_t AVX512VBM_FEATURE_BIT_ECX{1};
static const uint8_t AVX512VBM2_FEATURE_BIT_ECX{6};
static const uint8_t AVX512VNNI_FEATURE_BIT_ECX{11};
static const uint8_t AVX512BITALG_FEATURE_BIT_ECX{12};
static const uint8_t AVX512VPOPCNTDQ_FEATURE_BIT_ECX{14};
static const uint8_t AVX5124VNNIW_FEATURE_BIT_EDX{2};
static const uint8_t AVX5124FMAPS_FEATURE_BIT_EDX{3};
static const uint8_t AVX512VP2INTERSECT_FEATURE_BIT_EDX{8};
static const uint8_t AVX512FP16_FEATURE_BIT_EDX{23};
// Queried by cpuid EAX=07h, ECX=1
static const uint8_t AVX512BF16_FEATURE_BIT_EAX{05};

static const uint8_t XMM_FEATURE_BIT{1};
static const uint8_t YMM_FEATURE_BIT{2};
static const uint8_t AVX512_OPMASK_BIT{5};
static const uint8_t AVX512_HI256_BIT{6};
static const uint8_t AVX512_ZMM_HI256_BIT{7};

template <typename T>
static
bool is_bit_set(const T field, const uint8_t bit) noexcept
{
	assert((sizeof(T) * 8) > bit);

	return field & ((T)1 << bit);
}

CPUFeatures::CPUFeatures()
{
#if defined(ECHMET_COMPILER_GCC_LIKE) || defined(ECHMET_COMPILER_MINGW) || defined(ECHMET_COMPILER_MSYS) || defined(ECHMET_COMPILER_MSVC)
	m_supportedSIMD = fetch_supported_SIMD();
	m_cpu_name = fetch_cpu_name();
#else
	m_supportedSIMD = SupportedSIMD();
	m_cpu_name = "";
#endif // ECHMET_COMPILER_

}

std::string CPUFeatures::fetch_cpu_name()
{
	const size_t SZ = sizeof(uint32_t);

#if defined(ECHMET_COMPILER_GCC_LIKE) || defined(ECHMET_COMPILER_MINGW) || defined(ECHMET_COMPILER_MSYS)
	auto fetch_string_part = [](char *str, uint32_t opcode) {
#else
	auto fetch_string_part = [SZ](char *str, uint32_t opcode) {
#endif // ECHMET_COMPILER_
		union {
			int32_t block[4];
			struct {
				uint32_t eax;
				uint32_t ebx;
				uint32_t ecx;
				uint32_t edx;
			} r;
		} regs;

	#if defined(ECHMET_COMPILER_GCC_LIKE) || defined(ECHMET_COMPILER_MINGW) || defined(ECHMET_COMPILER_MSYS)
		asm("cpuid;"
		    : "=a"(regs.r.eax), "=b"(regs.r.ebx), "=c"(regs.r.ecx), "=d"(regs.r.edx)
		    : "a"(opcode)
		    : );
	#elif defined(ECHMET_COMPILER_MSVC)
		__cpuidex(regs.block, opcode, 0);
	#endif // ECHMET_COMPILER_
		memcpy(str, &regs.r.eax, SZ);
		memcpy(str + SZ, &regs.r.ebx, SZ);
		memcpy(str + (2 * SZ), &regs.r.ecx, SZ);
		memcpy(str + (3 * SZ), &regs.r.edx, SZ);
	};

	char name_array[48];
	memset(name_array, 0, sizeof(name_array));

	uint32_t cpuid_mode = 0x80000000;
	union {
		int32_t block[4];
		uint32_t ret;
	} regs;
#if defined(ECHMET_COMPILER_GCC_LIKE) || defined(ECHMET_COMPILER_MINGW) || defined(ECHMET_COMPILER_MSYS)
	asm("cpuid;"
	    : "=a"(regs.ret)
	    : "a"(cpuid_mode)
	    : "ebx", "ecx", "edx");
#elif defined(ECHMET_COMPILER_MSVC)
	__cpuidex(regs.block, cpuid_mode, 0);
#endif // ECHMET_COMPILER_
	if (regs.ret < 0x80000004)
		return "";

	char *str = name_array;
	for (uint32_t opcode = 0x80000002; opcode <= 0x80000004; opcode++) {
		fetch_string_part(str, opcode);
		str += 4 * SZ;
	}

	ECHMET_DEBUG_CODE(fprintf(stderr, "CPU name string: %s\n", name_array));

	return name_array;
}

CPUFeatures::SupportedSIMD CPUFeatures::fetch_supported_SIMD()
{
	union {
		int32_t block[4];
		struct {
			uint32_t feature_flags_eax;
			uint32_t feature_flags_ebx;
			uint32_t feature_flags_ecx;
			uint32_t feature_flags_edx;
		} r;
	} regs;
	uint32_t cpuid_mode;
	const uint32_t xgetbv_mode = 0x0;

	/* Check if CPUID supports EAX = 1 mode */
	cpuid_mode = 0x0;
	#if defined(ECHMET_COMPILER_GCC_LIKE) || defined(ECHMET_COMPILER_MINGW) || defined(ECHMET_COMPILER_MSYS)
	asm("cpuid;"
	    : "=a"(regs.r.feature_flags_eax), "=b"(regs.r.feature_flags_ebx), "=c"(regs.r.feature_flags_ecx), "=d"(regs.r.feature_flags_edx)
	    : "a"(cpuid_mode)
	    : );
#elif defined(ECHMET_COMPILER_MSVC)
	__cpuidex(regs.block, cpuid_mode, 0);
#endif // ECHMET_COMPILER_
	if (regs.r.feature_flags_eax < 1) {
		ECHMET_DEBUG_CODE(fprintf(stderr, "CPU cannot be queried for its capabilities. Assuming no extended instruction sets\n"));
		return SupportedSIMD(); /* No EAX = 1 mode for CPUID */
	}

	/* Get info from CPUID */
	cpuid_mode = 0x1; /* Check support of older instruction sets */
#if defined(ECHMET_COMPILER_GCC_LIKE) || defined(ECHMET_COMPILER_MINGW) || defined(ECHMET_COMPILER_MSYS)
	asm("cpuid;"
	    : "=a"(regs.r.feature_flags_eax), "=b"(regs.r.feature_flags_ebx), "=c"(regs.r.feature_flags_ecx), "=d"(regs.r.feature_flags_edx)
	    : "a"(cpuid_mode)
	    : );
#elif defined(ECHMET_COMPILER_MSVC)
	__cpuidex(regs.block, cpuid_mode, 0);
#endif // ECHMET_COMPILER_

	const bool cpu_has_sse2 = is_bit_set(regs.r.feature_flags_edx, SSE2_FEATURE_BIT_EDX);
	const bool cpu_has_sse3 = is_bit_set(regs.r.feature_flags_ecx, SSE3_FEATURE_BIT_ECX);
	const bool cpu_has_ssse3 = is_bit_set(regs.r.feature_flags_ecx, SSSE3_FEATURE_BIT_ECX);
	const bool cpu_has_sse41 = is_bit_set(regs.r.feature_flags_ecx, SSE41_FEATURE_BIT_ECX);
	const bool cpu_has_sse42 = is_bit_set(regs.r.feature_flags_ecx, SSE42_FEATURE_BIT_ECX);
	const bool cpu_has_avx = is_bit_set(regs.r.feature_flags_ecx, AVX_FEATURE_BIT_ECX);
	const bool cpu_has_fma = is_bit_set(regs.r.feature_flags_ecx, FMA3_FEATURE_BIT_ECX);
	const bool cpu_has_xsave = is_bit_set(regs.r.feature_flags_ecx, XSAVE_FEATURE_BIT_ECX);

	ECHMET_DEBUG_CODE(
		fprintf(stderr, "CPUID base SIMD flags:\n"
				"SSE2: %d\n"
				"SSE3: %d\n"
				"SSSE3: %d\n"
				"SSE41: %d\n"
				"SSE42: %d\n"
				"AVX: %d\n"
				"FMA3: %d\n"
				"XSAVE %d\n",
				cpu_has_sse2,
				cpu_has_sse3,
				cpu_has_ssse3,
				cpu_has_sse41,
				cpu_has_sse42,
				cpu_has_avx,
				cpu_has_fma,
				cpu_has_xsave)
	);

	cpuid_mode = 0x7; /* Check support of AVX2 and AVX512 */
	uint32_t cpuid_mode_ecx = 0x0;
#if defined(ECHMET_COMPILER_GCC_LIKE) || defined(ECHMET_COMPILER_MINGW) || defined(ECHMET_COMPILER_MSYS)
	asm("cpuid;"
	    : "=a"(regs.r.feature_flags_eax), "=b"(regs.r.feature_flags_ebx), "=c"(regs.r.feature_flags_ecx), "=d"(regs.r.feature_flags_edx)
	    : "a"(cpuid_mode), "c"(cpuid_mode_ecx)
	    : );
#elif defined(ECHMET_COMPILER_MSVC)
	__cpuidex(regs.block, cpuid_mode, cpuid_mode_ecx);
#endif // ECHMET_COMPILER_


	const bool cpu_has_avx2 = is_bit_set(regs.r.feature_flags_ebx, AVX2_FEATURE_BIT_EBX);

	const bool cpu_has_avx512f       = is_bit_set(regs.r.feature_flags_ebx, AVX512F_FEATURE_BIT_EBX);
	const bool cpu_has_avx512cd      = is_bit_set(regs.r.feature_flags_ebx, AVX512CD_FEATURE_BIT_EBX);
	const bool cpu_has_avx512pf      = is_bit_set(regs.r.feature_flags_ebx, AVX512PF_FEATURE_BIT_EBX);
	const bool cpu_has_avx512er      = is_bit_set(regs.r.feature_flags_ebx, AVX512ER_FEATURE_BIT_EBX);
	const bool cpu_has_avx512vl      = is_bit_set(regs.r.feature_flags_ebx, AVX512VL_FEATURE_BIT_EBX);
	const bool cpu_has_avx512bw      = is_bit_set(regs.r.feature_flags_ebx, AVX512BW_FEATURE_BIT_EBX);
	const bool cpu_has_avx512dq      = is_bit_set(regs.r.feature_flags_ebx, AVX512DQ_FEATURE_BIT_EBX);
	const bool cpu_has_avx512ifma    = is_bit_set(regs.r.feature_flags_ebx, AVX512IFMA_FEATURE_BIT_EBX);
	const bool cpu_has_avx512vbm     = is_bit_set(regs.r.feature_flags_ecx, AVX512VBM_FEATURE_BIT_ECX);
	const bool cpu_has_avx512vbm2    = is_bit_set(regs.r.feature_flags_ecx, AVX512VBM2_FEATURE_BIT_ECX);
	const bool cpu_has_avx512vnni    = is_bit_set(regs.r.feature_flags_ecx, AVX512VNNI_FEATURE_BIT_ECX);
	const bool cpu_has_avx512bitalg  = is_bit_set(regs.r.feature_flags_ecx, AVX512BITALG_FEATURE_BIT_ECX);
	const bool cpu_has_avx512vpopcntdq = is_bit_set(regs.r.feature_flags_ecx, AVX512VPOPCNTDQ_FEATURE_BIT_ECX);
	const bool cpu_has_avx5124vnniw = is_bit_set(regs.r.feature_flags_edx, AVX5124VNNIW_FEATURE_BIT_EDX);
	const bool cpu_has_avx5124fmaps = is_bit_set(regs.r.feature_flags_edx, AVX5124FMAPS_FEATURE_BIT_EDX);
	const bool cpu_has_avx512vp2intersect = is_bit_set(regs.r.feature_flags_edx, AVX512VP2INTERSECT_FEATURE_BIT_EDX);
	const bool cpu_has_avx512fp16 = is_bit_set(regs.r.feature_flags_edx, AVX512FP16_FEATURE_BIT_EDX);

	ECHMET_DEBUG_CODE(
		fprintf(stderr, "CPUID advanced SIMD flags:\n"
				"AVX2: %d\n",
				cpu_has_avx2)
	);

	cpuid_mode = 0x7; /* Check support for more AVX512 instructions */
	cpuid_mode_ecx = 0x1;
#if defined(ECHMET_COMPILER_GCC_LIKE) || defined(ECHMET_COMPILER_MINGW) || defined(ECHMET_COMPILER_MSYS)
	asm("cpuid;"
	    : "=a"(regs.r.feature_flags_eax), "=b"(regs.r.feature_flags_ebx), "=c"(regs.r.feature_flags_ecx), "=d"(regs.r.feature_flags_edx)
	    : "a"(cpuid_mode), "c"(cpuid_mode_ecx)
	    : );
#elif defined(ECHMET_COMPILER_MSVC)
	__cpuidex(regs.block, cpuid_mode, cpuid_mode_ecx);
#endif // ECHMET_COMPILER_

	const bool cpu_has_avx512bf16 = is_bit_set(regs.r.feature_flags_eax, AVX512BF16_FEATURE_BIT_EAX);

	/* xgetbv instruction is not available so OS level support cannot be checked.
	 * Fall back to safe defaults */
	if (!cpu_has_xsave) {
		ECHMET_DEBUG_CODE(fprintf(stderr, "No XSAVE flag, OS support cannot be queried. Assuming no extended instruction sets"));
		return SupportedSIMD();
	} else {
		/* Check what level of support was enabled by OS through XCR0 register */
	#if defined(ECHMET_COMPILER_GCC_LIKE) || defined(ECHMET_COMPILER_MINGW) || defined(ECHMET_COMPILER_MSYS)
		asm("xgetbv;"
		    : "=a"(regs.r.feature_flags_eax), "=d"(regs.r.feature_flags_edx)
		    : "c"(xgetbv_mode)
		    : );
		const bool os_xmm_aware = is_bit_set(regs.r.feature_flags_eax, XMM_FEATURE_BIT); /* 128-bit long FPU registers available */
		const bool os_avx_aware = is_bit_set(regs.r.feature_flags_eax, YMM_FEATURE_BIT) & os_xmm_aware; /* 256-bit and 128-bit long FPU registers available */
		const bool os_avx512_aware = is_bit_set(regs.r.feature_flags_eax, AVX512_OPMASK_BIT) &&
					     is_bit_set(regs.r.feature_flags_eax, AVX512_HI256_BIT) &&
					     is_bit_set(regs.r.feature_flags_eax, AVX512_ZMM_HI256_BIT);
	#elif defined(ECHMET_COMPILER_MSVC)
		const uint64_t xcr0 = _xgetbv(0);
		const bool os_xmm_aware = is_bit_set(xcr0, XMM_FEATURE_BIT);
		const bool os_avx_aware = is_bit_set(xcr0, YMM_FEATURE_BIT) & os_xmm_aware;
		const bool os_avx512_aware = is_bit_set(xcr0, AVX512_OPMASK_BIT) &&
					     is_bit_set(xcr0, AVX512_HI256_BIT) &&
					     is_bit_set(xcr0, AVX512_ZMM_HI256_BIT);
	#endif // ECHMET_COMPILER_

		ECHMET_DEBUG_CODE(
			fprintf(stderr, "OS SIMD awareness:\n"
					"XMM: %d\n"
					"YMM: %d\n"
					"512REGS: %d\n",
					os_xmm_aware,
					os_avx_aware,
					os_avx512_aware)
		);

		AVX512Sets avx512 = !os_avx512_aware
			? AVX512Sets{}
			: AVX512Sets{
				cpu_has_avx512f,
				cpu_has_avx512cd,
				cpu_has_avx512pf,
				cpu_has_avx512er,
				cpu_has_avx512vl,
				cpu_has_avx512bw,
				cpu_has_avx512dq,
				cpu_has_avx512ifma,
				cpu_has_avx512vbm,
				cpu_has_avx512vbm2,
				cpu_has_avx512vnni,
				cpu_has_avx512bitalg,
				cpu_has_avx512vpopcntdq,
				cpu_has_avx5124vnniw,
				cpu_has_avx5124fmaps,
				cpu_has_avx512vp2intersect,
				cpu_has_avx512fp16,
				cpu_has_avx512bf16
			};

		ECHMET_DEBUG_CODE(
			fprintf(stderr, "AVX512 availability:\n"
					"F   : %d\n"
					"CD  : %d\n"
					"PF  : %d\n"
					"ER  : %d\n"
					"VL  : %d\n"
					"BW  : %d\n"
					"DQ  : %d\n"
					"IFMA: %d\n"
					"VBM : %d\n",
					avx512.F, avx512.CD, avx512.PF,
					avx512.ER, avx512.VL, avx512.BW,
					avx512.DQ, avx512.IFMA, avx512.VBM);
		);

		return SupportedSIMD(cpu_has_sse2 & os_xmm_aware,
				     cpu_has_sse3 & os_xmm_aware,
				     cpu_has_ssse3 & os_xmm_aware,
				     cpu_has_sse41 & os_xmm_aware,
				     cpu_has_sse42 & os_xmm_aware,
				     cpu_has_avx & os_avx_aware,
				     cpu_has_avx2 & os_avx_aware,
				     cpu_has_fma & os_avx_aware,
				     avx512);
	}
}

} // namespace ECHMET
