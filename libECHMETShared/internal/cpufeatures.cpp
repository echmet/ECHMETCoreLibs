#include "cpufeatures.h"

#include <cassert>
#include <cstring>

#ifdef ECHMET_COMPILER_MSVC
	#include <intrin.h>
	#include <immintrin.h>
#endif // ECHMET_COMPILER_MSVC

namespace ECHMET {

std::mutex CPUFeatures::s_init_lock;
CPUFeatures * CPUFeatures::s_instance{nullptr};

template <typename T>
static
bool is_bit_set(const T field, const uint8_t bit) noexcept
{
	assert((sizeof(T) * 8) > bit);

	return field & ((T)1 << bit);
}

CPUFeatures::SupportedSIMD::SupportedSIMD() noexcept :
	SSE2(false),
	SSE3(false),
	SSSE3(false),
	SSE41(false),
	SSE42(false),
	AVX(false),
	AVX2(false),
	AVX512(false),
	FMA3(false)
{}

CPUFeatures::SupportedSIMD::SupportedSIMD(const bool SSE2, const bool SSE3, const bool SSSE3,
					  const bool SSE41, const bool SSE42,
					  const bool AVX, const bool AVX2, const bool AVX512,
					  const bool FMA3) noexcept :
	SSE2(SSE2),
	SSE3(SSE3),
	SSSE3(SSSE3),
	SSE41(SSE41),
	SSE42(SSE42),
	AVX(AVX),
	AVX2(AVX2),
	AVX512(AVX512),
	FMA3(FMA3)
{}

CPUFeatures::SupportedSIMD::SupportedSIMD(const SupportedSIMD &other) noexcept :
	SSE2(other.SSE2),
	SSE3(other.SSE3),
	SSSE3(other.SSSE3),
	SSE41(other.SSE41),
	SSE42(other.SSE42),
	AVX(other.AVX),
	AVX2(other.AVX2),
	AVX512(other.AVX512),
	FMA3(other.FMA3)
{}

CPUFeatures::SupportedSIMD & CPUFeatures::SupportedSIMD::operator=(const SupportedSIMD &other) noexcept
{
	const_cast<bool&>(SSE2) = other.SSE2;
	const_cast<bool&>(SSE3) = other.SSE3;
	const_cast<bool&>(SSSE3) = other.SSSE3;
	const_cast<bool&>(SSE41) = other.SSE41;
	const_cast<bool&>(SSE42) = other.SSE42;
	const_cast<bool&>(AVX) = other.AVX;
	const_cast<bool&>(AVX2) = other.AVX2;
	const_cast<bool&>(AVX512) = other.AVX512;
	const_cast<bool&>(FMA3) = other.FMA3;

	return *this;
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


const std::string & CPUFeatures::name()
{
	initialize();

	return s_instance->m_cpu_name;
}

const CPUFeatures::SupportedSIMD & CPUFeatures::SIMD()
{
	initialize();

	return s_instance->m_supportedSIMD;
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

	cpuid_mode = 0x7; /* Check support of AVX2 and AVX512 */
	const uint32_t cpuid_mode_ecx = 0x0;
#if defined(ECHMET_COMPILER_GCC_LIKE) || defined(ECHMET_COMPILER_MINGW) || defined(ECHMET_COMPILER_MSYS)
	asm("cpuid;"
	    : "=a"(regs.r.feature_flags_eax), "=b"(regs.r.feature_flags_ebx), "=c"(regs.r.feature_flags_ecx), "=d"(regs.r.feature_flags_edx)
	    : "a"(cpuid_mode), "c"(cpuid_mode_ecx)
	    : );
#elif defined(ECHMET_COMPILER_MSVC)
	__cpuidex(regs.block, cpuid_mode, cpuid_mode_ecx);
#endif // ECHMET_COMPILER_


	const bool cpu_has_avx2 = is_bit_set(regs.r.feature_flags_ebx, AVX2_FEATURE_BIT_EBX);
	const bool cpu_has_avx512 = is_bit_set(regs.r.feature_flags_ebx, AVX512_FEATURE_BIT_EBX);

	/* xgetbv instruction is not available so OS level support cannot be checked.
	 * Fall back to safe defaults */
	if (!cpu_has_xsave)
		return SupportedSIMD();
	else {
		/* Check what level of support was enabled by OS through XCR0 register */
	#if defined(ECHMET_COMPILER_GCC_LIKE) || defined(ECHMET_COMPILER_MINGW) || defined(ECHMET_COMPILER_MSYS)
		asm("xgetbv;"
		    : "=a"(regs.r.feature_flags_eax), "=d"(regs.r.feature_flags_edx)
		    : "c"(xgetbv_mode)
		    : );
		const bool os_xmm_aware = is_bit_set(regs.r.feature_flags_eax, XMM_FEATURE_BIT); /* 128-bit long FPU registers available */
		const bool os_avx_aware = is_bit_set(regs.r.feature_flags_eax, YMM_FEATURE_BIT) & os_xmm_aware; /* 256-bit and 128-bit long FPU registers available */
		const bool os_avx512_aware = is_bit_set(regs.r.feature_flags_eax, AVX512_OPMASK_BIT) &
					     is_bit_set(regs.r.feature_flags_eax, AVX512_HI256_BIT) &
					     is_bit_set(regs.r.feature_flags_eax, AVX512_ZMM_HI256_BIT);
	#elif defined(ECHMET_COMPILER_MSVC)
		const uint64_t xcr0 = _xgetbv(0);
		const bool os_xmm_aware = is_bit_set(xcr0, XMM_FEATURE_BIT);
		const bool os_avx_aware = is_bit_set(xcr0, YMM_FEATURE_BIT) & os_xmm_aware;
		const bool os_avx512_aware = is_bit_set(xcr0, AVX512_OPMASK_BIT) &
					     is_bit_set(xcr0, AVX512_HI256_BIT) &
					     is_bit_set(xcr0, AVX512_ZMM_HI256_BIT);
	#endif // ECHMET_COMPILER_

		return SupportedSIMD(cpu_has_sse2 & os_xmm_aware,
				     cpu_has_sse3 & os_xmm_aware,
				     cpu_has_ssse3 & os_xmm_aware,
				     cpu_has_sse41 & os_xmm_aware,
				     cpu_has_sse42 & os_xmm_aware,
				     cpu_has_avx & os_avx_aware,
				     cpu_has_avx2 & os_avx_aware,
				     cpu_has_avx512 & os_avx512_aware,
				     cpu_has_fma & os_avx_aware);
	}
}

void CPUFeatures::initialize()
{
	s_init_lock.lock();

	if (s_instance == nullptr)
		s_instance = new CPUFeatures{};

	s_init_lock.unlock();
}

} // namespace ECHMET
