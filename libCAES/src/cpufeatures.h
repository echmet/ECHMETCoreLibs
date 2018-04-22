#ifndef ECHMET_CAES_CPUFEATURES_H
#define ECHMET_CAES_CPUFEATURES_H

#include <cstdint>
#include <cassert>
#include <cstring>
#include <mutex>
#include <string>

namespace ECHMET {
namespace CAES {

class CPUFeatures {
public:
	class SupportedSIMD {
	public:
		SupportedSIMD() :
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

		SupportedSIMD(const bool SSE2, const bool SSE3, const bool SSSE3,
			      const bool SSE41, const bool SSE42,
			      const bool AVX, const bool AVX2, const bool AVX512,
			      const bool FMA3) :
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

		SupportedSIMD & operator=(const SupportedSIMD &other) noexcept
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

		const bool SSE2;
		const bool SSE3;
		const bool SSSE3;
		const bool SSE41;
		const bool SSE42;
		const bool AVX;
		const bool AVX2;
		const bool AVX512;
		const bool FMA3;
	};

	static const std::string & name()
	{
		initialize();

		return s_instance->m_cpu_name;
	}

	static const SupportedSIMD & SIMD()
	{
		initialize();

		return s_instance->m_supportedSIMD;
	}


private:
	static const uint8_t SSE2_FEATURE_BIT_EDX = 26;
	static const uint8_t SSE3_FEATURE_BIT_ECX = 0;
	static const uint8_t SSSE3_FEATURE_BIT_ECX = 9;
	static const uint8_t SSE41_FEATURE_BIT_ECX = 19;
	static const uint8_t SSE42_FEATURE_BIT_ECX = 20;
	static const uint8_t AVX_FEATURE_BIT_ECX = 28;
	static const uint8_t FMA3_FEATURE_BIT_ECX = 12;
	static const uint8_t XSAVE_FEATURE_BIT_ECX = 26;
	static const uint8_t AVX2_FEATURE_BIT_EBX = 5;
	static const uint8_t AVX512_FEATURE_BIT_EBX = 16;

	static const uint8_t XMM_FEATURE_BIT = 1;
	static const uint8_t YMM_FEATURE_BIT = 2;
	static const uint8_t AVX512_OPMASK_BIT = 5;
	static const uint8_t AVX512_HI256_BIT = 6;
	static const uint8_t AVX512_ZMM_HI256_BIT = 7;

	explicit CPUFeatures()
	{
		uint32_t feature_flags_eax;
		uint32_t feature_flags_ebx;
		uint32_t feature_flags_ecx;
		uint32_t feature_flags_edx;
		uint32_t cpuid_mode;
		const uint32_t xgetbv_mode = 0x0;

		/* Get info from CPUID */
		cpuid_mode = 0x1; /* Check support of older instruction sets */
		asm("cpuid;"
		    : "=a"(feature_flags_eax), "=b"(feature_flags_ebx), "=c"(feature_flags_ecx), "=d"(feature_flags_edx)
		    : "a"(cpuid_mode)
		    : );

		const bool cpu_has_sse2 = is_bit_set(feature_flags_edx, SSE2_FEATURE_BIT_EDX);
		const bool cpu_has_sse3 = is_bit_set(feature_flags_ecx, SSE3_FEATURE_BIT_ECX);
		const bool cpu_has_ssse3 = is_bit_set(feature_flags_ecx, SSSE3_FEATURE_BIT_ECX);
		const bool cpu_has_sse41 = is_bit_set(feature_flags_ecx, SSE41_FEATURE_BIT_ECX);
		const bool cpu_has_sse42 = is_bit_set(feature_flags_ecx, SSE42_FEATURE_BIT_ECX);
		const bool cpu_has_avx = is_bit_set(feature_flags_ecx, AVX_FEATURE_BIT_ECX);
		const bool cpu_has_fma = is_bit_set(feature_flags_ecx, FMA3_FEATURE_BIT_ECX);
		const bool cpu_has_xsave = is_bit_set(feature_flags_ecx, XSAVE_FEATURE_BIT_ECX);

		cpuid_mode = 0x7; /* Check support of AVX2 and AVX512 */
		const uint32_t cpuid_mode_ecx = 0x0;
		asm("cpuid;"
		    : "=a"(feature_flags_eax), "=b"(feature_flags_ebx), "=c"(feature_flags_ecx), "=d"(feature_flags_edx)
		    : "a"(cpuid_mode), "c"(cpuid_mode_ecx)
		    : );

		const bool cpu_has_avx2 = is_bit_set(feature_flags_ebx, AVX2_FEATURE_BIT_EBX);
		const bool cpu_has_avx512 = is_bit_set(feature_flags_ebx, AVX512_FEATURE_BIT_EBX);

		m_cpu_name = fetch_cpu_name();

		/* xgetbv instruction is not available so OS level support cannot be checked.
		 * Fall back to safe defaults */
		if (!cpu_has_xsave)
			m_supportedSIMD = SupportedSIMD();
		else {
			/* Check what level of support was enabled by OS through XCR0 register */
			asm("xgetbv;"
			    : "=a"(feature_flags_eax), "=d"(feature_flags_edx)
			    : "c"(xgetbv_mode)
			    : );
			const bool os_xmm_aware = is_bit_set(feature_flags_eax, XMM_FEATURE_BIT); /* 128-bit long FPU registers available */
			const bool os_avx_aware = is_bit_set(feature_flags_eax, YMM_FEATURE_BIT) & os_xmm_aware; /* 256-bit and 128-bit long FPU registers available */
			const bool os_avx512_aware = is_bit_set(feature_flags_eax, AVX512_OPMASK_BIT) &
						     is_bit_set(feature_flags_eax, AVX512_HI256_BIT) &
						     is_bit_set(feature_flags_eax, AVX512_ZMM_HI256_BIT);

			m_supportedSIMD = SupportedSIMD(
						cpu_has_sse2 & os_xmm_aware,
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

	std::string fetch_cpu_name()
	{
		const size_t SZ = sizeof(uint32_t);

		auto fetch_string_part = [](char *str, uint32_t opcode) {
			uint32_t eax;
			uint32_t ebx;
			uint32_t ecx;
			uint32_t edx;

			asm("cpuid;"
			    : "=a"(eax), "=b"(ebx), "=c"(ecx), "=d"(edx)
			    : "a"(opcode)
			    : );
			memcpy(str, &eax, SZ);
			memcpy(str + SZ, &ebx, SZ);
			memcpy(str + (2 * SZ), &ecx, SZ);
			memcpy(str + (3 * SZ), &edx, SZ);
		};

		char name_array[48];
		name_array[47] = 0;

		uint32_t cpuid_mode = 0x80000000;
		uint32_t ret;
		asm("cpuid;"
		    : "=a"(ret)
		    : "a"(cpuid_mode)
		    : "ebx", "ecx", "edx");

		if (ret < 0x80000004)
			return "";

		char *str = name_array;
		for (uint32_t opcode = 0x80000002; opcode <= 0x80000004; opcode++) {
			fetch_string_part(str, opcode);
			str += 4 * SZ;
		}

		return name_array;
	}

	static void initialize()
	{
		s_init_lock.lock();

		if (s_instance == nullptr)
			s_instance = new CPUFeatures{};

		s_init_lock.unlock();
	}

	template <typename T>
	bool is_bit_set(const T field, const uint8_t bit) const noexcept
	{
		assert((sizeof(T) * 8) > bit);

		return field & (1 << bit);
	}

	SupportedSIMD m_supportedSIMD;
	std::string m_cpu_name;

	static std::mutex s_init_lock;
	static CPUFeatures *s_instance;
};

} // namespace CAES
} // namespace ECHMET

#endif // ECHMET_CAES_CPUFEATURES_H
