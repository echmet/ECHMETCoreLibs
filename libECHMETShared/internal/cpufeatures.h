#ifndef ECHMET_CPUFEATURES_H
#define ECHMET_CPUFEATURES_H

#include <cstdint>
#include <mutex>
#include <string>

namespace ECHMET {

class CPUFeatures {
public:
	class SupportedSIMD {
	public:
		SupportedSIMD() noexcept;
		SupportedSIMD(const bool SSE2, const bool SSE3, const bool SSSE3,
			      const bool SSE41, const bool SSE42,
			      const bool AVX, const bool AVX2, const bool AVX512,
			      const bool FMA3) noexcept;
		SupportedSIMD(const SupportedSIMD &other) noexcept;
		SupportedSIMD & operator=(const SupportedSIMD &other) noexcept;

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

	static const std::string & name();
	static const SupportedSIMD & SIMD();

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

	explicit CPUFeatures();
	std::string fetch_cpu_name();

	static void initialize();
	static SupportedSIMD fetch_supported_SIMD();

	SupportedSIMD m_supportedSIMD;
	std::string m_cpu_name;

	static std::mutex s_init_lock;
	static CPUFeatures *s_instance;
};

} // namespace ECHMET

#endif // ECHMET_CPUFEATURES_H
