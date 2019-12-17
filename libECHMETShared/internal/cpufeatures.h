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
