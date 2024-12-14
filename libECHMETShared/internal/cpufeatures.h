#ifndef ECHMET_CPUFEATURES_H
#define ECHMET_CPUFEATURES_H

#include <cstdint>
#include <mutex>
#include <string>

namespace ECHMET {

class CPUFeatures {
public:
	class AVX512Sets {
	public:
		AVX512Sets() noexcept;
		AVX512Sets(const bool F, const bool CD, const bool PF,
			   const bool ER, const bool VL, const bool BW,
			   const bool DQ, const bool IFMA, const bool VBM,
			   const bool VBM2, const bool VNNI, const bool BITALG,
			   const bool VPOPCNTDQ, const bool _4VNNIW, const bool _4FMAPS,
			   const bool VP2INTERSECT, const bool FP16, const bool BF16) noexcept;
		AVX512Sets(const AVX512Sets &other) noexcept;
		AVX512Sets & operator=(const AVX512Sets &other) noexcept;

		const bool F;
		const bool CD;
		const bool PF;
		const bool ER;
		const bool VL;
		const bool BW;
		const bool DQ;
		const bool IFMA;
		const bool VBM;
		const bool VBM2;
		const bool VNNI;
		const bool BITALG;
		const bool VPOPCNTDQ;
		const bool _4VNNIW;
		const bool _4FMAPS;
		const bool VP2INTERSECT;
		const bool FP16;
		const bool BF16;
	};

	class SupportedSIMD {
	public:
		SupportedSIMD() noexcept;
		SupportedSIMD(const bool SSE2, const bool SSE3, const bool SSSE3,
			      const bool SSE41, const bool SSE42,
			      const bool AVX, const bool AVX2, const bool FMA3,
			      const AVX512Sets AVX512) noexcept;
		SupportedSIMD(const SupportedSIMD &other) noexcept;
		SupportedSIMD & operator=(const SupportedSIMD &other) noexcept;

		const bool SSE2;
		const bool SSE3;
		const bool SSSE3;
		const bool SSE41;
		const bool SSE42;
		const bool AVX;
		const bool AVX2;
		const bool FMA3;
		const AVX512Sets AVX512;
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
