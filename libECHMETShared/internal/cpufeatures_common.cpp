#include "cpufeatures.h"

namespace ECHMET {

std::mutex CPUFeatures::s_init_lock;
CPUFeatures * CPUFeatures::s_instance{nullptr};

CPUFeatures::AVX512Sets::AVX512Sets() noexcept :
	F(false),
	CD(false),
	PF(false),
	ER(false),
	VL(false),
	BW(false),
	DQ(false),
	IFMA(false),
	VBM(false),
	VBM2(false),
	VNNI(false),
	BITALG(false),
	VPOPCNTDQ(false),
	_4VNNIW(false),
	_4FMAPS(false),
	VP2INTERSECT(false),
	FP16(false),
	BF16(false)
{}

CPUFeatures::AVX512Sets::AVX512Sets(const bool F, const bool CD, const bool PF,
				    const bool ER, const bool VL, const bool BW,
				    const bool DQ, const bool IFMA, const bool VBM,
				    const bool VBM2, const bool VNNI, const bool BITALG,
				    const bool VPOPCNTDQ, const bool _4VNNIW, const bool _4FMAPS,
				    const bool VP2INTERSECT, const bool FP16, const bool BF16) noexcept :
	F(F),
	CD(CD),
	PF(PF),
	ER(ER),
	VL(VL),
	BW(BW),
	DQ(DQ),
	IFMA(IFMA),
	VBM(VBM),
	VBM2(VBM2),
	VNNI(VNNI),
	BITALG(BITALG),
	VPOPCNTDQ(VPOPCNTDQ),
	_4VNNIW(_4VNNIW),
	_4FMAPS(_4FMAPS),
	VP2INTERSECT(VP2INTERSECT),
	FP16(FP16),
	BF16(BF16)
{}

CPUFeatures::AVX512Sets::AVX512Sets(const AVX512Sets &other) noexcept :
	F(other.F),
	CD(other.CD),
	PF(other.PF),
	ER(other.ER),
	VL(other.VL),
	BW(other.BW),
	DQ(other.DQ),
	IFMA(other.IFMA),
	VBM(other.VBM),
	VBM2(other.VBM2),
	VNNI(other.VNNI),
	BITALG(other.BITALG),
	VPOPCNTDQ(other.VPOPCNTDQ),
	_4VNNIW(other._4VNNIW),
	_4FMAPS(other._4FMAPS),
	VP2INTERSECT(other.VP2INTERSECT),
	FP16(other.FP16),
	BF16(other.BF16)
{}

CPUFeatures::AVX512Sets & CPUFeatures::AVX512Sets::operator=(const AVX512Sets &other) noexcept
{
	const_cast<bool&>(F) = other.F;
	const_cast<bool&>(CD) = other.CD;
	const_cast<bool&>(PF) = other.PF;
	const_cast<bool&>(ER) = other.ER;
	const_cast<bool&>(VL) = other.VL;
	const_cast<bool&>(BW) = other.BW;
	const_cast<bool&>(DQ) = other.DQ;
	const_cast<bool&>(IFMA) = other.IFMA;
	const_cast<bool&>(VBM) = other.VBM;
	const_cast<bool&>(VBM2) = other.VBM2;
	const_cast<bool&>(VNNI) = other.VNNI;
	const_cast<bool&>(BITALG) = other.BITALG;
	const_cast<bool&>(VPOPCNTDQ) = other.VPOPCNTDQ;
	const_cast<bool&>(_4VNNIW) = other._4VNNIW;
	const_cast<bool&>(_4FMAPS) = other._4FMAPS;
	const_cast<bool&>(VP2INTERSECT) = other.VP2INTERSECT;
	const_cast<bool&>(FP16) = other.FP16;
	const_cast<bool&>(BF16) = other.BF16;

	return *this;
}

CPUFeatures::SupportedSIMD::SupportedSIMD() noexcept :
	SSE2(false),
	SSE3(false),
	SSSE3(false),
	SSE41(false),
	SSE42(false),
	AVX(false),
	AVX2(false),
	FMA3(false),
	AVX512(AVX512Sets{})
{}

CPUFeatures::SupportedSIMD::SupportedSIMD(const bool SSE2, const bool SSE3, const bool SSSE3,
					  const bool SSE41, const bool SSE42,
					  const bool AVX, const bool AVX2, const bool FMA3,
					  const AVX512Sets AVX512) noexcept :
	SSE2(SSE2),
	SSE3(SSE3),
	SSSE3(SSSE3),
	SSE41(SSE41),
	SSE42(SSE42),
	AVX(AVX),
	AVX2(AVX2),
	FMA3(FMA3),
	AVX512(AVX512)
{}

CPUFeatures::SupportedSIMD::SupportedSIMD(const SupportedSIMD &other) noexcept :
	SSE2(other.SSE2),
	SSE3(other.SSE3),
	SSSE3(other.SSSE3),
	SSE41(other.SSE41),
	SSE42(other.SSE42),
	AVX(other.AVX),
	AVX2(other.AVX2),
	FMA3(other.FMA3),
	AVX512(other.AVX512)
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
	const_cast<bool&>(FMA3) = other.FMA3;
	const_cast<AVX512Sets&>(AVX512) = other.AVX512;

	return *this;
}

void CPUFeatures::initialize()
{
	s_init_lock.lock();

	if (s_instance == nullptr)
		s_instance = new CPUFeatures{};

	s_init_lock.unlock();
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

} // namespace ECHMET
