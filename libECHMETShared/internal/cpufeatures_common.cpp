#include "cpufeatures.h"

namespace ECHMET {

std::mutex CPUFeatures::s_init_lock;
CPUFeatures * CPUFeatures::s_instance{nullptr};

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
