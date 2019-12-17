#include "cpufeatures.h"

namespace ECHMET {

std::mutex CPUFeatures::s_init_lock;
CPUFeatures * CPUFeatures::s_instance{nullptr};

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
