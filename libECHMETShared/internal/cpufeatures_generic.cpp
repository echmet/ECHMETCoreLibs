#include "cpufeatures.h"

namespace ECHMET {

CPUFeatures::CPUFeatures()
{
	m_supportedSIMD = SupportedSIMD();
	m_cpu_name = "";
}

} // namespace ECHMET
