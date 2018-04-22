#include "cpufeatures.h"

namespace ECHMET {
namespace CAES {

std::mutex CPUFeatures::s_init_lock;
CPUFeatures * CPUFeatures::s_instance{nullptr};

}
}
