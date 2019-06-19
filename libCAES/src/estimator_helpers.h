#ifndef ECHMET_CAES_ESTIMATOR_HELPERS_HPP
#define ECHMET_CAES_ESTIMATOR_HELPERS_HPP

#include "totalequilibrium.h"
#include "vecmath/vecmath.h"
#include <vector>

namespace ECHMET {
namespace CAES {

template <typename CAESReal, InstructionSet ISet, bool ThreadSafe>
inline
void calculateDistributionWithDerivative(const CAESReal &v,
					 CAESReal *const ECHMET_RESTRICT_PTR distribution,
					 CAESReal *const ECHMET_RESTRICT_PTR dDistdV,
					 std::vector<TotalEquilibrium<CAESReal, ISet, ThreadSafe>> &totalEquilibria,
					 const ECHMETReal *const ECHMET_RESTRICT_PTR acRaw,
					 const std::vector<CAESReal> &activityCoefficients);

template <>
void calculateDistributionWithDerivative<double, InstructionSet::AVX, false>
					(const double &v,
					 double *const ECHMET_RESTRICT_PTR distribution,
					 double *const ECHMET_RESTRICT_PTR dDistdV,
					 std::vector<TotalEquilibrium<double, InstructionSet::AVX, false>> &totalEquilibria,
					 const ECHMETReal *const ECHMET_RESTRICT_PTR acRaw,
					 const std::vector<double> &activityCoefficients);

template <>
void calculateDistributionWithDerivative<double, InstructionSet::AVX, true>
					(const double &v,
					 double *const ECHMET_RESTRICT_PTR distribution,
					 double *const ECHMET_RESTRICT_PTR dDistdV,
					 std::vector<TotalEquilibrium<double, InstructionSet::AVX, true>> &totalEquilibria,
					 const ECHMETReal *const ECHMET_RESTRICT_PTR acRaw,
					 const std::vector<double> &activityCoefficients);

template <>
void calculateDistributionWithDerivative<double, InstructionSet::FMA3, false>
					(const double &v,
					 double *const ECHMET_RESTRICT_PTR distribution,
					 double *const ECHMET_RESTRICT_PTR dDistdV,
					 std::vector<TotalEquilibrium<double, InstructionSet::FMA3, false>> &totalEquilibria,
					 const ECHMETReal *const ECHMET_RESTRICT_PTR acRaw,
					 const std::vector<double> &activityCoefficients);

template <>
void calculateDistributionWithDerivative<double, InstructionSet::FMA3, true>
					(const double &v,
					 double *const ECHMET_RESTRICT_PTR distribution,
					 double *const ECHMET_RESTRICT_PTR dDistdV,
					 std::vector<TotalEquilibrium<double, InstructionSet::FMA3, true>> &totalEquilibria,
					 const ECHMETReal *const ECHMET_RESTRICT_PTR acRaw,
					 const std::vector<double> &activityCoefficients);

} // namespace CAES
} // namespace ECHMET

#endif // ECHMET_CAES_ESTIMATOR_HELPERS_HPP
