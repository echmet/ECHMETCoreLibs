#ifndef ECHMET_CAES_ESTIMATOR_HELPERS_H
#define ECHMET_CAES_ESTIMATOR_HELPERS_H

#include "totalequilibrium.h"
#ifdef ECHMET_USE_X86_EXTENSIONS
#include "vecmath/vecmath.h"
#else
#include "genericmath.h"
#endif // ECHMET_USE_X86_EXTENSIONS

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

template <typename CAESReal, InstructionSet ISet, bool ThreadSafe>
inline
void calculateDistributionWithDerivative(const CAESReal &v,
					 CAESReal *const ECHMET_RESTRICT_PTR distribution,
					 CAESReal *const ECHMET_RESTRICT_PTR dDistdV,
					 std::vector<TotalEquilibrium<CAESReal, ISet, ThreadSafe>> &totalEquilibria,
					 const ECHMETReal *const ECHMET_RESTRICT_PTR acRaw);

#ifdef ECHMET_USE_X86_EXTENSIONS
template <>
void calculateDistributionWithDerivative<double, InstructionSet::AVX, false>
					(const double &v,
					 double *const ECHMET_RESTRICT_PTR distribution,
					 double *const ECHMET_RESTRICT_PTR dDistdV,
					 std::vector<TotalEquilibrium<double, InstructionSet::AVX, false>> &totalEquilibria,
					 const ECHMETReal *const ECHMET_RESTRICT_PTR acRaw,
					 const std::vector<double> &activityCoefficients);

template <>
void calculateDistributionWithDerivative<double, InstructionSet::AVX, false>
					(const double &v,
					 double *const ECHMET_RESTRICT_PTR distribution,
					 double *const ECHMET_RESTRICT_PTR dDistdV,
					 std::vector<TotalEquilibrium<double, InstructionSet::AVX, false>> &totalEquilibria,
					 const ECHMETReal *const ECHMET_RESTRICT_PTR acRaw);

template <>
void calculateDistributionWithDerivative<double, InstructionSet::AVX, true>
					(const double &v,
					 double *const ECHMET_RESTRICT_PTR distribution,
					 double *const ECHMET_RESTRICT_PTR dDistdV,
					 std::vector<TotalEquilibrium<double, InstructionSet::AVX, true>> &totalEquilibria,
					 const ECHMETReal *const ECHMET_RESTRICT_PTR acRaw,
					 const std::vector<double> &activityCoefficients);
template <>
void calculateDistributionWithDerivative<double, InstructionSet::AVX, true>
					(const double &v,
					 double *const ECHMET_RESTRICT_PTR distribution,
					 double *const ECHMET_RESTRICT_PTR dDistdV,
					 std::vector<TotalEquilibrium<double, InstructionSet::AVX, true>> &totalEquilibria,
					 const ECHMETReal *const ECHMET_RESTRICT_PTR acRaw);

template <>
void calculateDistributionWithDerivative<double, InstructionSet::FMA3, false>
					(const double &v,
					 double *const ECHMET_RESTRICT_PTR distribution,
					 double *const ECHMET_RESTRICT_PTR dDistdV,
					 std::vector<TotalEquilibrium<double, InstructionSet::FMA3, false>> &totalEquilibria,
					 const ECHMETReal *const ECHMET_RESTRICT_PTR acRaw,
					 const std::vector<double> &activityCoefficients);
template <>
void calculateDistributionWithDerivative<double, InstructionSet::FMA3, false>
					(const double &v,
					 double *const ECHMET_RESTRICT_PTR distribution,
					 double *const ECHMET_RESTRICT_PTR dDistdV,
					 std::vector<TotalEquilibrium<double, InstructionSet::FMA3, false>> &totalEquilibria,
					 const ECHMETReal *const ECHMET_RESTRICT_PTR acRaw);

template <>
void calculateDistributionWithDerivative<double, InstructionSet::FMA3, true>
					(const double &v,
					 double *const ECHMET_RESTRICT_PTR distribution,
					 double *const ECHMET_RESTRICT_PTR dDistdV,
					 std::vector<TotalEquilibrium<double, InstructionSet::FMA3, true>> &totalEquilibria,
					 const ECHMETReal *const ECHMET_RESTRICT_PTR acRaw,
					 const std::vector<double> &activityCoefficients);

template <>
void calculateDistributionWithDerivative<double, InstructionSet::FMA3, true>
					(const double &v,
					 double *const ECHMET_RESTRICT_PTR distribution,
					 double *const ECHMET_RESTRICT_PTR dDistdV,
					 std::vector<TotalEquilibrium<double, InstructionSet::FMA3, true>> &totalEquilibria,
					 const ECHMETReal *const ECHMET_RESTRICT_PTR acRaw);

#ifndef ECHMET_DISABLE_AVX512

template <>
void calculateDistributionWithDerivative<double, InstructionSet::AVX512, false>
					(const double &v,
					 double *const ECHMET_RESTRICT_PTR distribution,
					 double *const ECHMET_RESTRICT_PTR dDistdV,
					 std::vector<TotalEquilibrium<double, InstructionSet::AVX512, false>> &totalEquilibria,
					 const ECHMETReal *const ECHMET_RESTRICT_PTR acRaw,
					 const std::vector<double> &activityCoefficients);
template <>
void calculateDistributionWithDerivative<double, InstructionSet::AVX512, false>
					(const double &v,
					 double *const ECHMET_RESTRICT_PTR distribution,
					 double *const ECHMET_RESTRICT_PTR dDistdV,
					 std::vector<TotalEquilibrium<double, InstructionSet::AVX512, false>> &totalEquilibria,
					 const ECHMETReal *const ECHMET_RESTRICT_PTR acRaw);

template <>
void calculateDistributionWithDerivative<double, InstructionSet::AVX512, true>
					(const double &v,
					 double *const ECHMET_RESTRICT_PTR distribution,
					 double *const ECHMET_RESTRICT_PTR dDistdV,
					 std::vector<TotalEquilibrium<double, InstructionSet::AVX512, true>> &totalEquilibria,
					 const ECHMETReal *const ECHMET_RESTRICT_PTR acRaw,
					 const std::vector<double> &activityCoefficients);

template <>
void calculateDistributionWithDerivative<double, InstructionSet::AVX512, true>
					(const double &v,
					 double *const ECHMET_RESTRICT_PTR distribution,
					 double *const ECHMET_RESTRICT_PTR dDistdV,
					 std::vector<TotalEquilibrium<double, InstructionSet::AVX512, true>> &totalEquilibria,
					 const ECHMETReal *const ECHMET_RESTRICT_PTR acRaw);

#endif // ECHMET_DISABLE_AVX512

#endif // ECHMET_USE_X86_EXTENSIONS

template <typename CAESReal, typename OutputReal, InstructionSet ISet>
inline
void estimateComplexesDistribution(const CNVec<CAESReal> *const ECHMET_RESTRICT_PTR complexNuclei,
				   const LigandVec<CAESReal> *const ECHMET_RESTRICT_PTR allLigands,
				   const size_t totalLigandCopySize,
				   const CAESReal *const ECHMET_RESTRICT_PTR estConcentrations, const size_t LGBlockOffset,
				   OutputReal *estimatedConcentrations);

#ifdef ECHMET_USE_X86_EXTENSIONS
template <>
void estimateComplexesDistribution<double, double, InstructionSet::SSE2>
				(const CNVec<double> *const ECHMET_RESTRICT_PTR complexNuclei,
				 const LigandVec<double> *const ECHMET_RESTRICT_PTR allLigands,
				 const size_t totalLigandCopySize,
				 const double *const ECHMET_RESTRICT_PTR estConcentrations, const size_t LGBlockOffset,
				 double *estimatedConcentrations);

template <>
void estimateComplexesDistribution<double, double, InstructionSet::AVX>
				(const CNVec<double> *const ECHMET_RESTRICT_PTR complexNuclei,
				 const LigandVec<double> *const ECHMET_RESTRICT_PTR allLigands,
				 const size_t totalLigandCopySize,
				 const double *const ECHMET_RESTRICT_PTR estConcentrations, const size_t LGBlockOffset,
				 double *estimatedConcentrations);

template <>
void estimateComplexesDistribution<double, double, InstructionSet::FMA3>
				(const CNVec<double> *const ECHMET_RESTRICT_PTR complexNuclei,
				 const LigandVec<double> *const ECHMET_RESTRICT_PTR allLigands,
				 const size_t totalLigandCopySize,
				 const double *const ECHMET_RESTRICT_PTR estConcentrations, const size_t LGBlockOffset,
				 double *estimatedConcentrations);

#ifndef ECHMET_DISABLE_AVX512

template <>
void estimateComplexesDistribution<double, double, InstructionSet::AVX512>
				(const CNVec<double> *const ECHMET_RESTRICT_PTR complexNuclei,
				 const LigandVec<double> *const ECHMET_RESTRICT_PTR allLigands,
				 const size_t totalLigandCopySize,
				 const double *const ECHMET_RESTRICT_PTR estConcentrations, const size_t LGBlockOffset,
				 double *estimatedConcentrations);

#endif // ECHMET_DISABLE_AVX512

#endif // ECHMET_USE_X86_EXTENSIONS

template <typename CAESReal, typename OutputReal, InstructionSet ISet>
inline
void setDistributionFast(const CAESReal *const ECHMET_RESTRICT_PTR estConcentrations, const size_t count, OutputReal *estimatedConcentrations);

#ifdef ECHMET_USE_X86_EXTENSIONS
template <>
void setDistributionFast<double, double, InstructionSet::SSE2>
			(const double *const ECHMET_RESTRICT_PTR estConcentrations, const size_t count, double *estimatedConcentrations);

template <>
void setDistributionFast<double, double, InstructionSet::AVX>
			(const double *const ECHMET_RESTRICT_PTR estConcentrations, const size_t count, double *estimatedConcentrations);
template <>
void setDistributionFast<double, double, InstructionSet::FMA3>
			(const double *const ECHMET_RESTRICT_PTR estConcentrations, const size_t count, double *estimatedConcentrations);

#ifndef ECHMET_DISABLE_AVX512

template <>
void setDistributionFast<double, double, InstructionSet::AVX512>
			(const double *const ECHMET_RESTRICT_PTR estConcentrations, const size_t count, double *estimatedConcentrations);

#endif // ECHMET_DISABLE_AVX512

#endif // ECHMET_USE_X86_EXTENSIONS


} // namespace CAES
} // namespace ECHMET

#endif // ECHMET_CAES_ESTIMATOR_HELPERS_H
