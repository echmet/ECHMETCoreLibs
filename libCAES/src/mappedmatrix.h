#ifndef ECHMET_CAES_MAPPEDMATRIX_H
#define ECHMET_CAES_MAPPEDMATRIX_H

#include "types.h"
#ifdef ECHMET_USE_X86_EXTENSIONS
#include "vecmath/vecmath.h"
#else
#include "genericmath.h"
#endif // ECHMET_USE_X86_EXTENSIONS

namespace ECHMET {
namespace CAES {

template <typename NRReal, InstructionSet ISet>
class MMTypes {
public:
	typedef Eigen::Map<SolverMatrix<NRReal>, Eigen::AlignmentType::Aligned16> Matrix;
	typedef Eigen::Map<SolverVector<NRReal>, Eigen::AlignmentType::Aligned16> Vector;
};

template <>
class MMTypes<double, InstructionSet::SSE2> {
public:
	typedef Eigen::Map<SolverMatrix<double>, Eigen::AlignmentType::Aligned16> Matrix;
	typedef Eigen::Map<SolverVector<double>, Eigen::AlignmentType::Aligned16> Vector;
};

template <>
class MMTypes<double, InstructionSet::AVX> {
public:
	typedef Eigen::Map<SolverMatrix<double>, Eigen::AlignmentType::Aligned32> Matrix;
	typedef Eigen::Map<SolverVector<double>, Eigen::AlignmentType::Aligned32> Vector;
};

template <>
class MMTypes<double, InstructionSet::FMA3> {
public:
	typedef Eigen::Map<SolverMatrix<double>, Eigen::AlignmentType::Aligned32> Matrix;
	typedef Eigen::Map<SolverVector<double>, Eigen::AlignmentType::Aligned32> Vector;
};

template <>
class MMTypes<double, InstructionSet::AVX512> {
public:
	typedef Eigen::Map<SolverMatrix<double>, Eigen::AlignmentType::Aligned64> Matrix;
	typedef Eigen::Map<SolverVector<double>, Eigen::AlignmentType::Aligned64> Vector;
};

} // namespace CAES
} // namespace ECHMET

#endif // ECHMET_CAES_MAPPEDMATRIX_H
