#ifndef ECHMET_CAES_MAPPEDMATRIX_H
#define ECHMET_CAES_MAPPEDMATRIX_H

#include "types.h"
#include "vecmath/vecmath.h"

namespace ECHMET {
namespace CAES {

template <typename NRReal, InstructionSet ISet>
class MMTypes {
public:
	typedef Eigen::Map<SolverMatrix<NRReal>, Eigen::AlignmentType::Unaligned> Matrix;
	typedef Eigen::Map<SolverVector<NRReal>, Eigen::AlignmentType::Unaligned> Vector;
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

} // namespace CAES
} // namespace ECHMET

#endif // ECHMET_CAES_MAPPEDMATRIX_H
