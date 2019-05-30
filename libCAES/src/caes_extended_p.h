#ifndef ECHMET_CAES_CAES_EXTENDED_P_H
#define ECHMET_CAES_CAES_EXTENDED_P_H

#include <echmetcaes_extended.h>
#include "caes_p.h"
#include <functional>
#include <map>
#include <vector>

#define ECHMET_IMPORT_INTERNAL
#include <containers/echmetvec_p.h>
#undef ECHMET_IMPORT_INTERNAL

namespace ECHMET {
namespace CAES {

template <InstructionSet ISet, typename... EParams>
static
RetCode derivatorSkin(RealVec *&derivatives, const ECHMETReal &H, const bool isCorrection, const SysComp::ChemicalSystem &chemSystem, const RealVec *analyticalConcentrations, std::function<RetCode (RealVec *, const ECHMETReal &, SolverImpl<mpfr::mpreal, ISet> *, const SysComp::ChemicalSystem &, const RealVec *, const SolverVector<mpfr::mpreal>, EParams...)> &executor, EParams... params);

static RetCode calculateMixedConcentrationDerivatives(RealVec *derivatives, Solver *solver, const ECHMETReal &H, const SysComp::ChemicalSystem &chemSystem, const RealVec *analyticalConcentrations,
						      const SysComp::Constituent *perturbedConstituentJ, const SysComp::Constituent *perturbedConstituentK,
						      const ECHMETReal &inIonicStrength);
static RetCode calculateSecondConcentrationDerivatives(RealVec *derivatives, Solver *solver, const ECHMETReal &H, const SysComp::ChemicalSystem &chemSystem, const RealVec *analyticalConcentrations,
						       const SysComp::Constituent *perturbedConstituent,
						       const ECHMETReal &inIonicStrength);

} // namespace CAES
} // namespace ECHMET

#include "caes_extended.hpp"

#endif // ECHMET_CAES_CAES_EXTENDED_P_H
