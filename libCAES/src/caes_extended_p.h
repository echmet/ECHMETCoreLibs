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

typedef std::map<std::string, ECHMETReal> DissocDegreesDerivativesMap;
template <typename CAESReal>
using ERVec =  std::vector<CAESReal>;

template <typename CAESReal>
using TEVec = std::vector<TotalEquilibrium<CAESReal>>;

class DDSContextImpl : public DDSContext {
public:
	explicit DDSContextImpl(const DissocDegreesDerivativesMap &&ddsMapping);
	virtual void ECHMET_CC destroy() const noexcept override;
	virtual RetCode ECHMET_CC findDissocDegreeDerivative(ECHMETReal &value, const FixedString *ionicFormName) const noexcept override;

private:
	const DissocDegreesDerivativesMap m_ddsMapping;
};

template <typename CAESReal>
class pKaShiftedConstituent {
public:
	pKaShiftedConstituent(const SysComp::Constituent *constituent, const ERVec<CAESReal> &&shiftedpKas);

	const SysComp::Constituent *constituent;
	const ERVec<CAESReal> shiftedpKas;
};
template <typename CAESReal>
using pKaShiftedConstituentsVec = std::vector<pKaShiftedConstituent<CAESReal>>;

template <typename CAESReal>
class pBShiftedIonicForm {
public:
	pBShiftedIonicForm(const SysComp::IonicForm *ionicForm, const CAESReal &&shiftedpB, const CAESReal &ligandIFConcentration);

	const SysComp::IonicForm *ionicForm;
	const CAESReal shiftedpB;
	const CAESReal ligandIFConcentration;
};
template <typename CAESReal>
using pBShiftedIonicFormsVec = std::vector<pBShiftedIonicForm<CAESReal>>;

template <typename CAESReal>
void calculateDDSForFreeForms(const pKaShiftedConstituentsVec<CAESReal> &shCVec, const CAESReal &cH, DissocDegreesDerivativesMap &mapping, const RealVec *analyticalConcentrations);

template <typename CAESReal>
void calculateDDSForComplexes(const pBShiftedIonicFormsVec<CAESReal> &shIFsVec, DissocDegreesDerivativesMap &mapping, const RealVec *analyticalConcentrations);

template <typename CAESReal>
ERVec<CAESReal> calculatePeakMasterWitchcraft(const TotalEquilibrium<CAESReal> &te, const CAESReal &v);

template <typename... EParams>
RetCode derivatorSkin(RealVec *&derivatives, const ECHMETReal &H, const bool isCorrection, const SysComp::ChemicalSystem &chemSystem, const RealVec *analyticalConcentrations, std::function<RetCode (RealVec *, const ECHMETReal &, SolverImpl<mpfr::mpreal> *, const SysComp::ChemicalSystem &, const RealVec *, const SolverVector<mpfr::mpreal>, EParams...)> &executor, EParams... params);

RetCode calculateMixedConcentrationDerivatives(RealVec *derivatives, Solver *solver, const ECHMETReal &H, const SysComp::ChemicalSystem &chemSystem, const RealVec *analyticalConcentrations, const SysComp::Constituent *perturbedConstituentJ, const SysComp::Constituent *perturbedConstituentK);
RetCode calculateSecondConcentrationDerivatives(RealVec *derivatives, Solver *solver, const ECHMETReal &H, const SysComp::ChemicalSystem &chemSystem, const RealVec *analyticalConcentrations, const SysComp::Constituent *perturbedConstituent);

template <typename CAESReal>
ERVec<CAESReal> makeShiftedpKasVector(const SysComp::Constituent *c, const CAESReal &cH, const RealVec *icVec);

template <typename CAESReal>
pBShiftedIonicFormsVec<CAESReal> makepBShiftedIonicFormsVector(const SysComp::ConstituentVec *constituents, const RealVec *ionicConcentrations);

template <typename CAESReal>
pKaShiftedConstituentsVec<CAESReal> makeShiftedConstituentsVector(const SysComp::ConstituentVec *constituents, const CAESReal &cH, const RealVec *ionicConcentration);

template <typename CAESReal>
TEVec<CAESReal> makeShiftedTotalEquilibriaVector(const pKaShiftedConstituentsVec<CAESReal> &shCVec);

template <typename CAESReal>
void processChain(const std::vector<pBShiftedIonicForm<CAESReal>> &chain, DissocDegreesDerivativesMap &mapping);

template <typename CAESReal>
void walkConstituentsIonicForms(const SysComp::Constituent *c, pBShiftedIonicFormsVec<CAESReal> &cshIFsVec, const RealVec *icVec);


} // namespace CAES
} // namespace ECHMET

#include "caes_extended.hpp"

#endif // ECHMET_CAES_CAES_EXTENDED_P_H
