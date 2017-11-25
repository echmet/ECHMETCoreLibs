#ifndef ECHMET_CAES_CAES_EXTENDED_HPP
#define ECHMET_CAES_CAES_EXTENDED_HPP

#include <internal/echmetmath_internal.h>
#include <algorithm>
#include <cassert>
#include <map>

namespace ECHMET {
namespace CAES {

template <typename CAESReal>
pKaShiftedConstituent<CAESReal>::pKaShiftedConstituent(const SysComp::Constituent *constituent, const ERVec<CAESReal> &&shiftedpKas) :
	constituent(constituent),
	shiftedpKas(shiftedpKas)
{
}

template <typename CAESReal>
pBShiftedIonicForm<CAESReal>::pBShiftedIonicForm(const SysComp::IonicForm *ionicForm, const CAESReal &&shiftedpB, const CAESReal &ligandIFConcentration) :
	ionicForm(ionicForm),
	shiftedpB(shiftedpB),
	ligandIFConcentration(ligandIFConcentration)
{
}

/*!
 * Calculate dissocation degrees for complex forms. It is assumed that this is governed only by complexation equilibria
 * with the free ligand being the pivotal constituent.
 * If the function fails, the state of <tt>mapping</tt> variable is undefined.
 * <b>NOTE: This is more or less an educated guess that may yield values useless for diffusion coefficients calculations.</b>
 *
 * @param[in] shIFsVec Vector of constituents with shifted pBs.
 * @param[in] mapping Dissociation degrees mapping.
 */
template <typename CAESReal>
void calculateDDSForComplexes(const pBShiftedIonicFormsVec<CAESReal> &shIFsVec, DissocDegreesDerivativesMap &mapping, const RealVec *analyticalConcentrations)
{
	std::vector<pBShiftedIonicForm<CAESReal>> chain{};

	if (shIFsVec.empty())
		return;

	const SysComp::IonicForm *lastIF = shIFsVec.front().ionicForm;
	int32_t nucleusBaseCharge = lastIF->ancestor->totalCharge;

	assert(lastIF->ligandCount == 1);

	chain.emplace_back(shIFsVec.front());

	for (auto shIF = shIFsVec.cbegin() + 1; shIF != shIFsVec.cend(); shIF++) {
		/* Check if the current ionic form is a next element in the current chain of
		 * consecutive equilibria */
		auto isAnotherChain = [](const SysComp::IonicForm *current, const SysComp::IonicForm *last, const int32_t nucleusBaseCharge) {
			const int32_t currentNucleusBaseCharge = [](const SysComp::IonicForm *iF) {
				const SysComp::IonicForm *ancestor = iF->ancestor;

				while (ancestor->ligand != nullptr)
					ancestor = ancestor->ancestor;

				return ancestor->totalCharge;
			}(current);

			/* Brace yourself! */
			/* Check that the nucleus is the same */
			if (*current->nucleus->name != *last->nucleus->name)
				return true;

			/* Check that the nucleus charge is the same */
			if (currentNucleusBaseCharge != nucleusBaseCharge)
				return true;

			/* Check that the ligans is the same */
			if (*current->ligand->name != *last->ligand->name)
				return true;

			/* Check that the ligand charge is the same */
			if (current->ligandCharge != last->ligandCharge)
				return true;

			/* Check that the ligand count is greater by one than that of the previous form */
			if (current->ligandCount != last->ligandCount + 1)
				return true;

			return false;
		};

		if (isAnotherChain(shIF->ionicForm, lastIF, nucleusBaseCharge)) {
			assert(shIF->ionicForm->ligandCount == 1);

			processChain(chain, mapping, analyticalConcentrations);

			chain.clear();
			nucleusBaseCharge = shIF->ionicForm->ancestor->totalCharge;
		}

		lastIF = shIF->ionicForm;
		chain.emplace_back(*shIF);
	}
}

/*!
 * Calculate dissocation degrees for free forms. This is assumed to be driven by acidobazic equilibria only.
 * If the function fails, the state of <tt>mapping</tt> variable is undefined.
 *
 * @param[in] shCVec Vector of constituents with shifted pKas.
 * @param[in] cH Concentration on of \p H3O+ ion in <tt>mol/dm3</tt>.
 * @param[in] mapping Dissociation degrees mapping.
 */
template <typename CAESReal>
void calculateDDSForFreeForms(const pKaShiftedConstituentsVec<CAESReal> &shCVec, const CAESReal &cH, DissocDegreesDerivativesMap &mapping, const RealVec *analyticalConcentrations)
{
	for (const auto &shC : shCVec) {
		const auto inC = shC.constituent;
		const size_t anCIdx = inC->analyticalConcentrationIndex;
		const TotalEquilibrium<CAESReal> te{inC->chargeLow, inC->chargeHigh, shC.shiftedpKas, analyticalConcentrations->at(anCIdx) / 1000.0};

		const ERVec<CAESReal> results = calculatePeakMasterWitchcraft(te, cH);

		for (int num = te.numLow; num <= te.numHigh; num++) {
			const size_t idx = num - te.numLow;
			const SysComp::IonicForm *iF = [](const SysComp::Constituent *c, const int32_t charge) -> const SysComp::IonicForm * {
				for (size_t idx = 0; idx < c->ionicForms->size(); idx++) {
					const SysComp::IonicForm *iF = c->ionicForms->at(idx);

					if (iF->totalCharge == charge && iF->ligand == nullptr)
						return iF;
				}

				return nullptr;
			}(inC, num);

			assert(iF != nullptr);

			const std::string key(iF->name->c_str());
			assert(mapping.find(key) == mapping.end());

			mapping.emplace(key, CAESRealToECHMETReal(results.at(idx)));
		}
	}
}

/*!
 * Forget about sanity. Let's just reimplement what PeakMaster does to get to these numbers and get the hell outta here.
   I'm sorry, but _you_ get to pick up the pieces when this falls apart.
 *
 * @param[in] Total equilibrium for the given constituent, either acodibazic or complexation.
 * @param[in] Concentration of the pivotal constituent.
 *
 * @return Vector of dissociation degrees derivatives by the pivotal constituent concentration (or something like that...).
 */
template <typename CAESReal>
ERVec<CAESReal> calculatePeakMasterWitchcraft(const TotalEquilibrium<CAESReal> &te, const CAESReal &v)
{
	const CAESReal DELTA = 1.0e-25;	/* No, there is no actual reason why I picked this specific value for delta */
	const ERVec<CAESReal> dist = te.distribution(v);
	ERVec<CAESReal> results{};
	CAESReal GOneLow;
	CAESReal GOneHigh;
	CAESReal GOne;

	te.Ts(v - DELTA, GOneLow);
	te.Ts(v + DELTA, GOneHigh);
	te.Ts(v, GOne);

	const CAESReal GOneDiff = (GOneHigh - GOneLow) / (2.0 * DELTA);

	for (int num = te.numLow; num <= te.numHigh; num++) {
		const size_t idx = num - te.numLow;
		const CAESReal numerator = GOneDiff * v - GOne * num;
		const CAESReal denominator = te.Ls.at(idx) * VMath::pow<CAESReal>(v, num + 1);
		const CAESReal diffDissocDegree = -(numerator / denominator) * VMath::pow<CAESReal>(dist.at(idx), 2) / 1000.0; /* We are dividing by 1000 to scale back to mmol/dm3 domain */

		results.emplace_back(diffDissocDegree);
	}

	return results;
}

/*!
 * Wrapper function for the derivators.
 * The wrapper performs the necessary initialization and cleanup.
 *
 * @param[out] derivatives Vector to contain the results of derivations.
 * @param[in] H The H by which to perturb the system.
 * @param[in] isCorrection Correct the system composition for ionic strength.
 * @param[in] chemSystem The chemical system to solve.
 * @param[in] analyticalConcentrations Analytical concentrations of constituents.
 * @param[in] executor Function computing the derivation.
 * @param[in] params List of additional parameters for the \p executor.
 *
 * @retval RetCode::OK Success.
 * @retval RetCode::E_INVALID_ARGUMENT Invalid parameter was passed to the function.
 * @retval RetCode::E_NO_MEMORY Insufficient memory to perform the calculation.
 * @retval . Any other \p RetCode value returned by \p executor.
 */
template <typename... EParams>
RetCode derivatorSkin(RealVec *derivatives, const ECHMETReal &H, Solver *solver, const SysComp::ChemicalSystem &chemSystem, const RealVec *analyticalConcentrations, std::function<RetCode (RealVec *, const ECHMETReal &, SolverImpl<mpfr::mpreal> *, const SysComp::ChemicalSystem &, const RealVec *, const SolverVector<mpfr::mpreal> &, EParams...)> &executor, EParams... params)
{
	RetCode tRet;
	const int currentMpfrPrec = mpfr::mpreal::get_default_prec();

	if (H <= ECHMETReal{0})
		return RetCode::E_INVALID_ARGUMENT;

	if (derivatives != nullptr) {
		if (derivatives->size() != chemSystem.ionicForms->size())
			return RetCode::E_INVALID_ARGUMENT;
	}

	SolverImpl<mpfr::mpreal> *solverImpl = dynamic_cast<SolverImpl<mpfr::mpreal> *>(solver);
	if (solverImpl == nullptr)
		return RetCode::E_INVALID_ARGUMENT;

	mpfr::mpreal::set_default_prec(mpfr::digits2bits(200));

	SolverVector<mpfr::mpreal> estimatedConcentrations{};
	tRet = estimateDistributionInternal(solver->context(), analyticalConcentrations, estimatedConcentrations);
	if (tRet != RetCode::OK) {
		mpfr::mpreal::set_default_prec(currentMpfrPrec);
		return tRet;
	}

	tRet = executor(derivatives, H, solverImpl, chemSystem, analyticalConcentrations, estimatedConcentrations, params...);

	mpfr::mpreal::set_default_prec(currentMpfrPrec);

	return tRet;
}

/*!
 * Generates a vector of all constituents with pKa values shifted accordingly to the pH(pKa)-shift principle.
 *
 * @param[in] constituents Vector of all constituents.
 * @param[in] cH Concentration of H3O+ ion in mol/dm3.
 *
 * @return Vector of all constituents with recalculated pKa values
 */
template <typename CAESReal>
pKaShiftedConstituentsVec<CAESReal> makeShiftedConstituentsVector(const SysComp::ConstituentVec *constituents, const CAESReal &cH, const RealVec *ionicConcentrations)
{
	pKaShiftedConstituentsVec<CAESReal> shCVec{};

	for (size_t idx = 0; idx < constituents->size(); idx++) {
		const SysComp::Constituent *c = constituents->at(idx);

		shCVec.emplace_back(c, makeShiftedpKasVector(c, cH, ionicConcentrations));
	}

	return shCVec;
}

/*!
 * Generates a vector of all constituents with pB values shifted accordingly to the pH(pKa)-shift principle.
 *
 * @param[in] constituents Vector of all constituents.
 *
 * @return Vector of all constituents with recalculated pKa values
 */
template <typename CAESReal>
pBShiftedIonicFormsVec<CAESReal> makepBShiftedIonicFormsVector(const SysComp::ConstituentVec *constituents, const RealVec *ionicConcentrations)
{
	pBShiftedIonicFormsVec<CAESReal> shIFsVec;

	for (size_t idx = 0; idx < constituents->size(); idx++) {
		const SysComp::Constituent *c = constituents->at(idx);

		walkConstituentsIonicForms(c, shIFsVec, ionicConcentrations);
	}

	return shIFsVec;
}

/*!
 * Calculates pJ-shifted acidobazic dissociation constants.
 *
 * @param[in] c Constituent to process.
 * @param[in] cH Concentration of \p H3O+ ion in <tt>mol/dm3</tt>.
 *
 * @return Vector of pH-shifted pKa constants.
 */
template <typename CAESReal>
ERVec<CAESReal> makeShiftedpKasVector(const SysComp::Constituent *c, const CAESReal &cH, const RealVec *icVec)
{
	/* Constituent with only one charge does not undergo acidobazic equilibria */
	if (c->pKas->size() < 1)
		return ERVec<CAESReal>{};

	ERVec<CAESReal> shiftedpKas{};
	std::vector<const SysComp::IonicForm *> freeIonicFroms{};

	/* Find uncomplexed ionic forms of the constituent */
	for (size_t idx = 0; idx < c->ionicForms->size(); idx++) {
		const SysComp::IonicForm *iF = c->ionicForms->at(idx);

		if (iF->ligand == nullptr) {
			freeIonicFroms.emplace_back(iF);

			/* We assume that the ionic forms have been generated consecutively from lowest
			 * to highest charge. If this assumption breaks, so do we here. */
			ECHMET_DEBUG_CODE(if (freeIonicFroms.size() > 1)
						assert((*(freeIonicFroms.cend() - 1))->totalCharge > (*(freeIonicFroms.cend() - 2))->totalCharge););
		}
	}

	/* Calculate shifted pKas */
	for (size_t idx = 1; idx < freeIonicFroms.size(); idx++) {
		const auto &prevIF = freeIonicFroms.at(idx - 1);
		const auto &currentIF = freeIonicFroms.at(idx);

		const CAESReal Ka = (cH * icVec->at(prevIF->ionicConcentrationIndex) / 1000.0) / (icVec->at(currentIF->ionicConcentrationIndex) / 1000.0);

		shiftedpKas.emplace_back(-VMath::log10(Ka));
	}

	return shiftedpKas;
}

/*!
 * Calculates dissociation degree derivatives for a given complexation equilibira chain.
 *
 * @param[in] chain Chain containing consecutive reciprocal stability constants.
 * @param[in,out] mapping Dissociation degrees mapping.
 */
template <typename CAESReal>
void processChain(const std::vector<pBShiftedIonicForm<CAESReal>> &chain, DissocDegreesDerivativesMap &mapping, const RealVec *analyticalConcentrations)
{
	const size_t ncCIdx = chain.front().ionicForm->nucleus->analyticalConcentrationIndex;
	const CAESReal nucleusC = analyticalConcentrations->at(ncCIdx);
	const CAESReal ligandIFc = chain.front().ligandIFConcentration;
	ERVec<CAESReal> reciprocalpBs{};

	for (const auto &form : chain)
		reciprocalpBs.emplace_back(form.shiftedpB);

	TotalEquilibrium<CAESReal> te{0, static_cast<int>(chain.size()), reciprocalpBs, nucleusC};

	const ERVec<CAESReal> results = calculatePeakMasterWitchcraft(te, ligandIFc);

	for (size_t idx = 1; idx < results.size(); idx++)
		mapping.emplace(std::string(chain.at(idx - 1).ionicForm->name->c_str()), CAESRealToECHMETReal(results.at(idx)));
}

/*!
 * Walks through all ionic forms of a given constituent and calculates pH-shifted reciprocal
 * complexation stability constants.
 *
 * @param[in] c The given constituent.
 * @param[in,out] shIFVec Vector of shifted reciprocal stability constants.
 */
template <typename CAESReal>
void walkConstituentsIonicForms(const SysComp::Constituent *c, pBShiftedIonicFormsVec<CAESReal> &shIFsVec, const RealVec *icVec)
{
	for (size_t idx = 0; idx < c->ionicForms->size(); idx++) {
		const SysComp::IonicForm *iF = c->ionicForms->at(idx);

		if (iF->ancestor == nullptr)
			continue;

		const CAESReal cAncestor = icVec->at(iF->ancestor->ionicConcentrationIndex) / 1000.0;
		const CAESReal cLigand = [icVec](const SysComp::Constituent *c, const int32_t charge) {
			for (size_t idx = 0; idx < c->ionicForms->size(); idx++) {
				const SysComp::IonicForm *iF = c->ionicForms->at(idx);
				if (iF->ancestor == nullptr && iF->totalCharge == charge)
					return icVec->at(iF->ionicConcentrationIndex);
			}

			throw std::logic_error{"No free ligand ionic form found"};
		}(iF->ligand, iF->ligandCharge) / 1000.0;

		/* Use the reciprocal value of stability constant to stay consitent with TotalEquilibrium behavior */
		const CAESReal rB = 1.0 / (icVec->at(iF->ionicConcentrationIndex) / 1000.0) / (cLigand * cAncestor);

		shIFsVec.emplace_back(iF, -VMath::log10(rB), cLigand);
	}
}


} // namespace CAES
} // namespace ECHMET

#endif // ECHMET_CAES_CAES_EXTENDED_HPP
