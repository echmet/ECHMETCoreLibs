#ifndef ECHMET_IONPROPS_IONPROPS_WORKERS_HPP
#define ECHMET_IONPROPS_IONPROPS_WORKERS_HPP

#include <echmetphchconsts.h>
#include <internal/phchconsts_calcs.hpp>
#include <internal/echmetmath_internal.h>
#include <functional>

namespace ECHMET {
namespace IonProps {

template <typename IPReal>
class OFPackReleaser {
public:
	OFPackReleaser(typename ComputationContextImpl<IPReal>::OnsagerFuossPack *pack) :
		m_pack(pack)
	{}

	~OFPackReleaser()
	{
		m_pack->finalize();
	}

private:
	typename ComputationContextImpl<IPReal>::OnsagerFuossPack *m_pack;
};

template <typename IPReal>
class StdVectorRetriever {
public:
	static const IPReal & get(const std::vector<IPReal> &icConcs, const size_t idx)
	{
		return icConcs[idx];
	}
};

template <typename IPReal>
class RealVecRetriever {
public:
	static_assert(std::is_same<IPReal, ECHMETReal>::value, "IPReal and ECHMETReal types must be the same when using RealVecRetriever");

	static const ECHMETReal & get(const RealVec *icConcs, const size_t idx)
	{
		return icConcs->elem(idx);
	}
};

/*!
 * Corrects ionic mobilities using Onsager-Fuoss law.
 * Internal implementation.
 *
 * @param[in] ions Vector of all ions in the system.
 * @param[in] onsFuoPack Preallocated arrays needed for Onsager-Fuoss calculations.
 * @param[in] icConcs Vector if ionic concentrations.
 * @param[in] ionicStrength Ionic strength of the system in <tt>mol/dm<sup>3</sup></tt>.
 * @param[in,out] calcProps. Corresponding \p SysComp\::CalculatedProperties struct where the corrected ionic mobilities will be stored.
 */
template <typename IPReal, typename ConcentrationsVec, template <typename> class Retriever>
static
void correctIonicMobilities(typename ComputationContextImpl<IPReal>::OnsagerFuossPack *onsFuoPack, const std::vector<Ion<IPReal>> &ions, const ConcentrationsVec &icConcs, const IPReal &ionicStrength, SysComp::CalculatedProperties &calcProps)
{
	const IPReal dk = PhChConsts::dkWat * PhChConsts::eAbsVac;
	const IPReal b1 = cxpow<IPReal>(PhChConsts::e, 3) / (12.0 * M_PI) * VMath::sqrt(PhChConsts::NA / cxpow<IPReal>(dk * PhChConsts::bk * PhChConsts::Tlab, 3));
	const IPReal b2 = cxpow<IPReal>(PhChConsts::e, 2) / (6.0 * M_PI * PhChConsts::th) * VMath::sqrt(PhChConsts::NA / (dk * PhChConsts::bk * PhChConsts::Tlab));
	const std::function<IPReal (const Ion<IPReal> &)> omega = [](const Ion<IPReal> &ion) {
		return ion.limitMobility / (IPReal(PhChConsts::e) * VMath::abs(ion.charge));
	};

	const size_t N = ions.size();
	IPReal gammaTotal = 2000.0 * ionicStrength;
	typename ComputationContextImpl<IPReal>::Matrix &R = onsFuoPack->R;		/* R vectors for all constituents calculated up to sixth level */
	typename ComputationContextImpl<IPReal>::Matrix &H = onsFuoPack->H;		/* H matrix */
	std::vector<IPReal> &mu_I = onsFuoPack->mu_I;					/* mu_I = gamma_I / gammaTotal */

	/* Calculate all mu_Is */
	for (size_t idx = 0; idx < N; idx++) {
		const IPReal &conc = Retriever<IPReal>::get(icConcs, ions[idx].ionicConcentrationIndex);
		const IPReal g = conc * ions[idx].charge * ions[idx].charge;
		mu_I[idx] = g / gammaTotal;
	}
	/* Calculate H matrix */
	for (size_t j = 0; j < N; j++) {
		const IPReal w_J = omega(ions[j]);

		for (size_t i = 0; i < N; i++) {
			IPReal k = 0.0; /* KrDelta * (Sum(mu_I) * (w_I / (w_I + w_J))) */
			IPReal l; /* w_I / (w_I + w_J) */

			/* KrDelta = 1 if (j == i), 0 otherwise */
			if (j == i) {
				for (size_t idx = 0; idx < N; idx++) {
					const IPReal w_I = omega(ions[idx]);
					k += mu_I[idx] * w_I / (w_I + w_J);
				}
			}

			const IPReal w_I = omega(ions[i]);

			l = mu_I[i] * w_I / (w_I + w_J);

			H(j, i) = k + l;

			ECHMET_DEBUG_CODE(fprintf(stderr, "[%zu, %zu]: w_I = %g, w_J = %g\n", j, i, IPRealToDouble(w_I), IPRealToDouble(w_J)));
			ECHMET_DEBUG_CODE(fprintf(stderr, "H(%lu, %lu) = %g [mu_I[%zu] = %g, k = %g, l = %g]\n", j, i, IPRealToDouble(k + l), i, IPRealToDouble(mu_I[i]), IPRealToDouble(k), IPRealToDouble(l)));
		}
	}

	/* Calculate R(0) for all constituents. */
	{
		IPReal sum_zIuI = 0.0;
		IPReal sum_zIuILim = 0.0;

		for (size_t idx = 0; idx < N; idx++) {
			const Ion<IPReal> &ion = ions[idx];

			sum_zIuI += ion.charge * mu_I[idx];
			sum_zIuILim += ECHMET::VMath::abs(ion.charge / ion.limitMobility) * mu_I[idx];
		}

		const IPReal sums_ratio = sum_zIuI / sum_zIuILim;

		for (size_t idx = 0; idx < N; idx++) {
		      const Ion<IPReal> &ion = ions[idx];

		      R(idx, 0) = ion.charge - sums_ratio * VMath::abs(ion.charge / ion.limitMobility);
		}
	}

	/* Calculate higher Rs for all constituents */
	for (size_t rdx = 1; rdx < PhChConsts::RCsCount; rdx++) {
		for (size_t idx = 0; idx < N; idx++) {
			IPReal sum = 0.0;

			for (size_t sdx = 0; sdx < N; sdx++) {
				IPReal c = 2.0 * H(idx, sdx);

				if (idx == sdx)
					c -= 1.0; /* Unit matrix has "ones" in the up-left -> bottom-right diagonal */

				sum += c * R(sdx, rdx - 1);

			}

			R(idx, rdx) = sum;
		}
	}

	/* Calculate the corrected mobilities */
	{
		const IPReal sqrtGammaT = VMath::sqrt(gammaTotal);
		const IPReal denominator = (1.0 + 4.743e-2 * VMath::sqrt(gammaTotal / 2));

		for (size_t idx = 0; idx < N; idx++) {
			IPReal Rtot = 0.0;

			for (size_t n = 0; n < PhChConsts::RCsCount; n++)
				Rtot += IPReal(PhChConsts::RCs[n]) * R(idx, n);

			IPReal v = (b1 * ions[idx].limitMobility * ions[idx].charge * Rtot) + (b2 * VMath::abs(ions[idx].charge));
			v *= sqrtGammaT / denominator;

			const IPReal mobility = (ions[idx].limitMobility - v) * 1.0e9;
			if (sgn(mobility) == sgn(ions[idx].limitMobility))
				(*calcProps.ionicMobilities)[ions[idx].ionicMobilityIndex] = IPRealToECHMETReal(mobility);
		}
	}
}

template <typename IPReal>
static
void correctIonicMobilitiesViscosity(const SysComp::ChemicalSystem &chemSystem, SysComp::CalculatedProperties &calcProps, const RealVec *analyticalConcentrations)
{
	/*
	 * This implements somewhat rudimentary viscosity corrections base of the following logic:
	 * From Stokes' law:
	 *
	 * n1 / n0 = v0 / v1
	 *
	 * where n is dynamic viscosity and v is velocity.
	 * This can be directly applied to mobilities because
	 *
	 * v = u . E
	 *
	 * in electrophoresis.
	 *
	 * The ratio of n1 / n0 can be expressed as
	 *
	 * k = f * c + 1
	 *
	 * Where f is a viscosity  effect of a compound and c its concentration.
	 * Corrected mobility the becomes.
	 *
	 * u1 = u0 / f
	 *
	 * We assume that viscosity effects of all compounds in the system are additive, therefore:
	 *
	 * k = SUM(k_i * c_i) + 1
	 *
	 * Further we assume that all compounds with non-zero f are affected to the same degree and
	 * compounds with zero f are not affected at all by increased viscosity.
	 */
	auto calcIonicFormViscosityCoefficient = [](const SysComp::IonicForm *iF) {
		IPReal viscosityCoefficient = 0.0;
		const SysComp::IonicForm *next = iF;

		while (next->ancestor != nullptr) {
			viscosityCoefficient += next->ligand->viscosityCoefficient;
			next = next->ancestor;
		}

		viscosityCoefficient += iF->nucleus->viscosityCoefficient;

		return viscosityCoefficient;
	};

	const IPReal ONE{1.0};
	IPReal totalViscosityCoeff = IPReal(0.0);

	for (size_t idx = 0; idx < chemSystem.constituents->size(); idx++) {
		const SysComp::Constituent *ctuent = chemSystem.constituents->at(idx);
		const IPReal c = IPReal{analyticalConcentrations->at(ctuent->analyticalConcentrationIndex)};

		totalViscosityCoeff += IPReal(ctuent->viscosityCoefficient) * c;
	}

	const IPReal k = totalViscosityCoeff + ONE;

	for (size_t idx = 0; idx < chemSystem.ionicForms->size(); idx++) {
		const SysComp::IonicForm *iF = chemSystem.ionicForms->at(idx);

		if (iF->ifType != SysComp::IonicFormType::CONSTITUENT)
			continue;

		const IPReal viscosityCoefficient = calcIonicFormViscosityCoefficient(iF);
		if (viscosityCoefficient > 0.0) {
			const IPReal uFree = (*calcProps.ionicMobilities)[iF->ionicMobilityIndex];
			const IPReal u = uFree / k;

			(*calcProps.ionicMobilities)[iF->ionicMobilityIndex] = IPRealToECHMETReal(u);
		}
	}
}

/*!
 * Internal implementation of \p calculateConductivity() working with \p ECHMETReal type.
 *
 * @param[in] chemSystem Corresponding chemical system.
 * @param[in,out] calcProps Corresponding \p SysComp\::CalculatedProperties where the calculated conductivity will be stored.
 *
 * @return Electric conductivity of the system in <tt>S/m</tt>.
 */
static
ECHMETReal calculateConductivityWorker(const SysComp::ChemicalSystem &chemSystem, SysComp::CalculatedProperties &calcProps) noexcept
{
	const RealVec *icVec = calcProps.ionicConcentrations;
	const RealVec *imVec = calcProps.ionicMobilities;
	const SysComp::IonicFormVec *ifVec = chemSystem.ionicForms;
	ECHMETReal conductivity = 0.0;

	for (size_t idx = 0; idx < icVec->size(); idx++) {
		const SysComp::IonicForm *iF = ifVec->at(idx);

		if (iF->totalCharge == 0)
			continue;

		const ECHMETReal conductivityContrib = icVec->at(idx) * VMath::abs(iF->totalCharge) * VMath::abs(imVec->at(idx) * 1.0e-9) * PhChConsts::F;

#ifdef ECHMET_HACKS_MD
		{
			std::string name;
			if (idx == 0)
				name = "H+";
			else if (idx == 1)
				name = "OH-";
			else
				name = iF->name->c_str();
	#ifdef ECHMET_USE_HIGH_PRECISION
			fprintf(stderr, "Conductivity of %s: %s\n", name.c_str(), conductivityContrib.toString().c_str());
	#else
			fprintf(stderr, "Conductivity of %s: %.17g\n", name.c_str(), conductivityContrib);
	#endif // ECHMET_USE_HIGH_PRECISION
		}
#endif // ECHMET_HACKS_MD


		conductivity += conductivityContrib;
	}

	calcProps.conductivity = IPRealToECHMETReal(conductivity);
	return conductivity;
}

/*!
 * Templated internal implementation of \p calculateConductivity() working with \p IPReal type.
 *
 * @param[in] icConcs Ionic concentrations of all ionic forms (including uncharged ones) in the system as \p IPReal .
 * @param[in] chemSystem Corresponding chemical system.
 * @param[in,out] calcProps Corresponding \p SysComp\::CalculatedProperties where the calculated conductivity will be stored.
 *
 * @return Electric conductivity of the system in <tt>S/m</tt>.
 */
template <typename IPReal>
static
IPReal calculateConductivityWorker(const std::vector<IPReal> &icConcs, const SysComp::ChemicalSystem &chemSystem, SysComp::CalculatedProperties &calcProps) noexcept
{
	const SysComp::IonicFormVec *ifVec = chemSystem.ionicForms;
	const RealVec *imVec = calcProps.ionicMobilities;
	IPReal conductivity = 0.0;

	for (size_t idx = 0; idx < ifVec->size(); idx++) {
		const IPReal &ic = icConcs.at(idx);
		const SysComp::IonicForm *iF = ifVec->at(idx);;

		const IPReal conductivityContrib = ic * VMath::abs(iF->totalCharge) * VMath::abs(imVec->at(idx) / 1.0e9) * PhChConsts::F;

#ifdef ECHMET_HACKS_MD
		{
			std::string name;
			if (idx == 0)
				name = "H+";
			else if (idx == 1)
				name = "OH-";
			else
				name = iF->name->c_str();
			fprintf(stderr, "Conductivity of %s: %s\n", name.c_str(), conductivityContrib.toString().c_str());
		}
#endif // ECHMET_HACKS_MD

		if (iF->totalCharge == 0)
			continue;

		conductivity += conductivityContrib;
	}

	calcProps.conductivity = IPRealToECHMETReal(conductivity);
	return conductivity;
}

/*!
 * Internal \p calculateEffectiveMobilities() implementation working with \p ECHMETReal type.
 *
 * @param[in] chemSystem Corresponding chemical system.
 * @param[in] analyticalConcentrations Pointer to vector of all analytical concentrations of all compounds in the system.
 * @param[in,out] calcProps Corresponding \p SysComp\::CalculatedProperties struct where the calculated effective mobilities will be stored.
 *
 * @retval RetCode::OK Success.
 * @retval RetCode::E_NO_MEMORY Insufficient memory to perform the calculation.
 * @retval RetCode::E_BAD_INPUT Invalid chemical system.
 */
static
RetCode calculateEffectiveMobilitiesWorker(const SysComp::ChemicalSystem &chemSystem, const RealVec *analyticalConcentrations, SysComp::CalculatedProperties &calcProps) noexcept
{
	/* NOTE:
	 * There used to be a different calculation of effective mobilities.
	 * It was removed because its implementation seemed incorrect and it might
	 * have been a relic from an older implementation of ECHMETCoreLibs
	 * which used different data organization.
	 * Should the effective mobility calculation ever seem broken, see
	 * commit "fd1cd1fcdce514e2c38d17b98d62e926335f98f4" for the
	 * original implementation.
	 */

	const RealVec *icVec = calcProps.ionicConcentrations;
	const RealVec *imVec = calcProps.ionicMobilities;
	RealVec *emVec = calcProps.effectiveMobilities;

	for (size_t idx = 0; idx < chemSystem.constituents->size(); idx++) {
		const SysComp::Constituent *c = chemSystem.constituents->elem(idx);
		const size_t acIdx = c->analyticalConcentrationIndex;
		const size_t emIdx = c->effectiveMobilityIndex;
		ECHMETReal efm = 0.0;

		if (analyticalConcentrations->at(acIdx) == 0.0) {
			(*emVec)[emIdx] = 0.0;
			continue;

			/* Doing this is safer than returning inf */
		}

		for (size_t jdx = 0; jdx < c->ionicForms->size(); jdx++) {
			const SysComp::IonicForm *iF = c->ionicForms->elem(jdx);
			const size_t icIdx = iF->ionicConcentrationIndex;
			const size_t imIdx = iF->ionicMobilityIndex;
			int xfrMult;

			const auto tRet = getTransferMultiplier(c, iF, xfrMult);
			if (tRet != RetCode::OK)
				return tRet;

			ECHMET_DEBUG_CODE(fprintf(stderr, "[%s] - IF:%s: Cnt = %d, TC = %d, u =  %.4g, ic = %g, ac = %g\n", c->name->c_str(), iF->name->c_str(), xfrMult, iF->totalCharge,
															    IPRealToDouble(imVec->at(imIdx)),
															    IPRealToDouble(icVec->at(icIdx)),
															    IPRealToDouble(analyticalConcentrations->at(acIdx))));
			efm += xfrMult * sgn(iF->totalCharge) * imVec->elem(imIdx) * icVec->elem(icIdx);

		}

		(*emVec)[emIdx] = efm / analyticalConcentrations->elem(acIdx);
	}

	return RetCode::OK;
}

/*!
 * Internal \p calculatepH_direct() implementation working with \p ECHMETReal.
 *
 * @param[in] cH Concentration of H<sub>3</sub>O<sup>+</sup> ions.
 * @param[in] ionicStrength Ionic strength of the system in <tt>mol/dm<sup>3</sup></tt>.
 */
static
ECHMETReal calculatepH_directWorker(const ECHMETReal &cH, const ECHMETReal &ionicStrength) noexcept
{
	const ECHMETReal gamma = activityCoefficientInternal(ionicStrength, VMath::sqrt(ionicStrength), 1);

	return -VMath::log10(cH * gamma) + 3.0;
}

/*!
 * Templated internal \p calculatepH_direct() implementation working with \p IPReal.
 * This overload is available only if \p ECHMETReal is typedefed to an IEEE754 floating
 *
 * @param[in] cH Concentration of H<sub>3</sub>O<sup>+</sup> ions.
 * @param[in] ionicStrength Ionic strength of the system in <tt>mol/dm<sup>3</sup></tt>.
 */
template <typename IPReal>
static
IPReal calculatepH_directWorker(const IPReal &cH, const IPReal &ionicStrength) noexcept
{
	const IPReal gamma = activityCoefficientInternal(ionicStrength, VMath::sqrt(ionicStrength), 1);

	return -VMath::log10(cH * gamma) + 3.0;
}

/*!
 * Internal \p calculatepH() implementation working with \p ECHMETReal.
 *
 * @param[in] isCorrection Correct for ionic strength.
 * @param[in] calcProps Corresponding \p SysComp\::CalculatedProperties struct.
 */
static
ECHMETReal calculatepHWorker(const NonidealityCorrections corrections, const SysComp::CalculatedProperties &calcProps) noexcept
{
	const bool isCorrection = nonidealityCorrectionIsSet(corrections, NonidealityCorrectionsItems::CORR_DEBYE_HUCKEL);

	return calculatepH_directWorker(calcProps.ionicConcentrations->at(0), isCorrection ? calcProps.ionicStrength : 0.0);
}

/*!
 * Templated internal \p calculatepH() implementation working with \p IPReal.
 *
 * @param[in] isCorrection Correct for ionic strength.
 * @param[in] calcProps Corresponding \p SysComp\::CalculatedProperties struct.
 */
template <typename IPReal>
static
IPReal calculatepHWorker(const std::vector<IPReal> &icConcs, const SysComp::CalculatedProperties &calcProps, const NonidealityCorrections corrections) noexcept
{
	const bool isCorrection = nonidealityCorrectionIsSet(corrections, NonidealityCorrectionsItems::CORR_DEBYE_HUCKEL);

	return calculatepH_directWorker<IPReal>(icConcs.at(0), isCorrection ? IPReal(calcProps.ionicStrength) : 0);
}

/*!
 * Internal implentation of cumulative mobility corrections.
 *
 * @param[in] onsFuoPack Preallocated arrays needed for Onsager-Fuoss calculations
 * @param[in] chemSystem Chemical system to operate on
 * @param[in] calcProps Calculated properties of chemical system
 * @param[in] analyticalConcentrations Vector of analytical concentrations
 * @param[in] corrections Nonideality corrections to perform
 *
 * @retval \p RetCode::OK On success
 * @retval \p RetCode::E_MEMORY Insufficient memory to perform the correcions
 * @retval \p RetCode::E_DATA_TOO_LARGE Size of the data exceeds maximum size of \p std::vector
 * @retval \p RetCode::E_INVALID_ARGUMENT \p nullptr passed as \p analyticalConcentrations when viscosity correction was requested
 */
static
RetCode correctMobilitiesWorker(typename ComputationContextImpl<ECHMETReal>::OnsagerFuossPack *onsFuoPack, const std::vector<Ion<ECHMETReal>> &ions, const SysComp::ChemicalSystem &chemSystem,
				SysComp::CalculatedProperties &calcProps, const RealVec *analyticalConcentrations,
				const NonidealityCorrections corrections) noexcept
{
	OFPackReleaser<ECHMETReal> ofpRel(onsFuoPack);

	try {
		if (nonidealityCorrectionIsSet(corrections, NonidealityCorrectionsItems::CORR_VISCOSITY)) {
			if (analyticalConcentrations == nullptr)
				return RetCode::E_INVALID_ARGUMENT;
			correctIonicMobilitiesViscosity<ECHMETReal>(chemSystem, calcProps, analyticalConcentrations);
		}
		if (nonidealityCorrectionIsSet(corrections, NonidealityCorrectionsItems::CORR_ONSAGER_FUOSS))
			correctIonicMobilities<ECHMETReal, RealVec *, RealVecRetriever>(onsFuoPack, ions, calcProps.ionicConcentrations, calcProps.ionicStrength, calcProps);
	} catch (std::bad_alloc &) {
		return RetCode::E_NO_MEMORY;
	} catch (std::length_error &) {
		return RetCode::E_DATA_TOO_LARGE;
	}

	return RetCode::OK;
}

/*!
 * Templated internal implentation of cumulative mobility corrections.
 *
 * @param[in] icConcs Vector of ionic concentrations. Values must be ordered in the order used by
 *		      \p ionicConcentrations vector in \p SysComp::CalculatedProperties
 * @param[in] onsFuoPack Preallocated arrays needed for Onsager-Fuoss calculations
 * @param[in] chemSystem Chemical system to operate on
 * @param[in] calcProps Calculated properties of chemical system
 * @param[in] analyticalConcentrations Vector of analytical concentrations
 * @param[in] corrections Nonideality corrections to perform
 *
 * @retval \p RetCode::OK On success
 * @retval \p RetCode::E_MEMORY Insufficient memory to perform the correcions
 * @retval \p RetCode::E_DATA_TOO_LARGE Size of the data exceeds maximum size of \p std::vector
 * @retval \p RetCode::E_INVALID_ARGUMENT \p nullptr passed as \p analyticalConcentrations when viscosity correction was requested
 */
template <typename IPReal>
static
RetCode correctMobilitiesWorker(const std::vector<IPReal> &icConcs, typename ComputationContextImpl<IPReal>::OnsagerFuossPack *onsFuoPack, const std::vector<Ion<IPReal>> &ions,
				const SysComp::ChemicalSystem &chemSystem, SysComp::CalculatedProperties &calcProps, const RealVec *analyticalConcentrations,
			        const NonidealityCorrections corrections) noexcept
{
	OFPackReleaser<IPReal> ofpRel(onsFuoPack);

	try {
		if (nonidealityCorrectionIsSet(corrections, NonidealityCorrectionsItems::CORR_VISCOSITY)) {
			if (analyticalConcentrations == nullptr)
				return RetCode::E_INVALID_ARGUMENT;
			correctIonicMobilitiesViscosity<IPReal>(chemSystem, calcProps, analyticalConcentrations);
		}
		if (nonidealityCorrectionIsSet(corrections, NonidealityCorrectionsItems::CORR_ONSAGER_FUOSS))
			correctIonicMobilities<IPReal, std::vector<IPReal>, StdVectorRetriever>(onsFuoPack, ions, icConcs, calcProps.ionicStrength, calcProps);
	} catch (std::bad_alloc &) {
		return RetCode::E_NO_MEMORY;
	} catch (std::length_error &) {
		return RetCode::E_DATA_TOO_LARGE;
	}

	return RetCode::OK;
}

} // IonProps
} // ECHMET

#endif // ECHMET_IONPROPS_IONPROPS_WORKERS_HPP
