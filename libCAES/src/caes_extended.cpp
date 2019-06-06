#include "caes_extended_p.h"
#define ECHMET_IMPORT_INTERNAL
#include <echmetionprops_extended.h>
#undef ECHMET_IMPORT_INTERNAL

namespace ECHMET {
namespace CAES {

class PrecisionResetter {
public:
	explicit PrecisionResetter(const long prec) :
		m_prec{prec}
	{}
	~PrecisionResetter()
	{
		mpfr::mpreal::set_default_prec(m_prec);
	}

private:
	const long m_prec;
};

template<typename TN, typename TH>
static
TN firstDerivativeCalculator2O(const TN &low, const TN &high, const TH &H)
{
#ifdef ECHMET_PLATFORM_WIN32
	const TH TWO{2.0};
#else
	thread_local static const TH TWO{2.0};
#endif // ECHMET_PLATFORM_WIN32

	const TN t = (high - low) / (TWO * H);

	return t;
}

template<typename TN, typename TH>
static
TN mixedDerivativeCalculator2O(const TN &LL, const TN &LH, const TN &HL, const TN &HH, const TH &h1, const TH &h2)
{
	const TN TOne = HH + LL;
	const TN TTwo = HL + LH;

	const TN r = (TOne - TTwo) / (4.0 * h1 * h2);

	return r;
}

template<typename TN, typename TH>
static
TN secondDerivativeCalculator2O(const TN &low, const TN &center, const TN &high, const TH &H)
{
	const TN t = (low - 2.0 * center + high) / (H * H);

	return t;
}

static
RetCode calculatepHResponse(ECHMETReal &bufferCapacity, const ECHMETReal &H, const NonidealityCorrections corrections, const SysComp::ChemicalSystem &chemSystem, const SysComp::CalculatedProperties &calcProps,
			    const RealVec *analyticalConcentrations, const SysComp::Constituent *perturbedConstituent)
{
	typedef std::tuple<int32_t, ECHMETReal> ContributingForm;
	typedef std::function<RetCode (RealVec *, const ECHMETReal &, SolverImpl<mpfr::mpreal, InstructionSet::GENERIC, true> *, const SysComp::ChemicalSystem &, const RealVec *, const MPRealVecWrap &, const SysComp::Constituent *)> ExFunc;

	const size_t N = chemSystem.ionicForms->size();
	const ECHMETReal &cAc = analyticalConcentrations->at(perturbedConstituent->analyticalConcentrationIndex);
	const ECHMETReal inIonicStrength = calcProps.ionicStrength;
	mpfr::mpreal pHHigh;
	mpfr::mpreal pHLow;
	SolverContext *solverCtx;
	Solver *solver;
	RetCode tRet;
	const std::vector<ContributingForm> contributingForms = [&cAc, &calcProps](const SysComp::IonicFormVec *iFVec) {
		std::vector<ContributingForm> forms{};

		for (size_t idx = 0; idx < iFVec->size(); idx++) {
			std::function<const SysComp::IonicForm * (const SysComp::IonicForm *)> getFreeForm = [&getFreeForm](const SysComp::IonicForm *iF) {
				if (iF->ancestor == nullptr)
					return iF;
				return getFreeForm(iF->ancestor);
			};

			const SysComp::IonicForm *iF = iFVec->at(idx);
			const SysComp::IonicForm *baseIF;

			if (iF->ancestor != nullptr)
				baseIF = getFreeForm(iF);
			else
				baseIF = iF;

			const ECHMETReal &ifC = calcProps.ionicConcentrations->at(iF->ionicConcentrationIndex);
			forms.emplace_back(baseIF->totalCharge, ifC / cAc);
		}

		return forms;
	}(perturbedConstituent->ionicForms);

	const mpfr::mpreal cHDelta = [&contributingForms, H]() {
		mpfr::mpreal total = 0.0;
		for (const auto & cf : contributingForms) {
			const int32_t z = std::get<0>(cf);
			const mpfr::mpreal x = std::get<1>(cf);

			total += z * x * 2.0 * H;
		}

		return total;
	}();

	tRet = createSolverContextInternal<mpfr::mpreal>(solverCtx, chemSystem);
	if (tRet != RetCode::OK)
		return tRet;

	try {
		const Solver::Options opts = Solver::defaultOptions();
		solver = new SolverImpl<mpfr::mpreal, InstructionSet::GENERIC, true>(static_cast<SolverContextImpl<mpfr::mpreal> *>(solverCtx), opts, corrections);
	} catch (std::bad_alloc &) {
		solverCtx->destroy();

		return RetCode::E_NO_MEMORY;
	}

	ExFunc executor = [&pHHigh, &pHLow, N, corrections, &inIonicStrength](RealVec *, const ECHMETReal &H, SolverImpl<mpfr::mpreal, InstructionSet::GENERIC, true> *solver, const SysComp::ChemicalSystem &chemSystem, const RealVec *analyticalConcentrations, const MPRealVecWrap &estimatedConcentrations,
			     const SysComp::Constituent *perturbedConstituent) {
		const size_t anCIdx = perturbedConstituent->analyticalConcentrationIndex;
		const mpfr::mpreal HH{H};
		RetCode tRet;

		auto solvePerturbed = [corrections, &inIonicStrength](SolverVector<mpfr::mpreal> &perturbedConcentrations, mpfr::mpreal &pH, SolverImpl<mpfr::mpreal, InstructionSet::GENERIC, true> *solver, const SolverVector<mpfr::mpreal> *inConcentrations, const SysComp::ChemicalSystem &chemSystem, const MPRealVecWrap &estimatedConcentrations) {
			mpfr::mpreal ionicStrength = mpfr::mpreal(inIonicStrength);
			SysComp::CalculatedProperties calcProps;

			std::vector<mpfr::mpreal> ics{};
			try {
				ics.resize(perturbedConcentrations.rows());
			} catch (std::bad_alloc &) {
				return RetCode::E_NO_MEMORY;
			}

			RetCode tRet = SysComp::initializeCalculatedProperties(calcProps, chemSystem);
			if (tRet != RetCode::OK)
				return tRet;

			tRet = solver->solveRaw(perturbedConcentrations, ionicStrength, inConcentrations, estimatedConcentrations.get(), 2000);
			if (tRet != RetCode::OK) {
				releaseCalculatedProperties(calcProps);

				return tRet;
			}

			for (int row = 0; row < perturbedConcentrations.rows(); row++)
				ics[row] = perturbedConcentrations(row);
			calcProps.ionicStrength = CAESRealToECHMETReal(ionicStrength);

			IonProps::ComputationContext *ionCtx = IonProps::makeComputationContextExtended<mpfr::mpreal>(chemSystem, ics, IonProps::ComputationContext::NONE);
			if (ionCtx == nullptr) {
				releaseCalculatedProperties(calcProps);

				return RetCode::E_NO_MEMORY;
			}

			pH = IonProps::calculatepHInternal<mpfr::mpreal>(ionCtx, corrections, calcProps);

			releaseCalculatedProperties(calcProps);
			ionCtx->destroy();

			return tRet;
		};

		SolverVector<mpfr::mpreal> icConcs{N};
		SolverVector<mpfr::mpreal> inConcentrations{};
		inConcentrations.resize(analyticalConcentrations->size());
		for (size_t idx = 0; idx < analyticalConcentrations->size(); idx++)
			inConcentrations(idx) = analyticalConcentrations->at(idx);



		inConcentrations(anCIdx) -= H;
		tRet = solvePerturbed(icConcs, pHLow, solver, &inConcentrations, chemSystem, estimatedConcentrations);
		if (tRet != RetCode::OK)
			return tRet;

		inConcentrations(anCIdx) += 2.0 * H;
		return solvePerturbed(icConcs, pHHigh, solver, &inConcentrations, chemSystem, estimatedConcentrations);
	};

	tRet = derivatorSkin<InstructionSet::GENERIC, const SysComp::Constituent *>(nullptr, H, solver, chemSystem, analyticalConcentrations, inIonicStrength, executor, perturbedConstituent);
	solverCtx->destroy();
	solver->destroy();
	if (tRet != RetCode::OK)
		return tRet;

	const mpfr::mpreal pHDelta = mpfr::abs(pHHigh - pHLow);

	/* NOTE:
	 * We are calculating the pH delta from pH values
	 * of the perturbed solutions.
	 * If the IS correction is enabled, the pH is calculated
	 * as -log10(a(H+)) and the high and low values are then subtracted.
	 * This is not the same thing as doing
	 * -log10([H+]_high - [H+]_low)) - keep this in mind!
	 */
	bufferCapacity = std::abs((cHDelta / pHDelta).toDouble());

	mpfr_free_cache();

	return RetCode::OK;
}

static
RetCode prepareDerivatorContextInternal(RealVec *&derivatives, Solver *&solver, const SysComp::ChemicalSystem &chemSystem, const NonidealityCorrections corrections,
					const bool freeCache) noexcept
{
	const size_t N = chemSystem.ionicForms->size();

	SolverContext *solverCtx;
	RetCode tRet = createSolverContextInternal<mpfr::mpreal>(solverCtx, chemSystem);
	if (tRet != RetCode::OK)
		return tRet;

	try {
		const Solver::Options opts = Solver::defaultOptions();
		solver = new SolverImpl<mpfr::mpreal, InstructionSet::GENERIC, true>(static_cast<SolverContextImpl<mpfr::mpreal> *>(solverCtx), opts, corrections);
	} catch (std::bad_alloc &) {
		solverCtx->destroy();

		return RetCode::E_NO_MEMORY;
	}

	derivatives = createRealVec(N);
	if (derivatives == nullptr) {
		solver->destroy();
		solverCtx->destroy();

		return RetCode::E_NO_MEMORY;
	}

	tRet = static_cast<VecImpl<ECHMETReal, false> *>(derivatives)->resize(N);
	if (tRet != RetCode::OK) {
		derivatives->destroy();
		solver->destroy();
		solverCtx->destroy();

		return tRet;
	}

	if (freeCache)
		mpfr_free_cache();

	return RetCode::OK;
}

RetCode ECHMET_CC calculateBufferCapacity(ECHMETReal &bufferCapacity, const NonidealityCorrections corrections, const SysComp::ChemicalSystem &chemSystem, const SysComp::CalculatedProperties &calcProps,
					  const RealVec *analyticalConcentrations) noexcept
{
	const double H{1.0e-30};
	const SysComp::Constituent *perturbedConstituent = nullptr;

	for (size_t idx = 0; idx < chemSystem.constituents->size(); idx++) {
		const SysComp::Constituent *ctuent = chemSystem.constituents->at(idx);
		const ECHMETReal &cAc = analyticalConcentrations->at(ctuent->analyticalConcentrationIndex);

		if ((ctuent->chargeLow < 0 || ctuent->chargeHigh > 0) && (cAc > 10.0 * H)) {
			perturbedConstituent = ctuent;
			break;
		}
	}

	if (perturbedConstituent == nullptr)
		return RetCode::E_BUFFER_CAPACITY_UNSOLVABLE;

	try {
		return calculatepHResponse(bufferCapacity, H, corrections, chemSystem, calcProps, analyticalConcentrations, perturbedConstituent);
	} catch (std::bad_alloc &) {
		return RetCode::E_NO_MEMORY;
	}
}

RetCode ECHMET_CC calculateFirstConcentrationDerivatives(RealVec *&derivatives, ECHMETReal &conductivityDerivative, const ECHMETReal &H, const NonidealityCorrections corrections, const SysComp::ChemicalSystem &chemSystem, const RealVec *analyticalConcentrations, const SysComp::Constituent *perturbedConstituent, const ECHMETReal &inIonicStrength) noexcept
{
	Solver *solver;
	RetCode tRet = prepareDerivatorContextInternal(derivatives, solver, chemSystem, corrections, false);
	if (tRet != RetCode::OK)
		return tRet;

	tRet = calculateFirstConcentrationDerivatives_prepared(derivatives, conductivityDerivative, solver, H, corrections, chemSystem, analyticalConcentrations, perturbedConstituent, inIonicStrength);

	auto ctx = solver->context();
	solver->destroy();
	ctx->destroy();

	return tRet;
}

RetCode ECHMET_CC calculateFirstConcentrationDerivatives_prepared(RealVec *derivatives, ECHMETReal &conductivityDerivative, Solver *solver, const ECHMETReal &H, const NonidealityCorrections corrections, const SysComp::ChemicalSystem &chemSystem, const RealVec *analyticalConcentrations, const SysComp::Constituent *perturbedConstituent, const ECHMETReal &inIonicStrength) noexcept
{
	typedef std::function<RetCode (RealVec *, const ECHMETReal &, SolverImpl<mpfr::mpreal, InstructionSet::GENERIC, true> *, const SysComp::ChemicalSystem &, const RealVec *, const MPRealVecWrap &, ECHMETReal &, const SysComp::Constituent *)> ExFunc;

	ExFunc executor = [corrections, &inIonicStrength](RealVec *derivatives, const ECHMETReal &H, SolverImpl<mpfr::mpreal, InstructionSet::GENERIC, true> *solver, const SysComp::ChemicalSystem &chemSystem, const RealVec *analyticalConcentrations, const MPRealVecWrap &estimatedConcentrations,
			     ECHMETReal &conductivityDerivative, const SysComp::Constituent *perturbedConstituent) {
		const size_t N = chemSystem.ionicForms->size();
		const size_t anCIdx = perturbedConstituent->analyticalConcentrationIndex;
		const mpfr::mpreal HH{H};
		RetCode tRet;
		mpfr::mpreal conductivityLow;
		mpfr::mpreal conductivityHigh;
		SolverVector<mpfr::mpreal> perturbedLow{N};
		SolverVector<mpfr::mpreal> perturbedHigh{N};

		auto solvePerturbed = [corrections, &inIonicStrength](SolverVector<mpfr::mpreal> &perturbedConcentrations, mpfr::mpreal &perturbedConductivity, SolverImpl<mpfr::mpreal, InstructionSet::GENERIC, true> *solver, const SolverVector<mpfr::mpreal> *inConcentrations, const SysComp::ChemicalSystem &chemSystem, const RealVec *acVec, const MPRealVecWrap &estimatedConcentrations) {
			mpfr::mpreal ionicStrength = mpfr::mpreal(inIonicStrength);
			SysComp::CalculatedProperties calcProps;


			std::vector<mpfr::mpreal> ics;
			try {
				ics.resize(perturbedConcentrations.rows());
			} catch (std::bad_alloc &) {
				return RetCode::E_NO_MEMORY;
			}

			RetCode tRet = SysComp::initializeCalculatedProperties(calcProps, chemSystem);
			if (tRet != RetCode::OK)
				return tRet;

			tRet = solver->solveRaw(perturbedConcentrations, ionicStrength, inConcentrations, estimatedConcentrations.get(), 2000);
			if (tRet != RetCode::OK) {
				releaseCalculatedProperties(calcProps);

				return tRet;
			}

			for (int row = 0; row < perturbedConcentrations.rows(); row++)
				ics[row] = perturbedConcentrations(row);
			calcProps.ionicStrength = CAESRealToECHMETReal(ionicStrength);

			IonProps::ComputationContext *ionCtx = IonProps::makeComputationContextExtended<mpfr::mpreal>(chemSystem, ics, IonProps::ComputationContext::NONE);
			if (ionCtx == nullptr) {
				releaseCalculatedProperties(calcProps);

				return RetCode::E_NO_MEMORY;
			}

			tRet = IonProps::correctMobilitiesInternal<mpfr::mpreal>(ionCtx, corrections, acVec, calcProps);
			if (tRet != RetCode::OK) {
				releaseCalculatedProperties(calcProps);
				ionCtx->destroy();

				return tRet;
			}

			perturbedConductivity = IonProps::calculateConductivityInternal<mpfr::mpreal>(ionCtx, calcProps);
			releaseCalculatedProperties(calcProps);
			ionCtx->destroy();

			return tRet;
		};

		SolverVector<mpfr::mpreal> inConcentrations{analyticalConcentrations->size()};
		for (size_t idx = 0; idx < analyticalConcentrations->size(); idx++)
			inConcentrations(idx) = analyticalConcentrations->at(idx);

		/* Calculate 2nd-order numerical derivative of ionic concentrations */
		inConcentrations(anCIdx) -= HH;
		tRet = solvePerturbed(perturbedLow, conductivityLow, solver, &inConcentrations, chemSystem, analyticalConcentrations, estimatedConcentrations);
		if (tRet != RetCode::OK)
			return tRet;

		inConcentrations(anCIdx) += 2.0 * HH;
		tRet = solvePerturbed(perturbedHigh, conductivityHigh, solver, &inConcentrations, chemSystem, analyticalConcentrations, estimatedConcentrations);
		if (tRet != RetCode::OK)
			return tRet;

		for (size_t idx = 0; idx < N; idx++) {
			const mpfr::mpreal der = firstDerivativeCalculator2O(perturbedLow(idx), perturbedHigh(idx), HH);

			(*derivatives)[idx] = CAESRealToECHMETReal(der);
		}

		/* DEBUG */
		ECHMET_DEBUG_CODE(
			fprintf(stderr, "cL %s\ncH %s\n", conductivityLow.toString().c_str(),
							  conductivityHigh.toString().c_str())
		)

		conductivityDerivative = CAESRealToECHMETReal(firstDerivativeCalculator2O(conductivityLow, conductivityHigh, HH));

		return tRet;
	};

	return derivatorSkin<InstructionSet::GENERIC, ECHMETReal &, const SysComp::Constituent *>(derivatives, H, solver, chemSystem, analyticalConcentrations, inIonicStrength, executor, conductivityDerivative, perturbedConstituent);
}

RetCode ECHMET_CC calculateCrossConcentrationDerivatives(RealVec *&derivatives, const ECHMETReal &H, const NonidealityCorrections corrections, const SysComp::ChemicalSystem &chemSystem, const RealVec *analyticalConcentrations, const SysComp::Constituent *perturbedConstituentJ, const SysComp::Constituent *perturbedConstituentK, const ECHMETReal &inIonicStrength) noexcept
{
	Solver *solver;
	RetCode tRet = prepareDerivatorContextInternal(derivatives, solver, chemSystem, corrections, false);
	if (tRet != RetCode::OK)
		return tRet;

	tRet = calculateCrossConcentrationDerivatives_prepared(derivatives, solver, H, chemSystem, analyticalConcentrations, perturbedConstituentJ, perturbedConstituentK, inIonicStrength);

	auto ctx = solver->context();
	solver->destroy();
	ctx->destroy();

	return tRet;
}

RetCode ECHMET_CC calculateCrossConcentrationDerivatives_prepared(RealVec *derivatives, Solver *solver, const ECHMETReal &H, const SysComp::ChemicalSystem &chemSystem, const RealVec *analyticalConcentrations, const SysComp::Constituent *perturbedConstituentJ, const SysComp::Constituent *perturbedConstituentK, const ECHMETReal &inIonicStrength) noexcept
{
	if (*(perturbedConstituentJ->name) == *(perturbedConstituentK->name))
		return calculateSecondConcentrationDerivatives(derivatives, solver, H, chemSystem, analyticalConcentrations, perturbedConstituentJ, inIonicStrength);
	else
		return calculateMixedConcentrationDerivatives(derivatives, solver, H, chemSystem, analyticalConcentrations, perturbedConstituentJ, perturbedConstituentK, inIonicStrength);
}

static
RetCode calculateSecondConcentrationDerivatives(RealVec *derivatives, Solver *solver, const ECHMETReal &H, const SysComp::ChemicalSystem &chemSystem, const RealVec *analyticalConcentrations, const SysComp::Constituent *perturbedConstituent, const ECHMETReal &inIonicStrength)
{
	typedef std::function<RetCode (RealVec *, const ECHMETReal &, SolverImpl<mpfr::mpreal, InstructionSet::GENERIC, true> *, const SysComp::ChemicalSystem &, const RealVec *, const MPRealVecWrap &, const SysComp::Constituent *)> ExFunc;

	ExFunc executor = [&inIonicStrength](RealVec *derivatives, const ECHMETReal &H, SolverImpl<mpfr::mpreal, InstructionSet::GENERIC, true> *solver, const SysComp::ChemicalSystem &chemSystem, const RealVec *analyticalConcentrations, const MPRealVecWrap &estimatedConcentrations,
			     const SysComp::Constituent *perturbedConstituent) {
		const size_t N = chemSystem.ionicForms->size();
		const size_t anCIdx = perturbedConstituent->analyticalConcentrationIndex;
		const mpfr::mpreal HH(H);
		RetCode tRet;
		SolverVector<mpfr::mpreal> perturbedLow(N);
		SolverVector<mpfr::mpreal> perturbedHigh(N);
		SolverVector<mpfr::mpreal> center(N);

		auto solvePerturbed = [&inIonicStrength](SolverVector<mpfr::mpreal> &perturbedConcentrations, SolverImpl<mpfr::mpreal, InstructionSet::GENERIC, true> *solver, const SolverVector<mpfr::mpreal> *inConcentrations,
					 const MPRealVecWrap &estimatedConcentrations) {
			 mpfr::mpreal ionicStrength = mpfr::mpreal(inIonicStrength);

			return solver->solveRaw(perturbedConcentrations, ionicStrength, inConcentrations, estimatedConcentrations.get(), 1000);
		};

		SolverVector<mpfr::mpreal> inConcentrations(analyticalConcentrations->size());
		for (size_t idx = 0; idx < analyticalConcentrations->size(); idx++)
			inConcentrations(idx) = analyticalConcentrations->at(idx);

		/* Re-solve the system again to have all concenctrations with the same accuracy */
		tRet = solvePerturbed(center, solver, &inConcentrations, estimatedConcentrations);
		if (tRet != RetCode::OK)
			return tRet;

		inConcentrations(anCIdx) -= HH;
		tRet = solvePerturbed(perturbedLow, solver, &inConcentrations, estimatedConcentrations);
		if (tRet != RetCode::OK)
			return tRet;

		inConcentrations(anCIdx) += 2.0 * HH;
		tRet = solvePerturbed(perturbedHigh, solver, &inConcentrations, estimatedConcentrations);
		if (tRet != RetCode::OK)
			return tRet;

		for (size_t idx = 0; idx < N; idx++) {
			const mpfr::mpreal der = secondDerivativeCalculator2O(perturbedLow(idx), center(idx), perturbedHigh(idx), HH);

			(*derivatives)[idx] = CAESRealToECHMETReal(der);
		}

		return tRet;
	};

	return derivatorSkin(derivatives, H, solver, chemSystem, analyticalConcentrations, inIonicStrength, executor, perturbedConstituent);
}

static
RetCode calculateMixedConcentrationDerivatives(RealVec *derivatives, Solver *solver, const ECHMETReal &H, const SysComp::ChemicalSystem &chemSystem, const RealVec *analyticalConcentrations, const SysComp::Constituent *perturbedConstituentJ, const SysComp::Constituent *perturbedConstituentK, const ECHMETReal &inIonicStrength)
{
	typedef std::function<RetCode (RealVec *, const ECHMETReal &, SolverImpl<mpfr::mpreal, InstructionSet::GENERIC, true> *, const SysComp::ChemicalSystem &, const RealVec *, const MPRealVecWrap &, const SysComp::Constituent *, const SysComp::Constituent *)> ExFunc;

	ExFunc executor = [&inIonicStrength](RealVec *derivatives, const ECHMETReal &H, SolverImpl<mpfr::mpreal, InstructionSet::GENERIC, true> *solver, const SysComp::ChemicalSystem &chemSystem, const RealVec *analyticalConcentrations, const MPRealVecWrap &estimatedConcentrations,
			     const SysComp::Constituent *perturbedConstituentJ, const SysComp::Constituent *perturbedConstituentK) {
		const size_t N = chemSystem.ionicForms->size();
		const size_t anCIdxJ = perturbedConstituentJ->analyticalConcentrationIndex;
		const size_t anCIdxK = perturbedConstituentK->analyticalConcentrationIndex;
		const mpfr::mpreal HH(H);
		RetCode tRet;
		SolverVector<mpfr::mpreal> perturbedLjLk(N);
		SolverVector<mpfr::mpreal> perturbedLjHk(N);
		SolverVector<mpfr::mpreal> perturbedHjHk(N);
		SolverVector<mpfr::mpreal> perturbedHjLk(N);

		auto solvePerturbed = [&inIonicStrength](SolverVector<mpfr::mpreal> &perturbedConcentrations, SolverImpl<mpfr::mpreal, InstructionSet::GENERIC, true> *solver, const SolverVector<mpfr::mpreal> *inConcentrations,
					 const MPRealVecWrap &estimatedConcentrations) {
			mpfr::mpreal ionicStrength = mpfr::mpreal(inIonicStrength);

			return solver->solveRaw(perturbedConcentrations, ionicStrength, inConcentrations, estimatedConcentrations.get(), 1000);
		};

		SolverVector<mpfr::mpreal> inConcentrations(analyticalConcentrations->size());
		for (size_t idx = 0; idx < analyticalConcentrations->size(); idx++)
			inConcentrations(idx) = analyticalConcentrations->at(idx);

		/* Low J, Low K */
		inConcentrations(anCIdxJ) -= HH;
		inConcentrations(anCIdxK) -= HH;
		tRet = solvePerturbed(perturbedLjLk, solver, &inConcentrations, estimatedConcentrations);
		if (tRet != RetCode::OK)
			return tRet;

		/* Low J, High K */
		inConcentrations(anCIdxK) += 2.0 * HH;
		tRet = solvePerturbed(perturbedLjHk, solver, &inConcentrations, estimatedConcentrations);
		if (tRet != RetCode::OK)
			return tRet;

		/* High J, High K */
		inConcentrations(anCIdxJ) += 2.0 * HH;
		tRet = solvePerturbed(perturbedHjHk, solver, &inConcentrations, estimatedConcentrations);
		if (tRet != RetCode::OK)
			return tRet;

		/* High J, Low K */
		inConcentrations(anCIdxK) -= 2.0 * HH;
		tRet = solvePerturbed(perturbedHjLk, solver, &inConcentrations, estimatedConcentrations);
		if (tRet != RetCode::OK)
			return tRet;

		for (size_t idx = 0; idx < N; idx++) {
			const mpfr::mpreal der = mixedDerivativeCalculator2O(perturbedLjLk(idx), perturbedLjHk(idx), perturbedHjLk(idx), perturbedHjHk(idx), HH, HH);

			(*derivatives)[idx] = CAESRealToECHMETReal(der);
		}

		return tRet;
	};

	return derivatorSkin(derivatives, H, solver, chemSystem, analyticalConcentrations, inIonicStrength, executor, perturbedConstituentJ, perturbedConstituentK);
}

RetCode ECHMET_CC prepareDerivatorContext(RealVec *&derivatives, Solver *&solver, const SysComp::ChemicalSystem &chemSystem, const NonidealityCorrections corrections) noexcept
{
	return prepareDerivatorContextInternal(derivatives, solver, chemSystem, corrections, true);
}

DDSContext::~DDSContext() noexcept {}

} // namespace CAES

IonProps::ComputationContext::~ComputationContext() noexcept {}

} // namespace ECHMET
