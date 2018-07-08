#include "syscomp_p.h"
#include <cmath>
#include <cstring>
#include <limits>
#include <string>

#define ECHMET_IMPORT_INTERNAL
#include <echmetphchconsts.h>
#include <containers/echmetskmap_p.h>

#ifdef ECHMET_DEBUG_CODE
#include <cstdio>
#include <iostream>
#endif // ECHMET_DEBUG_CODE

#define MAKE_WATER_ION(iF, u, charge) \
	iF->nucleus = nullptr; \
	iF->name = nullptr; \
	iF->nucleusCharge = std::numeric_limits<int32_t>::min(); \
	iF->limitMobility = u; \
	iF->totalCharge = charge; \
	makeNonComplex(iF);

namespace ECHMET {
namespace SysComp {

/*!
 * Build all complex forms for all ligand groups a given constituent
 *
 * @param[in,out] c The \p Constituent
 * @param[in,out] ifVec Vector of all \p IonicForm s in the chemical system
 * @param[in] ic Input representation of the constituent
 * @param[in] ligandsVec Vector of all ligands in the chemical system
 *
 * @retval RetCode::OK Success
 * @retval RetCode::E_NO_MEMORY Insufficient memory to perform the initialization
 * @retval RetCode::E_INVALID_COMPLEXATION The system contains nonsensical complexation relation
 */
static
RetCode buildComplexForms(Constituent *c, IonicFormVec *ifVec, const InConstituent &ic, const ConstituentVec *ligandsVec) noexcept
{
	RetCode tRet;

	for (int charge = c->chargeLow; charge <= c->chargeHigh; charge++) {
		const InComplexForm *inCF = nullptr;
		const ECHMETReal limitMobility = ic.mobilities->at(charge - ic.chargeLow);
		IonicForm *baseIF;

		if (ic.complexForms == nullptr)
			return RetCode::E_INVALID_CONSTITUENT;

		tRet = makeIonicForm(&baseIF, c, charge, limitMobility, ic.name);
		if (tRet != RetCode::OK)
			return tRet;

		try {
			ifVec->push_back(baseIF);
		} catch (std::bad_alloc &) {
			delete baseIF;
		}


		/* Add the base ionic form first */
		tRet = baseIF->nucleus->ionicForms->push_back(baseIF);
		if (tRet != RetCode::OK)
			return tRet;

		for (size_t cfIdx = 0; cfIdx < ic.complexForms->size(); cfIdx++) {
			const int _charge = ic.complexForms->at(cfIdx).nucleusCharge;

			if (_charge < c->chargeLow || _charge > c->chargeHigh)
				return RetCode::E_INVALID_COMPLEXATION;

			if (_charge == charge) {
				inCF = &ic.complexForms->at(cfIdx);
				break;
			}
		}

		if (inCF == nullptr)
			continue; /* This constituent does not form any complexes at this charge */

		/* Fork here for all ligand groups */
		for (size_t lgIdx = 0; lgIdx < inCF->ligandGroups->size(); lgIdx++) {
			int32_t totalComplexes = 0;
			const InLigandGroup *lg = &inCF->ligandGroups->at(lgIdx);
			IonicFormVec *groupIFVec = createECHMETVec<IonicForm *, false>(0);

			calculateMaximumVariants(0, lg->ligands, totalComplexes, 1);

			ECHMET_DEBUG_CODE(fprintf(stderr, "Maximum variants for constituent %s(%d): %u\n", c->name->c_str(), charge, totalComplexes));

			try {
				tRet = generateComplexForms(baseIF, groupIFVec, ifVec, lg->ligands, 0, totalComplexes, ligandsVec);

				if (tRet != RetCode::OK) {
					groupIFVec->destroy();

					return tRet;
				}
			} catch (std::bad_alloc &) {
				groupIFVec->destroy();

				return RetCode::E_NO_MEMORY;
			}

			tRet = baseIF->nucleus->ionicForms->append_vector(groupIFVec);
			if (tRet != RetCode::OK) {
				groupIFVec->destroy();

				return tRet;
			}
			groupIFVec->destroy();
		}
	}

	return RetCode::OK;
}

/*!
 * Builds vectors of all \p Constituent s \p and IonicForm s
 *
 * @param[in,out] cVec Vector of all \p Constituent s
 * @param[in,out] ifVec Vector of all \p IonicForm s
 * @param[in] Input description of the system
 *
 * @retval RetCode::OK Success
 * @retval RetCode::E_NO_MEMORY Insufficient memory to perform the initialization
 * @retval RetCode::E_INVALID_ARGUMENT Invalid object was passed as \p InConstituentVec
 * @retval RetCode::E_INVALID_CONSTITUENT The system contains a constituent with invalid properties
 * @retval RetCode::E_DATA_TOO_LARGE The system is too large to proces
 * @retval RetCode::E_INVALID_COMPLEXATION The system contains nonsensical complexation relation
 * @retval RetCode::E_DUPLICIT_CONSTITUENTS The system contains multiple constituents with the same name
 */
static
RetCode buildConstituentVec(ConstituentVec *cVec, IonicFormVec *ifVec, const InConstituentVec *inputData) noexcept
{
	RetCode tRet;

	const VecImpl<InConstituent, false> *dataImpl = dynamic_cast<const VecImpl<InConstituent, false> *>(inputData);
	if (dataImpl == nullptr)
		return RetCode::E_INVALID_ARGUMENT;

	if (dataImpl->size() == SIZE_MAX)
		return RetCode::E_DATA_TOO_LARGE;

	IonicFormVec *ligandIFVec = createECHMETVec<IonicForm *, false>(0);
	if (ligandIFVec == nullptr)
		return RetCode::E_NO_MEMORY;

	ConstituentVec *ligandsVec = createECHMETVec<Constituent *, false>(0);
	if (ligandsVec == nullptr) {
		ligandIFVec->destroy();

		return RetCode::E_NO_MEMORY;
	}

	/* This is just a more or less dumb conversion from the input
	 * data format to the internal representation of all constituents.
	 * All the magic happens later */
	auto processConstituentCommon = [&](Constituent *&c, const InConstituent &ic) {
		RetCode tRet;
		size_t pKasCount;
		size_t mobilitiesCount;

		pKasCount = ic.chargeHigh - ic.chargeLow;
		mobilitiesCount = pKasCount + 1;

		/* Pedantic sanity checks */
		if (pKasCount != ic.pKas->size())
			return RetCode::E_INVALID_CONSTITUENT;
		if (mobilitiesCount != ic.mobilities->size())
			return RetCode::E_INVALID_CONSTITUENT;
		if (ic.name->length() < 1)
			return RetCode::E_INVALID_CONSTITUENT;
		/* Check that a form with non-zero charge has a non-zero mobility
		 * and vice versa */
		for (int charge = ic.chargeLow; charge <= ic.chargeHigh; charge++) {
			const ECHMETReal mobility = ic.mobilities->at(charge - ic.chargeLow);
			if (mobility == 0.0 && charge != 0)
				return RetCode::E_INVALID_CONSTITUENT;
			if (mobility != 0.0 && charge == 0)
				return RetCode::E_INVALID_CONSTITUENT;
		}
		/* Negative viscosity coefficients are not allowed */
		if (ic.viscosityCoefficient < 0.0)
			return RetCode::E_INVALID_CONSTITUENT;

		try {
			c = new Constituent();
			c->ionicForms = nullptr;
			c->pKas = nullptr;
			c->limitMobilities = nullptr;
		} catch (std::bad_alloc &) {
			releaseConstituentVec(ligandsVec);
			releaseIonicFormVec(ligandIFVec);
			ligandsVec->destroy();
			ligandIFVec->destroy();

			return RetCode::E_NO_MEMORY;
		}

		try {
			c->ionicForms = createECHMETVec<IonicForm *, false>(0);
			c->pKas = createRealVec(pKasCount);
			c->limitMobilities = createRealVec(mobilitiesCount);
		} catch (std::bad_alloc &) {
			releaseConstituentVec(ligandsVec);
			releaseIonicFormVec(ligandIFVec);
			ligandsVec->destroy();
			ligandIFVec->destroy();
			c->limitMobilities->destroy();
			c->pKas->destroy();
			c->ionicForms->destroy();
			delete c;

			return RetCode::E_NO_MEMORY;
		} catch (std::length_error &) {
			releaseConstituentVec(ligandsVec);
			releaseIonicFormVec(ligandIFVec);
			ligandsVec->destroy();
			ligandIFVec->destroy();
			c->limitMobilities->destroy();
			c->pKas->destroy();
			c->ionicForms->destroy();
			delete c;

			return RetCode::E_DATA_TOO_LARGE;
		}

		/* Naïve item-by-item vector copying */
		for (const ECHMETReal &d : static_cast<const VecImpl<ECHMETReal, false> *>(ic.pKas)->STL()) {
			tRet = c->pKas->push_back(d);
			if (tRet != RetCode::OK) {
				releaseConstituentVec(ligandsVec);
				releaseIonicFormVec(ligandIFVec);
				ligandsVec->destroy();
				ligandIFVec->destroy();
				c->limitMobilities->destroy();
				c->pKas->destroy();
				c->ionicForms->destroy();
				delete c;

				return tRet;
			}

		}
		for (const ECHMETReal &d : static_cast<const VecImpl<ECHMETReal, false> *>(ic.mobilities)->STL()) {
			tRet = c->limitMobilities->push_back(d);
			if (tRet != RetCode::OK) {
				releaseConstituentVec(ligandsVec);
				releaseIonicFormVec(ligandIFVec);
				ligandsVec->destroy();
				ligandIFVec->destroy();
				c->limitMobilities->destroy();
				c->pKas->destroy();
				c->ionicForms->destroy();
				delete c;

				return tRet;
			}
		}

		c->name = createFixedString(ic.name->c_str());
		if (c->name == nullptr) {
			releaseConstituentVec(ligandsVec);
			releaseIonicFormVec(ligandIFVec);
			ligandsVec->destroy();
			ligandIFVec->destroy();
			c->limitMobilities->destroy();
			c->pKas->destroy();
			c->ionicForms->destroy();
			delete c;

			return RetCode::E_NO_MEMORY;
		}
		if (isConstituentContained(cVec, ligandsVec, c->name)) {
			releaseConstituentVec(ligandsVec);
			releaseIonicFormVec(ligandIFVec);
			ligandsVec->destroy();
			ligandIFVec->destroy();
			c->limitMobilities->destroy();
			c->pKas->destroy();
			c->ionicForms->destroy();
			c->name->destroy();
			delete c;

			return RetCode::E_DUPLICIT_CONSTITUENTS;
		}

		c->ctype = ic.ctype;
		c->chargeLow = ic.chargeLow;
		c->chargeHigh = ic.chargeHigh;
		c->viscosityCoefficient = ic.viscosityCoefficient;

		return RetCode::OK;
	};

	/* We need to know about all ligands in advance,
	 * process ligands first and store them in a separate vector */
	for (const InConstituent &ic : dataImpl->STL()) {
		Constituent *c;

		switch (ic.ctype) {
		case ConstituentType::INVALID:
			return RetCode::E_INVALID_CONSTITUENT;
		case ConstituentType::NUCLEUS:
			continue; /* Ignore nuclei in this pass */
		case ConstituentType::LIGAND:
			break;
		}

		tRet = processConstituentCommon(c, ic);
		if (tRet != RetCode::OK)
			return tRet;

		tRet = buildLigandIonicForms(c, ligandIFVec, ic.name);
		if (tRet != RetCode::OK) {
			releaseConstituentVec(ligandsVec);
			releaseIonicFormVec(ligandIFVec);
			ligandsVec->destroy();
			ligandIFVec->destroy();
			c->limitMobilities->destroy();
			c->pKas->destroy();
			c->ionicForms->destroy();
			c->name->destroy();
			delete c;

			return tRet;
		}

		tRet = ligandsVec->push_back(c);
		if (tRet != RetCode::OK) {
			releaseConstituentVec(ligandsVec);
			releaseIonicFormVec(ligandIFVec);
			ligandsVec->destroy();
			ligandIFVec->destroy();
			c->limitMobilities->destroy();
			c->pKas->destroy();
			c->ionicForms->destroy();
			c->name->destroy();
			delete c;

			return tRet;
		}
	}

	/* The ligands are processed, do the nuclei now */
	for (const InConstituent &ic : dataImpl->STL()) {
		Constituent *c;

		switch (ic.ctype) {
		case ConstituentType::INVALID:
			return RetCode::E_INVALID_CONSTITUENT;
		case ConstituentType::NUCLEUS:
			break;
		case ConstituentType::LIGAND:
			continue; /* Ignore ligands in this pass */
		}

		tRet = processConstituentCommon(c, ic);
		if (tRet != RetCode::OK)
			return tRet;

		/* The basic constituent is built now.
		 * If the constituent is a complex nucleus, build all complex forms for it */
		tRet = buildComplexForms(c, ifVec, ic, ligandsVec);
		if (tRet != RetCode::OK) {
			releaseConstituentVec(ligandsVec);
			releaseIonicFormVec(ligandIFVec);
			ligandsVec->destroy();
			ligandIFVec->destroy();
			c->limitMobilities->destroy();
			c->pKas->destroy();
			c->ionicForms->destroy();
			c->name->destroy();
			delete c;

			return tRet;
		}

		tRet = cVec->push_back(c);
		if (tRet != RetCode::OK) {
			releaseConstituentVec(ligandsVec);
			releaseIonicFormVec(ligandIFVec);
			ligandsVec->destroy();
			ligandIFVec->destroy();
			c->limitMobilities->destroy();
			c->pKas->destroy();
			c->ionicForms->destroy();
			c->name->destroy();
			delete c;
		}
	}

	tRet = ifVec->append_vector(ligandIFVec);
	if (tRet != RetCode::OK) {
		releaseConstituentVec(ligandsVec);
		releaseIonicFormVec(ligandIFVec);
		ligandsVec->destroy();
		ligandIFVec->destroy();

		return tRet;
	}

	if (tRet != RetCode::OK) {
		releaseConstituentVec(ligandsVec);
		releaseIonicFormVec(ligandIFVec);
		ligandsVec->destroy();
		ligandIFVec->destroy();

		return tRet;
	}
	tRet = cVec->append_vector(ligandsVec);
	if (tRet != RetCode::OK) {
		releaseConstituentVec(ligandsVec);
		releaseIonicFormVec(ligandIFVec);
		ligandsVec->destroy();
		ligandIFVec->destroy();

		return tRet;
	}

	/* Add H+ ahd OH- to the ionicConcentrations vector */
	{
		IonicForm *hydroxoniumIF;
		IonicForm *hydroxyleIF;

		try {
			hydroxyleIF = new IonicForm();
			MAKE_WATER_ION(hydroxyleIF, PhChConsts::mobilityOH, -1);
			hydroxyleIF->ifType = IonicFormType::OH;
		} catch (std::bad_alloc &) {
			releaseConstituentVec(ligandsVec);
			releaseIonicFormVec(ligandIFVec);
			ligandsVec->destroy();
			ligandIFVec->destroy();

			return RetCode::E_NO_MEMORY;
		}

		tRet = ifVec->push_front(hydroxyleIF);
		if (tRet != RetCode::OK) {
			delete hydroxyleIF;

			releaseConstituentVec(ligandsVec);
			releaseIonicFormVec(ligandIFVec);
			ligandsVec->destroy();
			ligandIFVec->destroy();

			return tRet;
		}

		try {
			hydroxoniumIF = new IonicForm();
			MAKE_WATER_ION(hydroxoniumIF, PhChConsts::mobilityH3O, 1);
			hydroxoniumIF->ifType = IonicFormType::H;
		} catch (std::bad_alloc &) {
			releaseConstituentVec(ligandsVec);
			releaseIonicFormVec(ligandIFVec);
			ligandsVec->destroy();
			ligandIFVec->destroy();

			return RetCode::E_NO_MEMORY;
		}

		tRet = ifVec->push_front(hydroxoniumIF);
		if (tRet != RetCode::OK) {
			delete hydroxoniumIF;

			releaseConstituentVec(ligandsVec);
			releaseIonicFormVec(ligandIFVec);
			ligandsVec->destroy();
			ligandIFVec->destroy();

			return tRet;
		}
	}

	ligandIFVec->destroy();
	ligandsVec->destroy();

	return tRet;
}

/*!
 * Builds all ionic forms for a given ligand
 *
 * @param[in] c Constituent representing the ligand
 * @param[in, out] ligandIFVec Vector of \p IonicForm s of ligands
 * @param[in] name Ligand base name
 *
 * @retval RetCode::OK Success
 * @retval RetCode::E_NO_MEMORY Insufficient memory to create ionic forms
 */
static
RetCode buildLigandIonicForms(const Constituent *c, IonicFormVec *ligandIFVec, const FixedString *name) noexcept
{
	for (int charge = c->chargeLow; charge <= c->chargeHigh; charge++) {
		const ECHMETReal limitMobility = c->limitMobilities->at(charge - c->chargeLow);
		IonicForm *lIF;
		RetCode tRet = makeIonicForm(&lIF, c, charge, limitMobility, name);

		if (tRet != RetCode::OK) {
			delete lIF;

			return tRet;
		}

		tRet = ligandIFVec->push_back(lIF);
		if (tRet != RetCode::OK) {
			delete lIF;

			return RetCode::E_NO_MEMORY;
		}

		tRet = c->ionicForms->push_back(lIF);
		if (tRet != RetCode::OK)
			return RetCode::E_NO_MEMORY;
	}

	return RetCode::OK;
}

/*!
 * Calculates the maximum number of complex form that a given complexing
 * component can form with its ligands. Note that the procedure may call itself recursively.
 *
 * @param[in] i Index in the vector of all ligand ionic forms for the given complexing component to start from
 * @param[in] ligandIFs Vector of all ligand ionic forms for the given complexing component
 * @param[in,out] total Total number of forms counted so far. This must be initialized to zero by the first callee.
 * @param[in] accum Total number of forms propagating through one branch of recursive calls
 */
static
void calculateMaximumVariants(const size_t i, const InLFVec *ligandIFs, int32_t &total, int32_t accum) noexcept
{
	const size_t N = ligandIFs->size();

	for (size_t j = i; j < N; j++) {
		uint32_t m = ligandIFs->at(j).maxCount;

		total += m * accum;
		if (N > 1)
			calculateMaximumVariants(j + 1, ligandIFs, total, (m * accum));
	}
}

/*!
 * Deep-copies \p InComplexForms vector of \p InConstituent
 *
 * @param[in] src \p InComplexForms vector to be copied.
 *
 * @return Pointer to the duplicated vector, \p NULL on failure.
 */
static
InCFVec * duplicateComplexForms(const InCFVec *src) noexcept
{
	const auto releaseLigandForm = [](InLigandForm &lf) {
		lf.ligandName->destroy();
		lf.pBs->destroy();
		lf.mobilities->destroy();
	};

	const auto releaseLigandsVec = [&releaseLigandForm](InLFVec *vec) {
		for (size_t idx = 0; idx < vec->size(); idx++)
			releaseLigandForm(vec->operator[](idx));

		vec->destroy();
	};

	const auto releaseLigandGroupsVec = [&releaseLigandsVec](InLGVec *vec) {
		for (size_t idx = 0; idx < vec->size(); idx++) {
			const auto &lgGroup = vec->operator[](idx);

			releaseLigandsVec(lgGroup.ligands);
		}

		vec->destroy();
	};

	const auto releaseComplexForms = [&releaseLigandGroupsVec](InCFVec *cfVec) {
		for (size_t idx = 0; idx < cfVec->size(); idx++)
			releaseLigandGroupsVec(cfVec->operator[](idx).ligandGroups);
	};

	InCFVec *dupCfVec = createInCFVec(src->size());

	if (dupCfVec == nullptr)
		return nullptr;

	for (size_t idx = 0; idx < src->size(); idx++) {
		const auto &cForm = src->at(idx);

		InComplexForm dupCForm;
		InLGVec *dupLgVec = createInLGVec(cForm.ligandGroups->size());
		if (dupLgVec == nullptr)
			goto err_out;

		dupCForm.nucleusCharge = cForm.nucleusCharge;

		/* Copy ligand groups */
		size_t jdx;
		for (jdx = 0; jdx < cForm.ligandGroups->size(); jdx++) {
			const auto lgGrp = cForm.ligandGroups->at(jdx);

			InLFVec *dupLigands = createInLFVec(lgGrp.ligands->size());
			if (dupLigands == nullptr)
				break;

			size_t kdx;
			for (kdx = 0; kdx < lgGrp.ligands->size(); kdx++) {
				const auto &lf = lgGrp.ligands->at(kdx);
				InLigandForm dupLf;

				dupLf.ligandName = createFixedString(lf.ligandName->c_str());
				if (dupLf.ligandName == nullptr)
					break;

				dupLf.pBs = lf.pBs->duplicate();
				if (dupLf.pBs == nullptr) {
					dupLf.ligandName->destroy();

					break;
				}

				dupLf.mobilities = lf.mobilities->duplicate();
				if (dupLf.mobilities == nullptr) {
					dupLf.pBs->destroy();
					dupLf.ligandName->destroy();

					break;
				}

				if (dupLigands->push_back(dupLf) != RetCode::OK) {
					releaseLigandForm(dupLf);

					break;
				}
			}
			if (kdx < lgGrp.ligands->size()) { /* Not all ligand forms were duplicated successfully */
				releaseLigandsVec(dupLigands);

				break;
			}

			if (dupLgVec->push_back({ dupLigands }) != RetCode::OK) {
				releaseLigandsVec(dupLigands);

				break;
			}
		}
		if (jdx < cForm.ligandGroups->size()) { /* Not all ligand groups were duplicated successfully */
			releaseLigandGroupsVec(dupLgVec);

			goto err_out;
		}

		dupCForm.ligandGroups = dupLgVec;

		if (dupCfVec->push_back(dupCForm) != RetCode::OK) {
			releaseLigandGroupsVec(dupLgVec);

			goto err_out;
		}
	}

	return dupCfVec;

err_out:
	releaseComplexForms(dupCfVec);

	return nullptr;
}

/*!
 * Generates all complex forms for one ligand group of a given Nucleus
 *
 * @param[in] baseIF Ionic form acting as the Nucleus
 * @param[in] Vector of all \p IonicForm s in the entire group
 * @param[in] Vector if all \p IonicForm s in the entire chemical system
 * @param[in] Vector of all ligands to form a complex with
 * @param[in] formShift Read the code...
 * @param[in] toGenerate Number of complex ionic forms left to generate
 * @param[in] ligandsVec Vector of all \p Constituent s of Ligand type
 *
 * @retval RetCode::OK Success
 * @retval RetCode::E_NO_MEMORY Insufficient memory to generate all complex forms
 * @retval RetCode::E_INVALID_COMPLEXATION Nonsensical complexation scheme
 */
static
RetCode generateComplexForms(IonicForm *baseIF, IonicFormVec *groupIFVec, IonicFormVec *ifVec, InLFVec *ligandIFs,
			     size_t formShift, int32_t toGenerate, const ConstituentVec *ligandsVec) noexcept(false)
{
	RetCode tRet;
	const size_t n = ligandIFs->size();
	size_t lShift = 0;

	/* Find ligand shift */
	for (size_t idx = 0; idx < n; idx++) {
		const InLigandForm *lF = &ligandIFs->at(idx);
		if (isLigandFormContained(baseIF->containedLigandIFs, lF, lF->charge))
			lShift = idx + 1;
	}

	for (size_t idx = lShift; idx < n; idx++) {
		const InLigandForm *lF = &ligandIFs->at(idx);

		const std::string name = std::string(baseIF->name->c_str()) + std::string(lF->ligandName->c_str());
		Constituent *ligand = findLigand(lF, ligandsVec);

		if (ligand == nullptr) {
			ECHMET_DEBUG_CODE(fprintf(stderr, "nullptr to ligand"));
			return RetCode::E_INVALID_COMPLEXATION;
		}
		if (ligand->chargeLow > lF->charge || ligand->chargeHigh < lF->charge) {
			ECHMET_DEBUG_CODE(fprintf(stderr, "Invalid charge congifuration, expected low=%d, high=%d, got charge=%d\n", ligand->chargeLow, ligand->chargeHigh, lF->charge));
			return RetCode::E_INVALID_COMPLEXATION;
		}

		IonicForm *previousIF = nullptr;
		for (uint32_t lCnt = 1; lCnt <= lF->maxCount; lCnt++) {
			IonicForm *newIF;

			newIF = new IonicForm();
			newIF->name = createFixedString((name + "(" + std::to_string(lF->charge) + ")" + (lCnt == 1 ? "" : std::to_string(lCnt))).c_str());
			newIF->ifType = IonicFormType::CONSTITUENT;

			if (newIF->name == nullptr) {
				delete newIF;

				return RetCode::E_NO_MEMORY;
			}

			if (isIonicFormContained(baseIF->nucleus->ionicForms, newIF)) {
				newIF->name->destroy();
				delete newIF;

				return RetCode::E_INVALID_COMPLEXATION;
			}

			newIF->containedLigandIFs = createECHMETVec<ContainedLigandIonicForm, false>(0);
			if (newIF->containedLigandIFs == nullptr) {
				newIF->name->destroy();
				delete newIF;

				return RetCode::E_NO_MEMORY;
			}

			/* Check if we have a mixed complex form */
			if (baseIF->ligand == nullptr) {
				/* We have a free nucleus */
				newIF->pB = lF->pBs->at(lCnt - 1);
				newIF->limitMobility = lF->mobilities->at(lCnt - 1);
			} else {
				/* We have a mixed complex, leave the data empty and have
				 * the user enter it later */
				newIF->pB = nan("");
				newIF->limitMobility = nan("");
			}

			/* Fill out the rest of the complexation data */
			newIF->ligand = ligand;
			newIF->ligandCharge = lF->charge;
			if (previousIF == nullptr)
				newIF->ancestor = baseIF;
			else
				newIF->ancestor = previousIF;

			newIF->nucleus = newIF->ancestor->nucleus;
			newIF->nucleusCharge = baseIF->nucleusCharge;
			newIF->totalCharge = baseIF->totalCharge + (lF->charge * lCnt);
			newIF->ligandCount = lCnt;

			/* Check that a complex form with overall non-zero charge has a non-zero mobility */
			if (newIF->limitMobility == 0.0 && newIF->totalCharge != 0) {
				newIF->name->destroy();
				newIF->containedLigandIFs->destroy();
				delete newIF;

				return RetCode::E_INVALID_COMPLEXATION;
			}
			/* Check that a complex with overall zero charge has zero mobility */
			if (newIF->limitMobility != 0.0 && newIF->totalCharge == 0) {
				newIF->name->destroy();
				newIF->containedLigandIFs->destroy();
				delete newIF;

				return RetCode::E_INVALID_COMPLEXATION;
			}

			if (baseIF->containedLigandIFs != nullptr) {
				tRet = newIF->containedLigandIFs->append_vector(baseIF->containedLigandIFs);
				if (tRet != RetCode::OK) {
					newIF->name->destroy();
					newIF->containedLigandIFs->destroy();
					delete newIF;

					return tRet;
				}
			}

			tRet = newIF->containedLigandIFs->push_back(makeContainedLigandIonicForm(ligand, lF->charge));
			if (tRet != RetCode::OK) {
				newIF->name->destroy();
				newIF->containedLigandIFs->destroy();
				delete newIF;

				return tRet;
			}

			tRet = groupIFVec->push_back(newIF);
			if (tRet != RetCode::OK) {
				newIF->name->destroy();
				newIF->containedLigandIFs->destroy();
				delete newIF;

				return tRet;
			}

			tRet = ifVec->push_back(newIF);
			if (tRet != RetCode::OK) {
				groupIFVec->pop_back();

				newIF->name->destroy();
				newIF->containedLigandIFs->destroy();
				delete newIF;

				return tRet;
			}

			tRet = ligand->ionicForms->push_back(newIF);
			if (tRet != RetCode::OK) {
				groupIFVec->pop_back();
				ifVec->pop_back();

				newIF->name->destroy();
				newIF->containedLigandIFs->destroy();
				delete newIF;

				return tRet;
			}

			previousIF = newIF;
			toGenerate--;
			ECHMET_DEBUG_CODE(fprintf(stderr, "baseIF: %s, lF: %s, toGenerate: %d, formShift %zu\n", baseIF->name->c_str(), newIF->name->c_str(), toGenerate, formShift));
		}
	}


	formShift++;
	if (toGenerate <= 0)
		return RetCode::OK;

	return generateComplexForms(groupIFVec->at(formShift - 1), groupIFVec, ifVec, ligandIFs, formShift, toGenerate, ligandsVec);
}

/*!
 * Tries to find a matching \p Constituent for a given \p InLigandForm
 *
 * @param[in] lF The \p InLigandForm
 * @param[in] ligandVec Vector of \p Constituent s of Ligand type
 *
 * @retval Pointer to \p Constituent if a match is found
 * @retval nullptr No match was found
 */
static
Constituent * findLigand(const InLigandForm *lF, const ConstituentVec *ligandsVec) noexcept
{
	if (lF == nullptr || ligandsVec == nullptr)
		return nullptr;

	for (size_t idx = 0; idx < ligandsVec->size(); idx++) {
		ECHMET_DEBUG_CODE(fprintf(stderr, "Target=%s, Current=%s\n", lF->ligandName->c_str(), ligandsVec->at(idx)->name->c_str()));

		FixedString *other = ligandsVec->at(idx)->name;
		if (*(lF->ligandName) == *other)
			return ligandsVec->at(idx);
	}

	return nullptr;
}

/*!
 * Initializes mapping of indices into the arrays of computed
 * properties to names of constituents and ionic forms
 *
 * @param[in] chemSystem Chemical system to create the mapping for
 */
static
void initializeChemicalSystemMapping(ChemicalSystem &chemSystem) noexcept(false)
{
	const size_t NIF = chemSystem.ionicForms->size();
	const size_t NCO = chemSystem.constituents->size();

	std::map<std::string, size_t> &icMapping = static_cast<SKMapImpl<size_t> *>(chemSystem.ionicConcentrationsByName)->STL();
	std::map<std::string, size_t> &imMapping = static_cast<SKMapImpl<size_t> *>(chemSystem.ionicMobilitiesByName)->STL();
	std::map<std::string, size_t> &emMapping = static_cast<SKMapImpl<size_t> *>(chemSystem.effectiveMobilitiesByName)->STL();
	std::map<std::string, size_t> &ancMapping = static_cast<SKMapImpl<size_t> *>(chemSystem.analyticalConcentrationsByName)->STL();

	(*chemSystem.ionicForms)[0]->ionicConcentrationIndex = 0;
	(*chemSystem.ionicForms)[0]->ionicMobilityIndex = 0;
	(*chemSystem.ionicForms)[1]->ionicConcentrationIndex = 1;
	(*chemSystem.ionicForms)[1]->ionicMobilityIndex = 1;

	for (size_t idx = 2; idx < NIF; idx++) {
		IonicForm *iF = chemSystem.ionicForms->at(idx);
		const std::string name = std::string(iF->name->c_str());

		icMapping.emplace(name, idx);
		imMapping.emplace(name, idx);
		iF->ionicConcentrationIndex = idx;
		iF->ionicMobilityIndex = idx;
	}

	for (size_t idx = 0; idx < NCO; idx++) {
		Constituent *c = chemSystem.constituents->at(idx);
		const std::string name(c->name->c_str());

		emMapping.emplace(name, idx);
		ancMapping.emplace(name, idx);
		c->analyticalConcentrationIndex = idx;
		c->effectiveMobilityIndex = idx;
	}
}

/*!
 * Checks if a given constituent is contained in a vector of \p Constituent s
 *
 * @param[in] cVec Vector of Nucleus type \p Constituent s
 * @param[in] ligandsVec Vector of Ligand type \p Constituent s
 * @param[in] newName Internal name of the constituent to be looked for
 *
 * @retval true The constituent is contained
 * @retval false The constituent is not contained
 */
static
bool isConstituentContained(const ConstituentVec *cVec, const ConstituentVec *ligandsVec, const FixedString *newName) noexcept
{
	auto _internal = [](const ConstituentVec *vec, const FixedString *name) {
		for (size_t idx = 0; idx < vec->size(); idx++) {
			const Constituent *c = vec->at(idx);

		if (*(c->name) == (*name))
			return true;
		}

		return false;
	};

	if (_internal(cVec, newName))
		return true;

	return _internal(ligandsVec, newName);
}

/*!
 * Checks if a given ionic form is contained in a vector of \p IonicForm s
 *
 * @param[in] ionicForms Vector of \p IonicForms s
 * @param[in] iF The ionic form
 *
 * @retval true The ionic form is contained
 * @retval false The ionic form is not contained
 */
static
bool isIonicFormContained(const IonicFormVec *ionicForms, const IonicForm *iF) noexcept
{
	for (size_t idx = 0; idx < ionicForms->size(); idx++) {
		const IonicForm *cIF = ionicForms->at(idx);

		if (*(cIF->name) == *(iF->name))
			return true;
	}

	return false;
}

/*!
 * Checks if a given ligand ionic form is contained in a list of \p ContainedLigandIonicForms
 *
 * @param[in] containedLigandIFs Vector of \p ContainedLigandIonicForms
 * @param[in] ligand The ligand whose presence is to be checked
 * @param[in] charge Charge of the checked ligand
 *
 * @retval true The ligand with the given charge is contained
 * @retval false The ligand iwth the given charge is not contained
 */
static
bool isLigandFormContained(const ContainedLigandIonicFormVec *containedLigandIFs, const InLigandForm *ligand, const int charge) noexcept
{
	if (containedLigandIFs == nullptr)
		return false;

	for (size_t idx = 0; idx < containedLigandIFs->size(); idx++) {
		const ContainedLigandIonicForm &cLF = containedLigandIFs->at(idx);

		if (*(cLF.ligand->name) == *(ligand->ligandName) && cLF.charge == charge)
			return true;
	}

	return false;
}

/*!
 * Creates a \p ContainedIonicForm object
 *
 * @param[in] ligand The ligand
 * @param[in] charge Charge of the ligand
 *
 * @return Initialized \p ContainedLigandIonicForm
 */
static
ContainedLigandIonicForm makeContainedLigandIonicForm(const Constituent *ligand, const int charge)
{
	ContainedLigandIonicForm cLF;
	cLF.ligand = ligand;
	cLF.charge = charge;

	return cLF;
}

/*!
 * Creates an \p IonicForm and sets its base properties
 *
 * @param[in,out] iF Pointer to the \p IonicForm to be created
 * @param[in] c \p Constituent to which the ionic form belongs to
 * @param[in] charge Electric charge of the ionic form
 * @param[in] limitMobility Limit electrophoretic mobility of the ionic form
 * @param[in] Internal base name of the \p IonicForm
 *
 * @retval RetCode::OK Success
 * @retval RetCode::E_NO_MEMORY Insufficient memory to create the ionic form
 */
static
RetCode makeIonicForm(IonicForm **iF, const Constituent *c, const int charge, const ECHMETReal &limitMobility, const FixedString *name) noexcept
{
	std::string _name;

	try {
		*iF = new IonicForm();

		_name = std::string(name->c_str());
		_name += "(" + std::to_string(charge) + ")";
	} catch (std::bad_alloc &) {
		return RetCode::E_NO_MEMORY;
	}

	(*iF)->name = createFixedString(_name.c_str());
	if ((*iF)->name == nullptr) {
		delete *iF;

		return RetCode::E_NO_MEMORY;
	}

	(*iF)->nucleus = c;
	(*iF)->nucleusCharge = charge;
	(*iF)->totalCharge = charge;
	(*iF)->limitMobility = limitMobility;
	(*iF)->ifType = IonicFormType::CONSTITUENT;
	makeNonComplex(*iF);

	return RetCode::OK;
}

/*!
 * Sets variables that describe complexation relation of \p IonicForm
 * to default values that denote no complexation
 *
 * @param[in] iF \p IonicForm to operate on
 */
static
void makeNonComplex(IonicForm *iF) noexcept
{
	iF->ligand = nullptr;
	iF->ligandCharge = std::numeric_limits<int32_t>::min();
	iF->ligandCount = 0;
	iF->ancestor = nullptr;
	iF->pB = nan("");
	iF->containedLigandIFs = nullptr;
}

/*!
 * Releases all resources claimed by \p ConstituentVec
 * including its contents
 *
 * @param[in] vec \p ConstituentVec to be released
 */
static
void releaseConstituentVec(const ConstituentVec *vec) noexcept
{
	for (size_t idx = 0; idx < vec->size(); idx++) {
		Constituent *c = vec->at(idx);

		c->name->destroy();
		c->ionicForms->destroy();
		c->pKas->destroy();
		c->limitMobilities->destroy();
		delete c;
	}
}

/*!
 * Releases all resources claimed by \p IonicFormVec
 * including its contents
 *
 * @param[in] \p IonicFormVec to be released
 */
static
void releaseIonicFormVec(const IonicFormVec *vec) noexcept
{
	for (size_t idx = 0; idx < vec->size(); idx++) {
		IonicForm *iF = vec->at(idx);

		if (iF->ifType == IonicFormType::CONSTITUENT)
			iF->name->destroy();
		if (iF->containedLigandIFs != nullptr)
			iF->containedLigandIFs->destroy();

		delete iF;
	}
}

/* Public interface functions */

InCFVec * ECHMET_CC createInCFVec(const size_t reserve) noexcept
{
	return ECHMET::createECHMETVec<InComplexForm, false>(reserve);
}

InConstituentVec * ECHMET_CC createInConstituentVec(const size_t reserve) noexcept
{
	return ECHMET::createECHMETVec<InConstituent, false>(reserve);
}

InLFVec * ECHMET_CC createInLFVec(const size_t reserve) noexcept
{
	return ECHMET::createECHMETVec<InLigandForm, false>(reserve);
}

InLGVec * ECHMET_CC createInLGVec(const size_t reserve) noexcept
{
	return ECHMET::createECHMETVec<InLigandGroup, false>(reserve);
}

RetCode ECHMET_CC initializeCalculatedProperties(CalculatedProperties &calcProps, const ChemicalSystem &chemSystem) noexcept
{
	const size_t NIF = chemSystem.ionicForms->size();
	const size_t NCO = chemSystem.constituents->size();

	calcProps.ionicConcentrations = createRealVec(NIF);
	if (calcProps.ionicConcentrations == nullptr)
		return RetCode::E_NO_MEMORY;

	calcProps.ionicMobilities = createRealVec(NIF);
	if (calcProps.ionicMobilities == nullptr) {
		calcProps.ionicConcentrations->destroy();
		return RetCode::E_NO_MEMORY;
	}

	calcProps.effectiveMobilities = createRealVec(NCO);
	if (calcProps.ionicMobilities == nullptr) {
		calcProps.ionicMobilities->destroy();
		calcProps.ionicConcentrations->destroy();
		return RetCode::E_NO_MEMORY;
	}


	try {
		calcProps.ionicConcentrations->resize(NIF);
		calcProps.ionicMobilities->resize(NIF);
		calcProps.effectiveMobilities->resize(NCO);

		for (size_t idx = 0; idx < NIF; idx++) {
			const IonicForm *iF = chemSystem.ionicForms->at(idx);

			(*calcProps.ionicMobilities)[idx] = iF->limitMobility;
		}

		return RetCode::OK;
	} catch (std::bad_alloc &) {
		calcProps.ionicConcentrations->destroy();
		calcProps.ionicMobilities->destroy();
		calcProps.effectiveMobilities->destroy();
		return RetCode::E_NO_MEMORY;
	}
}

RetCode ECHMET_CC duplicateInConstituent(InConstituent &dst, const InConstituent &src) noexcept
{
	FixedString *name;
	RealVec *pKas;
	RealVec *mobilities;
	InCFVec *complexForms = nullptr;

	name = createFixedString(src.name->c_str());
	if (name == nullptr)
		return RetCode::E_NO_MEMORY;

	pKas = src.pKas->duplicate();
	if (pKas == nullptr) {
		name->destroy();

		return RetCode::E_NO_MEMORY;
	}

	mobilities = src.mobilities->duplicate();
	if (mobilities == nullptr) {
		name->destroy();
		pKas->destroy();

		return RetCode::E_NO_MEMORY;
	}

	if (src.complexForms != nullptr) {
		complexForms = duplicateComplexForms(src.complexForms);

		if (complexForms == nullptr) {
			name->destroy();
			pKas->destroy();
			mobilities->destroy();

			return RetCode::E_NO_MEMORY;
		}
	}

	dst.ctype = src.ctype;
	dst.name = name;
	dst.chargeLow = src.chargeLow;
	dst.chargeHigh = src.chargeHigh;
	dst.pKas = pKas;
	dst.mobilities = mobilities;
	dst.complexForms = complexForms;
	dst.viscosityCoefficient = src.viscosityCoefficient;

	return RetCode::OK;
}

InConstituentVec * ECHMET_CC duplicateInConstituentVec(const InConstituentVec *src) noexcept
{
	InConstituentVec *dupVec = createInConstituentVec(src->size());

	if (dupVec == nullptr)
		return nullptr;

	for (size_t idx = 0; idx < src->size(); idx++) {
		const auto &inC = src->at(idx);

		InConstituent dupC;
		if (duplicateInConstituent(dupC, inC) != RetCode::OK) {
			releaseInputData(dupVec);

			return nullptr;
		}

		if (dupVec->push_back(dupC) != RetCode::OK) {
			releaseInConstituent(dupC);
			releaseInputData(dupVec);

			return nullptr;
		}
	}

	return dupVec;
}

RetCode ECHMET_CC makeAnalyticalConcentrationsVec(RealVec *&acVec, const ChemicalSystem &chemSystem) noexcept
{
	acVec = createECHMETVec<ECHMETReal, false>(0);

	if (acVec == nullptr)
		return RetCode::E_NO_MEMORY;

	return acVec->resize(chemSystem.constituents->size());
}

RetCode ECHMET_CC makeComposition(ChemicalSystem &chemSystem, CalculatedProperties &calcProps, const InConstituentVec *inputData) noexcept
{
	RetCode tRet;
	ConstituentVec *cVec = nullptr;
	IonicFormVec *ifVec = nullptr;
	ConcentrationSKMap *ancMap = nullptr;
	ConcentrationSKMap *icMap = nullptr;
	ConcentrationSKMap *imMap = nullptr;
	ConcentrationSKMap *emMap = nullptr;

	memset(&chemSystem, 0, sizeof(ChemicalSystem));
	memset(&calcProps, 0, sizeof(CalculatedProperties));

	try {
		cVec = ECHMET::createECHMETVec<Constituent *, false>(0);
		ifVec = ECHMET::createECHMETVec<IonicForm *, false>(0);
		ancMap = createSKMap<size_t>();
		icMap = createSKMap<size_t>();
		imMap = createSKMap<size_t>();
		emMap = createSKMap<size_t>();
	} catch (std::bad_alloc &) {
		emMap->destroy();
		imMap->destroy();
		icMap->destroy();
		ancMap->destroy();
		ifVec->destroy();
		cVec->destroy();

		memset(&chemSystem, 0, sizeof(ChemicalSystem));
		memset(&calcProps, 0, sizeof(CalculatedProperties));

		return RetCode::E_NO_MEMORY;
	}

	tRet = buildConstituentVec(cVec, ifVec, inputData);
	if (tRet != RetCode::OK) {
		emMap->destroy();
		imMap->destroy();
		icMap->destroy();
		ancMap->destroy();

		releaseConstituentVec(cVec);
		releaseIonicFormVec(ifVec);

		cVec->destroy();
		ifVec->destroy();

		memset(&chemSystem, 0, sizeof(ChemicalSystem));
		memset(&calcProps, 0, sizeof(CalculatedProperties));

		return tRet;
	}

	chemSystem.constituents = cVec;
	chemSystem.ionicForms = ifVec;
	chemSystem.analyticalConcentrationsByName = ancMap;
	chemSystem.ionicConcentrationsByName = icMap;
	chemSystem.ionicMobilitiesByName = imMap;
	chemSystem.effectiveMobilitiesByName = emMap;

	try {
		initializeChemicalSystemMapping(chemSystem);
	} catch (std::bad_alloc &) {
		emMap->destroy();
		imMap->destroy();
		icMap->destroy();
		ancMap->destroy();

		releaseConstituentVec(cVec);
		releaseIonicFormVec(ifVec);

		cVec->destroy();
		ifVec->destroy();

		memset(&chemSystem, 0, sizeof(ChemicalSystem));
		memset(&calcProps, 0, sizeof(CalculatedProperties));

		return RetCode::E_NO_MEMORY;
	}

	tRet = initializeCalculatedProperties(calcProps, chemSystem);
	if (tRet != RetCode::OK) {
		emMap->destroy();
		imMap->destroy();
		icMap->destroy();
		ancMap->destroy();

		releaseConstituentVec(cVec);
		releaseIonicFormVec(ifVec);

		ancMap->destroy();
		cVec->destroy();
		ifVec->destroy();

		memset(&chemSystem, 0, sizeof(ChemicalSystem));
		memset(&calcProps, 0, sizeof(CalculatedProperties));

		return tRet;
	}

	return RetCode::OK;
}

void ECHMET_CC releaseChemicalSystem(ChemicalSystem &chemSystem) noexcept
{
	if (chemSystem.ionicForms != nullptr) {
		releaseIonicFormVec(chemSystem.ionicForms);
		chemSystem.ionicForms->destroy();
	}

	if (chemSystem.constituents != nullptr) {
		releaseConstituentVec(chemSystem.constituents);
		chemSystem.constituents->destroy();
	}

	chemSystem.analyticalConcentrationsByName->destroy();
	chemSystem.effectiveMobilitiesByName->destroy();
	chemSystem.ionicConcentrationsByName->destroy();
	chemSystem.ionicMobilitiesByName->destroy();
}

void ECHMET_CC releaseCalculatedProperties(CalculatedProperties &calcProps) noexcept
{
	calcProps.effectiveMobilities->destroy();
	calcProps.ionicMobilities->destroy();
	calcProps.ionicConcentrations->destroy();
}

void ECHMET_CC releaseInConstituent(const InConstituent &inC) noexcept
{
	inC.name->destroy();
	inC.pKas->destroy();
	inC.mobilities->destroy();

	InCFVec *cForms = inC.complexForms;
	if (cForms == nullptr)
		return;

	for (size_t jdx = 0; jdx < cForms->size(); jdx++) {
		const InLGVec *lgGroups = cForms->at(jdx).ligandGroups;
		if (lgGroups == nullptr)
			return;

		for (size_t kdx = 0; kdx < lgGroups->size(); kdx++) {
			const InLigandGroup &lgGrp = lgGroups->at(kdx);

			for (size_t ldx = 0; ldx < lgGrp.ligands->size(); ldx++) {
				const InLigandForm &iLF = lgGrp.ligands->at(ldx);

				iLF.ligandName->destroy();
				iLF.pBs->destroy();
				iLF.mobilities->destroy();
			}

			lgGrp.ligands->destroy();
		}

		lgGroups->destroy();
	}

	cForms->destroy();
}

void ECHMET_CC releaseInputData(const InConstituentVec *inVec) noexcept
{
	for (size_t idx = 0; idx < inVec->size(); idx++) {
		const InConstituent &inC = inVec->at(idx);

		releaseInConstituent(inC);
	}

	inVec->destroy();
}

} // namespace SysComp

} // name ECHMET
