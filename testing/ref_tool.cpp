#include <cstring>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include "echmetelems.h"
#include "echmetionprops.h"
#include "echmetsyscomp.h"
#include "jsonloader/inputreader.h"
#include "jsonloader/expected_json_ldr.h"
#include "json_input_processor.h"
#include <echmetcaes.h>
#include <echmetcaes_extended.h>

#include <iomanip>

template <typename T>
bool tolCmp(T a, T b, T tol) {
	if (b == 0)
		return std::abs(a) <= tol;
	else
		return std::abs((a / b) - T(1)) <= tol;
}

static
void applyConcentrations(ECHMET::RealVec *acVec, const ECHMET::JsonInputProcessor::ConcentrationMap &acMap, const ECHMET::SysComp::ChemicalSystem &chemSystem)
{
	for (size_t idx = 0; idx < chemSystem.constituents->size(); idx++) {
		const ECHMET::SysComp::Constituent *ctuent = chemSystem.constituents->at(idx);

		(*acVec)[ctuent->analyticalConcentrationIndex] = acMap.at(ctuent->name->c_str());
	}
}

static
bool checkCorrectness(const char *inputDataFile, const ECHMET::SysComp::ChemicalSystem &chemSystem, const ECHMET::SysComp::CalculatedProperties &calcProps, const double bufferCapacity, const bool correctForDH, const bool correctForOF, const bool correctVS)
{
	const double TOL = 1.0e-12;
	bool failed = false;
	bool hasData = false;

	LoaderErrorCode err;
	expectation_array_t *expectations = ldr_load_expectations(inputDataFile, &err);
	if (!expectations) {
		if (err == JLDR_E_NOT_FOUND) {
			std::cout << "Test system contains to expected data to check against\n";
			return true;
		}
		std::cout << "Expected data in the test system are invalid\n";
		return false;
	}

	for (size_t idx = 0; idx < expectations->count; idx++) {
		const expectation_t *ex = &expectations->expectations[idx];

		if (ex->debyeHuckel != correctForDH || ex->onsagerFuoss != correctForOF || ex->viscosity != correctVS)
			continue;

		hasData = true;

		if (!tolCmp(bufferCapacity, ex->bufferCapacity, TOL)) {
			std::cout << "Unexpected buffer capacity, got " << bufferCapacity << ", expected " << ex->bufferCapacity << "\n";
			failed = true;
		}
		if (!tolCmp(calcProps.ionicStrength, ex->ionicStrength, TOL)) {
			std::cout << "Unexpected ionic strength, got " << calcProps.ionicStrength << ", expected " << ex->ionicStrength << "\n";
			failed = true;
		}

		std::cout << ex->ions.count << "\n";

		for (size_t jdx = 0; jdx < chemSystem.ionicForms->size(); jdx++) {
			const ECHMET::SysComp::IonicForm *iF = chemSystem.ionicForms->at(jdx);

			std::string iFName = "";
			if (iF->ifType == ECHMET::SysComp::H)
				iFName = "H+";
			else if (iF->ifType == ECHMET::SysComp::OH)
				iFName = "OH-";
			else
				iFName = iF->name->c_str();

			size_t kdx = 0;
			for (; kdx < ex->ions.count; kdx++) {
				const expected_ion_t *ion = &ex->ions.ions[kdx];

				if (!std::strcmp(ion->name, iFName.c_str())) {
					const double iFConc = calcProps.ionicConcentrations->at(iF->ionicConcentrationIndex);
					if (!tolCmp(iFConc, ion->concentration, TOL)) {
						std::cout << "Unexpected concentration of ionic form " << iFName << ", got " << iFConc << ", expected " << ion->concentration << "\n";
						failed = true;
					}

					const double iFMob = calcProps.ionicMobilities->at(iF->ionicMobilityIndex);
					if (!tolCmp(iFMob, ion->mobility, TOL)) {
						std::cout << "Unexpected mobility of ionic form " << iFName << ", got " << iFConc << ", expected " << ion->concentration << "\n";
						failed = true;
					}
					break;
				}
			}
			if (kdx == ex->ions.count) {
				std::cout << "Expected concentrations do not contain data for ionic form " << iFName << "\n";
				failed = true;
			}
		}
	}

	ldr_destroy_expectation_array(expectations);

	if (hasData) {
		if (failed)
			std::cout << "*** TEST FAILED !!! ***\n";
		else
			std::cout << "* Test OK *\n";
		return !failed;
	} else {
		std::cout << "*** Test system has no expected results for this combination of nonideality corrections\n";
		return true;
	}
}

static
ECHMET::NonidealityCorrections makeCorrections(const bool correctForDH, const bool correctForOF, const bool correctForVS)
{
	ECHMET::NonidealityCorrections corrs;

	if (correctForDH)
		ECHMET::nonidealityCorrectionSet(corrs, ECHMET::CORR_DEBYE_HUCKEL);
	if (correctForOF)
		ECHMET::nonidealityCorrectionSet(corrs, ECHMET::CORR_ONSAGER_FUOSS);
	if (correctForVS)
		ECHMET::nonidealityCorrectionSet(corrs, ECHMET::CORR_VISCOSITY);

	return corrs;
}

static
void printEquilibrium(const ECHMET::SysComp::ChemicalSystem &chemSystem, ECHMET::SysComp::CalculatedProperties &calcProps, ECHMET::NonidealityCorrections corrs, ECHMET::RealVec *acVec)
{
	ECHMET::IonProps::ComputationContext *ctx = ECHMET::IonProps::makeComputationContext(chemSystem, ECHMET::IonProps::ComputationContext::NONE);
	ECHMET::IonProps::correctMobilities(ctx, corrs, acVec, calcProps);
	std::cout << std::setprecision(13) << "Ionic strength (mM): " << calcProps.ionicStrength * 1000 << "\n";
	std::cout << "pH: " << ECHMET::IonProps::calculatepH(ctx, corrs, calcProps) << "\n";
	std::cout << "Ionic composition:\n";
	ctx->destroy();

	for (size_t idx = 0; idx < chemSystem.ionicForms->size(); idx++) {
		const ECHMET::SysComp::IonicForm *iF = chemSystem.ionicForms->at(idx);

		std::cout << "\t";

		switch (iF->ifType) {
		case ECHMET::SysComp::H:
			std::cout << "[H+]";
			break;
		case ECHMET::SysComp::OH:
			std::cout << "[OH-]";
			break;
		case ECHMET::SysComp::CONSTITUENT:
			std::cout << "[" << iF->name->c_str() << "]";
			break;
		}

		std::cout << " (mM): " << calcProps.ionicConcentrations->at(iF->ionicConcentrationIndex) << "\n";
	}

	std::cout << "Ionic mobilities:\n";
	for (size_t idx = 0; idx < chemSystem.ionicForms->size(); idx++) {
		const ECHMET::SysComp::IonicForm *iF = chemSystem.ionicForms->at(idx);

		std::cout << "\t";

		switch (iF->ifType) {
		case ECHMET::SysComp::H:
			std::cout << "[H+]";
			break;
		case ECHMET::SysComp::OH:
			std::cout << "[OH-]";
			break;
		case ECHMET::SysComp::CONSTITUENT:
			std::cout << "[" << iF->name->c_str() << "]";
			break;
		}

		std::cout << calcProps.ionicMobilities->at(iF->ionicMobilityIndex) << "\n";
	}
}

int launch(int argc, char **argv)
{
	const char *inputDataFile;
	bool correctForDH;
	bool correctForOF;
	bool correctForVS;
	ECHMET::JsonInputProcessor::InputDescription inputDesc;
	ECHMET::RetCode tRet;

	if (argc < 5) {
		std::cout << "Usage: inputFile DH_CORRECTION(number), OF_CORRECTION(Number), VS_CORRECTION(Number)\n";
		return EXIT_FAILURE;
	}

	inputDataFile = argv[1];

	try {
		correctForDH = std::atoi(argv[2]) >= 1;
		correctForOF = std::atoi(argv[3]) >= 1;
		correctForVS = std::atoi(argv[4]) >= 1;
	} catch (...) {
		std::cerr << "ERROR: Invalid input parameters\n";

		return EXIT_FAILURE;
	}

	try {
		ECHMET::JsonInputProcessor inputProc;
		InputReader reader;
		const constituent_array_t *ctarray = reader.read(inputDataFile);
		inputDesc = inputProc.process(ctarray);
	} catch (InputReader::MalformedInputException &ex) {
		std::cerr << ex.what();
		return EXIT_FAILURE;
	}

	ECHMET::SysComp::ChemicalSystem chemSystem;
	ECHMET::SysComp::CalculatedProperties calcProps;
	ECHMET::RealVec *acVec;
	ECHMET::CAES::SolverContext *solverCtx;
	ECHMET::CAES::Solver *solver;
	double bufferCap = -1;
	ECHMET::NonidealityCorrections corrs = makeCorrections(correctForDH, correctForOF, correctForVS);
	bool isCorrect = false;

	tRet = ECHMET::SysComp::makeComposition(chemSystem, calcProps, inputDesc.BGEComposition);
	if (tRet != ECHMET::OK) {
		std::cerr << "Cannot make composition " <<  ECHMET::errorToString(tRet) << std::endl;
		goto out;
	}
	acVec = ECHMET::createRealVec(chemSystem.constituents->size());
	if (acVec == ECHMET_NULLPTR) {
		std::cerr << "Cannot make analytical concentrations vector..." << std::endl;
		goto out_1;
	}
	acVec->resize(chemSystem.constituents->size());

	applyConcentrations(acVec, inputDesc.BGEConcentrations, chemSystem);

	tRet = ECHMET::CAES::createSolverContext(solverCtx, chemSystem);
	if (tRet != ECHMET::OK) {
		std::cerr << "Cannot create solver context " << ECHMET::errorToString(tRet) << std::endl;
		goto out_1;
	}

	solver = createSolver(solverCtx, ECHMET::CAES::Solver::defaultOptions(), corrs);
	if (solver == ECHMET_NULLPTR) {
		std::cerr << "Cannot create solver...";
		goto out_2;
	}

	tRet = solver->estimateDistributionSafe(acVec, calcProps);
	if (tRet != ECHMET::OK) {
		std::cerr << "Cannot estimate distribution " << ECHMET::errorToString(tRet);
		goto out_3;
	}

	tRet = solver->solve(acVec, calcProps, 5000, ECHMET_NULLPTR);
	if (tRet != ECHMET::OK) {
		std::cerr << "Cannot solve system " << ECHMET::errorToString(tRet) << std::endl;
		goto out_3;
	}

	printEquilibrium(chemSystem, calcProps, corrs, acVec);

	tRet = ECHMET::CAES::calculateBufferCapacity(bufferCap, corrs, chemSystem, calcProps, acVec);
	if (tRet != ECHMET::OK) {
		std::cerr << "Cannot calculate buffer capacity " << ECHMET::errorToString(tRet) << "\n";
		goto out_3;
	}

	std::cout << "Buffer capacity: " << bufferCap << "\n";

	isCorrect = checkCorrectness(inputDataFile, chemSystem, calcProps, bufferCap, correctForDH, correctForOF, correctForVS);

out_3:
	solver->destroy();
out_2:
	solverCtx->destroy();
	acVec->destroy();

out_1:
	ECHMET::SysComp::releaseChemicalSystem(chemSystem);
	ECHMET::SysComp::releaseCalculatedProperties(calcProps);
out:
	ECHMET::SysComp::releaseInputData(inputDesc.BGEComposition);

	return isCorrect ? EXIT_SUCCESS : EXIT_FAILURE;
}

int main(int argc, char **argv)
{
	return launch(argc, argv);
}

