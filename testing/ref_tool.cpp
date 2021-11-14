#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include "echmetelems.h"
#include "echmetionprops.h"
#include "jsonloader/inputreader.h"
#include "json_input_processor.h"
#include <echmetcaes.h>
#include <echmetcaes_extended.h>

void applyConcentrations(ECHMET::RealVec *acVec, const ECHMET::JsonInputProcessor::ConcentrationMap &acMap, const ECHMET::SysComp::ChemicalSystem &chemSystem)
{
	for (size_t idx = 0; idx < chemSystem.constituents->size(); idx++) {
		const ECHMET::SysComp::Constituent *ctuent = chemSystem.constituents->at(idx);

		(*acVec)[ctuent->analyticalConcentrationIndex] = acMap.at(ctuent->name->c_str());
	}
}

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

void printEquilibrium(const ECHMET::SysComp::ChemicalSystem &chemSystem, ECHMET::SysComp::CalculatedProperties &calcProps, ECHMET::NonidealityCorrections corrs)
{
	const ECHMET::IonProps::ComputationContext *ctx = ECHMET::IonProps::makeComputationContext(chemSystem, ECHMET::IonProps::ComputationContext::NONE);
	std::cout << "Ionic strength (mM): " << calcProps.ionicStrength * 1000 << "\n";
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
		const constituent_array_t ctarray = reader.read(inputDataFile);
		inputDesc = inputProc.process(&ctarray);
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

	printEquilibrium(chemSystem, calcProps, corrs);

	tRet = ECHMET::CAES::calculateBufferCapacity(bufferCap, corrs, chemSystem, calcProps, acVec);
	if (tRet != ECHMET::OK) {
		std::cerr << "Cannot calculate buffer capacity " << ECHMET::errorToString(tRet) << "\n";
		goto out_3;
	}

	std::cout << "Buffer capacity: " << bufferCap << "\n";

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

	return EXIT_SUCCESS;
}

int main(int argc, char **argv)
{
	return launch(argc, argv);
}

