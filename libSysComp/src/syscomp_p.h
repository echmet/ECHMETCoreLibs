#ifndef ECHMET_SYSCOMP_SYSCOMP_P_H
#define ECHMET_SYSCOMP_SYSCOMP_P_H

#include <echmetsyscomp.h>

#define ECHMET_IMPORT_INTERNAL
#include <containers/echmetvec_p.h>
#undef ECHMET_IMPORT_INTERNAL

#define _STRINGIFY(t) #t

namespace ECHMET {
namespace SysComp {

static RetCode buildComplexForms(Constituent *c, IonicFormVec *ifVec, const InConstituent &ic, const ConstituentVec *ligandsVec) noexcept;
static RetCode buildConstituentVec(ConstituentVec *cVec, IonicFormVec *ifVec, const InConstituentVec *inputData) noexcept;
static RetCode buildLigandIonicForms(const Constituent *c, IonicFormVec *ligandIFVec, const FixedString *name) noexcept;
static void calculateMaximumVariants(const size_t i, const InLFVec *ligandIFs, int32_t &total, int32_t accum) noexcept;
static Constituent * findLigand(const InLigandForm *lF, const ConstituentVec *cVec) noexcept;
static RetCode generateComplexForms(IonicForm *baseIF, IonicFormVec *groupIFVec, IonicFormVec *ifVec, InLFVec *ligandIFs, size_t formShift, int32_t toGenerate,
				    const ConstituentVec *ligandsVec) noexcept(false);
static void initializeChemicalSystemMapping(ChemicalSystem &chemSystem) noexcept(false);
static bool isConstituentContained(const ConstituentVec *cVec, const ConstituentVec *ligandsVec, const FixedString *newName) noexcept;
static bool isIonicFormContained(const IonicFormVec *ionicForms, const IonicForm *iF) noexcept;
static bool isLigandFormContained(const ContainedLigandIonicFormVec *containedLigandIFs, const InLigandForm *ligand, const int charge) noexcept;
static ContainedLigandIonicForm makeContainedLigandIonicForm(const Constituent *ligand, const int charge);
static RetCode makeIonicForm(IonicForm **iF, const Constituent *c, const int charge, const ECHMETReal &limitMobility, const FixedString *name) noexcept;
static void makeNonComplex(IonicForm *iF) noexcept;

template <typename T>
static
void releaseConcentrationVec(Vec<T> *vec) noexcept
{
	static_assert(std::is_pointer<T>::value, "Type " _STRINGIFY(T) " is not a pointer");

	for (size_t idx = 0; idx < vec->size(); idx++)
		delete vec->at(idx);
}

static void releaseConstituentVec(const ConstituentVec *vec) noexcept;
static void releaseIonicFormVec(const IonicFormVec *vec) noexcept;

} // namespace SysComp
} // ECHMET


#endif // ECHMET_SYSCOMP_SYSCOMP_P_H
