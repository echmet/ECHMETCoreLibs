#ifndef ECHMET_CAES_TYPES_HPP
#define ECHMET_CAES_TYPES_HPP

namespace ECHMET {
namespace CAES {

/*!
 * LigandIonicForm c-tor.
 *
 * @param[in] name Name of the ionic form.
 * @param[in] charge Electric charge of the ionic form.
 * @param[in] index Index of the ionic form (used in Jacobian generation).
 * @param[in] base Pointer to the ligand this form is based on.
 */
template <typename CAESReal>
LigandIonicForm<CAESReal>::LigandIonicForm(const std::string &name, const int charge, const size_t index, const Ligand<CAESReal> *base) :
	name(name),
	base(base),
	index(index),
	charge(charge),
	max(0)
{
}

/*!
 * LigandIonicForm c-tor.
 *
 * @param[in] other \param LigandIonicForm used as base of the "specialized" ligand ionic form.
 * @param[in] max Maximum number of this ligand ionic forms that can be bound to a given complex nucleus.
 * @param[in] pBs Vector of consecutive stability constants.
 */
template <typename CAESReal>
LigandIonicForm<CAESReal>::LigandIonicForm(const LigandIonicForm &other, const uint32_t max, const std::vector<CAESReal> &pBs) :
	name(other.name),
	base(other.base),
	index(other.index),
	charge(other.charge),
	max(max)
{
	if (pBs.size() != max)
		throw InvalidComplexationConfigurationException();

	const_cast<std::vector<CAESReal>&>(this->pBs) = pBs;
}

/*!
 * \p LigandIonicForm equality operator.
 *
 * @retval true The forms are considered equal.
 * @retval false The forms are considered not equal.
 */
template <typename CAESReal>
bool LigandIonicForm<CAESReal>::operator==(const LigandIonicForm &other) const
{
	return (this->name.compare(other.name) == 0);
}

/*!
 * Ligand c-tor.
 *
 * @param[in] name Name of the ligand.
 * @param[in] chargeLow Lowest electric charge the ligand can attain.
 * @param[in] chargeHigh Highest electric charge the ligand can attain.
 * @param[in] totalConcentrationIdx Index in the total concentrations vector for this ligand.
 * @param[in] pKas Vector of acidobazic dissociation constants for this ligand.
 */
template <typename CAESReal>
Ligand<CAESReal>::Ligand(const std::string &name, const int chargeLow, const int chargeHigh,
			 const size_t analyticalConcentrationIndex, const std::vector<CAESReal> &pKas) :
	name(name),
	chargeLow(chargeLow),
	chargeHigh(chargeHigh),
	analyticalConcentrationIndex(analyticalConcentrationIndex),
	pKas(pKas)
{
	ionicForms.reserve(chargeHigh - chargeLow + 1);
}

template <typename CAESReal>
using LigandVec =  std::vector<Ligand<CAESReal> *>;

/*!
 * ComplexNucleus c-tor.
 *
 * @param[in] name Name of the complex nucleus.
 * @param[in] chargeLow Lowest electric charge the complex nucleus can attain.
 * @param[in] chargeHigh Highest electric charge the complex nucleus can attain.
 * @param[in] totalConcentrationIdx Index in the total concentrations vector for this complex nuclues.
 * @param[in] exclusive Exclusive or non-exclusive complexation mode.
 * @param[in] ligandIFVec Vector of ligand ionic forms this complex nucleus complexes with.
 * @param[in] pKas Vector of acidobazic dissociation constants for this complex nucleus.
 */
template <typename CAESReal>
ComplexNucleus<CAESReal>::ComplexNucleus(const std::string &name, const int chargeLow, const int chargeHigh,
					 const size_t analyticalConcentrationIndex,
					 const std::vector<CAESReal> &pKas) :
	name(name),
	chargeLow(chargeLow),
	chargeHigh(chargeHigh),
	analyticalConcentrationIndex(analyticalConcentrationIndex),
	pKas(pKas)
{
}

/*!
 * ComplexNucleus d-tor.
 */
template <typename CAESReal>
ComplexNucleus<CAESReal>::~ComplexNucleus()
{
}

/*!
 * ContainedLigandIonicForm default c-tor.
 */
template <typename CAESReal>
ContainedLigandIonicForm<CAESReal>::ContainedLigandIonicForm() :
	count(0),
	lIF(nullptr)
{
}

/*!
 * ContainedLigandIonicForm c-tor.
 *
 * @param[in] count Count of the bound ionic forms.
 * @param[in] lIF The contained \p LigandIonicForm.
 */
template <typename CAESReal>
ContainedLigandIonicForm<CAESReal>::ContainedLigandIonicForm(const uint32_t count, const LigandIonicForm<CAESReal> *lIF) :
	count(count),
	lIF(lIF)
{
}

/*!
 * ContainedLigandIonicForm assignment operator.
 */
template <typename CAESReal>
ContainedLigandIonicForm<CAESReal> & ContainedLigandIonicForm<CAESReal>::operator=(const ContainedLigandIonicForm &other)
{
	const_cast<uint32_t&>(count) = other.count;
	const_cast<const LigandIonicForm<CAESReal> *&>(lIF) = other.lIF;

	return *this;
}

/*!
 * From c-tor.
 *
 * @param[in] name Name of the complex from.
 * @param[in] totalCharge Overall electric charge of the form.
 * @param[in] ligandIFIdx Index of the ligand ionic for used to make this complex form
 *            in the vector of all ligand ionic forms.
 * @param[in] ligandsContained Vector of all ligand ionic forms contained in the complex form.
 * @param[in] ancestorIdx Index of the complex form this complex form is based on in the vector
 *            of all complex forms.
 * @param[in] myIdx Index of this complex form in the vector of all complex forms.
 */
template <typename CAESReal>
Form<CAESReal>::Form(const std::string &name, const int totalCharge,
		     const size_t ligandIFIdx, const ContainedLigandIonicFormVec<CAESReal> &ligandsContained,
		     const size_t ancestorIdx, const size_t ancestorGlobalIdx,
		     const size_t myIdx, const CAESReal pB) :
	name(name),
	totalCharge(totalCharge),
	ligandIFIdx(ligandIFIdx),
	myIdx(myIdx),
	ancestorIdx(ancestorIdx),
	ancestorGlobalIdx(ancestorGlobalIdx),
	ligandsContained(ligandsContained),
	pB(pB)
{
}

} //namespace CAES
} //namespace ECHMET

#endif // ECHMET_CAES_TYPES_HPP
