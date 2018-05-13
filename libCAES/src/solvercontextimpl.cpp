#include "solvercontextimpl.hpp"

namespace ECHMET {
namespace CAES {

std::vector<int> makeChargesSquared(const int chMax)
{
	std::vector<int> chSq;

	chSq.resize(chMax + 1); /* Include zero charge */

	for (int ch = 0; ch <= chMax; ch++)
		chSq[ch] = ch * ch;

	return chSq;
}

} // namespace CAES
} // namespace ECHMET
