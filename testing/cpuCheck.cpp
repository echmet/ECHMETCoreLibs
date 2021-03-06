#include <echmetelems.h>
#include <algorithm>
#include <cctype>
#include <iostream>

static
std::string isSIMDAvailable(const std::string &name, bool ECHMET::CPUSIMD::* dp, ECHMET::CPUSIMD &simd)
{
	if (simd.*dp)
		return name + ": Available";
	return name + ": Not available";
}

static
std::string trimBoth(const char *cstr)
{
	std::string s{cstr};
	std::string tmp{};

	auto toTrim = [&](auto it, const auto stop) {
		for (; it != stop; it++) {
			if (!std::isspace(*it))
				break;
		}

		return it;
	};

	auto it = toTrim(s.cbegin(), s.cend());

	std::transform(it, s.cend(), std::back_inserter(tmp), [](const auto &c) { return c; });

	auto len = tmp.length();
	auto rit = toTrim(s.crbegin(), s.crend());

	len -= std::distance(s.crbegin(), rit);

	return tmp.substr(0, len);
}

int main()
{
	auto cpuIdent = ECHMET::cpuIdentifier();
	auto simd = ECHMET::cpuSupportedSIMD();

	if (cpuIdent != nullptr)
		std::cout << "CPU Identifier: " << trimBoth(cpuIdent->c_str()) << std::endl;

	std::cout << isSIMDAvailable("SSE2", &ECHMET::CPUSIMD::SSE2, simd) << std::endl;
	std::cout << isSIMDAvailable("SSE3", &ECHMET::CPUSIMD::SSE3, simd) << std::endl;
	std::cout << isSIMDAvailable("SSSE3", &ECHMET::CPUSIMD::SSSE3, simd) << std::endl;
	std::cout << isSIMDAvailable("SSE41", &ECHMET::CPUSIMD::SSE41, simd) << std::endl;
	std::cout << isSIMDAvailable("SSE42", &ECHMET::CPUSIMD::SSE42, simd) << std::endl;
	std::cout << isSIMDAvailable("AVX", &ECHMET::CPUSIMD::AVX, simd) << std::endl;
	std::cout << isSIMDAvailable("AVX2", &ECHMET::CPUSIMD::AVX2, simd) << std::endl;
	std::cout << isSIMDAvailable("AVX512", &ECHMET::CPUSIMD::AVX512, simd) << std::endl;
	std::cout << isSIMDAvailable("FMA3", &ECHMET::CPUSIMD::FMA3, simd) << std::endl;

	return 0;
}
	
