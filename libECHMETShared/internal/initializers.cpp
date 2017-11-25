/*!
 * Initializers to set default runtime parameters
 */

#include <echmetelems.h>
#ifdef ECHMET_USE_HIGH_PRECISION
#include "mpreal.h"
#endif // ECHMET_USE_HIGH_PRECISION

#define MPFR_DEFAULT_DIGITS 100

#ifdef ECHMET_PLATFORM_UNIX
void __attribute((constructor)) initializeLibrary()
{
#ifdef ECHMET_USE_HIGH_PRECISION
	mpfr::mpreal::set_default_prec(mpfr::digits2bits(MPFR_DEFAULT_DIGITS));
#endif // ECHMET_USE_HIGH_PRECISION
}
#elif defined ECHMET_PLATFORM_WIN32
// dllmain.cpp : Defines the entry point for the DLL application.
#include <windows.h>

BOOL APIENTRY DllMain(HMODULE hModule, DWORD ul_reason_for_call, LPVOID lpReserved)
{
	(void)hModule;
	(void)lpReserved;

	switch (ul_reason_for_call) {
	case DLL_PROCESS_ATTACH:
#ifdef ECHMET_USE_HIGH_PRECISION
		mpfr::mpreal::set_default_prec(mpfr::digits2bits(MPFR_DEFAULT_DIGITS));
#endif // ECHMET_USE_HIGH_PRECISION
		break;
	case DLL_THREAD_ATTACH:
	case DLL_THREAD_DETACH:
	case DLL_PROCESS_DETACH:
		break;
	}

	return TRUE;
}

#endif // ECHMET_PLATFORM_
