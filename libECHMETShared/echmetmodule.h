#include <echmetcorelibs_config.h>

/* Enforce calling convention */
#ifndef ECHMET_CC
	#if defined ECHMET_PLATFORM_WIN32
		#define ECHMET_CC __stdcall
	#elif defined ECHMET_PLATFORM_UNIX
		#ifdef ECHMET_COMPILER_GCC_LIKE
			#ifdef __i386__
				#define ECHMET_CC __attribute__((__cdecl__))
			#else
				#define ECHMET_CC
			#endif // __x86_64__
		#else
			#error "Unsupported or misdetected compiler"
		#endif // ECHMET_COMPILER_*
	#else
		#error "Unsupported or misdetected target platform"
	#endif // ECHMET_PLATFORM_*
#endif // ECHMET_CC

/* Define force-inlined functions */
#ifndef ECHMET_FORCE_INLINE
	#if defined ECHMET_PLATFORM_WIN32
		#if defined ECHMET_COMPILER_MINGW || defined ECHMET_COMPILER_MSYS
			#define ECHMET_FORCE_INLINE inline __attribute__((always_inline))
		#elif defined ECHMET_COMPILER_MSVC
			#define ECHMET_FORCE_INLINE __forceinline
		#else
			#error "Unsupported or misdetected compiler"
		#endif // ECHMET_COMPILER_*
	#elif defined ECHMET_PLATFORM_UNIX
		#ifdef ECHMET_COMPILER_GCC_LIKE
			#define ECHMET_FORCE_INLINE inline __attribute__((always_inline))
		#else
			#error "Unsupported or misdetected compiler"
		#endif // ECHMET_COMPILER_*
	#else
		#error "Unsupported or misdetected target platform"
	#endif // ECHMET_PLATFORM_*
#endif // ECHMET_FORCE_INLINE

/* Allow for redefinitions of ECHMET_API as needed */
#ifdef ECHMET_API
	#undef ECHMET_API
#endif // ECHMET_API

/* Export only symbols that are part of the public API */
#if defined ECHMET_PLATFORM_WIN32
	#if defined ECHMET_DLL_BUILD && !defined(ECHMET_IMPORT_INTERNAL)
		#if defined ECHMET_COMPILER_MINGW || defined ECHMET_COMPILER_MSYS
			#define ECHMET_API __attribute__ ((dllexport))
		#elif defined ECHMET_COMPILER_MSVC
			#define ECHMET_API __declspec(dllexport)
		#else
			#error "Unsupported or misdetected compiler"
		#endif // ECHMET_COMPILER_*
	#else
		#if defined ECHMET_COMPILER_MINGW || defined ECHMET_COMPILER_MSYS
			#define ECHMET_API __attribute__ ((dllimport))
		#elif defined ECHMET_COMPILER_MSVC
			#define ECHMET_API __declspec(dllimport)
		#else
			#error "Unsupported or misdetected compiler"
		#endif // ECHMET_COMPILER_*
	#endif // ECHMET_DLL_BUILD
#elif defined ECHMET_PLATFORM_UNIX
	#ifdef ECHMET_COMPILER_GCC_LIKE
		#if defined ECHMET_DLL_BUILD && !defined(ECHMET_IMPORT_INTERNAL)
			#define ECHMET_API __attribute__ ((visibility ("default")))
		#else
			#define ECHMET_API
		#endif // ECHMET_DLL_BUILD
	#else
		#error "Unsupported or misdetected compiler"
	#endif // ECHMET_COMPILER_*
#else
	#error "Unsupported or misdetected target platform"
#endif // ECHMET_PLATFORM_*

/*
 * EXPLANATION:
 * Okay, you have probably just finished examining this header file and now you
 * are asking yourself what the hell is this supposed to do and what was wrong
 * with the guy who wrote it this way? You can find the answer only to the
 * former here.
 *
 * All ECHMET modules are intended to work on various platforms and be buildable
 * with different compilers. This is, obviously, a pain in the ass to take care
 * of properly. This header attempts to take care of four issues in a reasonably
 * reusable way, namely:
 *
 * - Target platform detection
 * - Target architecture detection (x86 and x86_64 only at this point)
 * - Compiler detection
 * - Symbol import/export mode.
 *
 * We need to know the target platform and architecture in order to enforce
 * correct calling convention. This is done in the first block of ifdefs
 * and it is pretty easy.
 *
 * Compiler detection is closely related to symbol import/export mode because
 * different compilers understand different directives that govern this.
 * Any symbol tagged with ECHMET_API will be marked as global (=exported)
 * in the resulting binary as long as the following conditions are met:
 *
 *  - The binary is being build in DLL mode and ECHMET_DLL_BUILD is defined
 *  - Export is overridden by ECHMET_IMPORT_INTERNAL
 *
 * If these conditions are not met, the symbol is marked for import (on Windows)
 * or not marked at all (on UNIX).
 *
 * Notice that this file does not have include guards. This allows us to
 * switch the ECHMET_IMPORT_INTERNAL flag on and off as needed by
 * defining or undefining it and reincluding this header.
 * This - although a bit messy - is useful when we are building a library A
 * that exports symbols through the ECHMET_API tag but it also pulls in
 * other library B that uses the same mechanism to export symbols. If the symbols
 * exported by library B were not wrapped in an ECHMET_IMPORT_INTERNAL block,
 * library A would export symbols from both A and B which is not what we want.
 *
 * Clear? ...yeah, didn't think so.
 */

/* C++11 detection */
#ifndef ECHMET_CPP11_DETECTION_DONE
#define ECHMET_CPP11_DETECTION_DONE

#if __cplusplus >= 201103L
	#define ECHMET_HAVE_CPP_11
#endif

/* Not funny, MSVC guys!!! */
#ifndef ECHMET_HAVE_CPP_11
	#if _MSC_VER >= 1900
		#define ECHMET_HAVE_CPP_11
	#endif
#endif

#ifdef ECHMET_HAVE_CPP_11
	#include <type_traits>
	#define ECHMET_NOEXCEPT noexcept
	#define ECHMET_ST_ENUM(ename) enum class ename : int32_t
	#define ECHMET_WK_ENUM(ename) enum ename : int32_t
	#define ECHMET_NULLPTR nullptr
	#define ENUM_FORCE_INT32_SIZE(dummy)
	#define ECHMET_ST_ENUM_UTYPE(ename) std::underlying_type<ename>::type
#else /* C++11-unaware compiler */
	#define ECHMET_NOEXCEPT throw()
	#define ECHMET_ST_ENUM(ename) enum ename
	#define ECHMET_WK_ENUM(ename) enum ename
	#define ECHMET_NULLPTR 0x0
	#define ENUM_FORCE_INT32_SIZE(ename) , __ECHMET_ST_ENUM_PLACEHOLDER ## ename = (1 << 31)
	#define ECHMET_ST_ENUM_UTYPE(ename) int32_t
#endif // __cplusplus

#endif // ECHMET_CPP11_DETECTION_DONE
