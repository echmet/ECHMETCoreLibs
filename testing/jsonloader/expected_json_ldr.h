#ifndef EXPECTED_LDR_H
#define EXPECTED_LDR_H

#ifdef __cplusplus
extern "C" {
#endif

#include "json_ldr.h"
#include <stddef.h>

typedef struct expected_ion {
	char *name;
	double concentration;
	double mobility;
} expected_ion_t;

typedef struct expected_ion_array {
	expected_ion_t *ions;
	size_t count;
} expected_ion_array_t;

typedef struct expectation {
	int debyeHuckel;
	int onsagerFuoss;
	int viscosity;
	double ionicStrength;
	double bufferCapacity;
	expected_ion_array_t ions;
} expectation_t;

typedef struct expectation_array {
	expectation_t *expectations;
	size_t count;
} expectation_array_t;

void ldr_destroy_expectation(expectation_t *expectation);
void ldr_destroy_expectation_array(expectation_array_t *array);
expectation_array_t * ldr_load_expectations(const char *fileName, enum LoaderErrorCode *err);

#ifdef __cplusplus
}
#endif

#endif // EXPECTED_LDR_H
