#ifndef EXPECTED_LDR_H
#define EXPECTED_LDR_H

#ifdef __cplusplus
extern "C" {
#endif

#include "json_ldr.h"
#include <stddef.h>

struct expectation {
	int debyeHuckel;
	int onsagerFuoss;
	int viscosity;
	double ionicStrength;
	double bufferCapacity;
};
typedef struct expectation expectation_t;

struct expectation_array {
	expectation_t *expectations;
	size_t count;
};
typedef struct expectation_array expectation_array_t;

#ifdef __cplusplus
}
#endif

void ldr_destroy_expectation(expectation_t *expectation);
void ldr_destroy_expectation_array(expectation_array_t *expectations);
expectation_array_t ldr_load_expectations(const char *fileName, enum LoaderErrorCode *err);

#endif // EXPECTED_LDR_H
