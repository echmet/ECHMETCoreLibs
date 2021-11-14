#include <jansson.h>
#include <stdio.h>
#include <string.h>
#include "json_ldr_util.h"
#include "expected_json_ldr.h"

static
int assign_object(json_t **item, const json_t *root, const char *key)
{
	*item = json_object_get(root, key);
	if (*item == NULL) {
		fprintf(stderr, "Could not find key \"%s\" in JSON object\n", key);
		return 0;
	}
	return 1;
}

static
expected_ion_array_t parse_expected_ions(const json_t *array, enum LoaderErrorCode *err)
{
	size_t idx;
	json_t *value;
	size_t sz = json_array_size(array);

	expected_ion_array_t expected = {
		calloc(sz, sizeof(expected_ion_t)),
		0
	};

	memset(expected.ions, 0, sz * sizeof(expected_ion_t));

	json_array_foreach(array, idx, value) {
		json_t *jName;
		json_t *jConc;
		json_t *jMob;

		const char *auxName;
		char *name;

		/* Get ion name */
		if (!assign_object(&jName, value, "name")) {
			*err = JLDR_E_BAD_INPUT;
			return expected;
		}
		if (!json_is_string(jName)) {
			*err = JLDR_E_BAD_INPUT;
			return expected;
		}
		auxName = json_string_value(jName);
		name = (char *)malloc(strlen(auxName) + 1);
		strcpy(name, auxName);

		/* Get ion concentration */
		if (!assign_object(&jConc, value, "concentration")) {
			free(name);
			*err = JLDR_E_BAD_INPUT;
			return expected;
		}
		if (!json_is_real(jConc)) {
			free(name);
			*err = JLDR_E_BAD_INPUT;
			return expected;
		}

		/* Get ionic mobility */
		if (!assign_object(&jMob, value, "mobility")) {
			free(name);
			*err = JLDR_E_BAD_INPUT;
			return expected;
		}
		if (!json_is_real(jMob)) {
			free(name);
			*err = JLDR_E_BAD_INPUT;
			return expected;
		}

		expected.ions[idx].name = name;
		expected.ions[idx].concentration = json_real_value(jConc);
		expected.ions[idx].mobility = json_real_value(jMob);
		expected.count++;
	}

	*err = JLDR_OK;
	return expected;
}

static
expectation_array_t * parse_expected(const json_t *array, enum LoaderErrorCode *err)
{
	size_t idx;
	json_t *value;
	expectation_array_t *expectations = malloc(sizeof(expectation_array_t));
	expectations->expectations = calloc(json_array_size(array), sizeof(expectation_t));
	expectations->count = 0;

	json_array_foreach(array, idx, value) {
		/* JSON objects */
		json_t *jDebyeHuckel;
		json_t *jOnsagerFuoss;
		json_t *jViscosity;
		json_t *jBufferCapacity;
		json_t *jIonicStrength;
		json_t *jConcentrations;

		expectation_t *expectation = &expectations->expectations[idx];
		memset(expectation, 0, sizeof(expectation_t));

		if (!assign_object(&jDebyeHuckel, value, "debyeHuckel")) {
			*err = JLDR_E_BAD_INPUT;
			return expectations;
		}
		if (!json_is_boolean(jDebyeHuckel)) {
			fprintf(stderr, "Expected boolean for debyeHuckel\n");
			*err = JLDR_E_BAD_INPUT;
			return expectations;
		}
		expectation->debyeHuckel = json_boolean_value(jDebyeHuckel);

		if (!assign_object(&jOnsagerFuoss, value, "onsagerFuoss")) {
			*err = JLDR_E_BAD_INPUT;
			return expectations;
		}
		if (!json_is_boolean(jOnsagerFuoss)) {
			fprintf(stderr, "Expected boolean for onsagerFuoss\n");
			*err = JLDR_E_BAD_INPUT;
			return expectations;
		}
		expectation->onsagerFuoss = json_boolean_value(jOnsagerFuoss);

		if (!assign_object(&jViscosity, value, "viscosity")) {
			*err = JLDR_E_BAD_INPUT;
			return expectations;
		}
		if (!json_is_boolean(jViscosity)) {
			fprintf(stderr, "Expected boolean for viscosity\n");
			*err = JLDR_E_BAD_INPUT;
			return expectations;
		}
		expectation->viscosity = json_boolean_value(jViscosity);

		if (!assign_object(&jBufferCapacity, value, "bufferCapacity")) {
			*err = JLDR_E_BAD_INPUT;
			return expectations;
		}
		if (!json_is_real(jBufferCapacity)) {
			fprintf(stderr, "Expected real for bufferCapacity\n");
			*err = JLDR_E_BAD_INPUT;
			return expectations;
		}

		expectation->bufferCapacity = json_real_value(jBufferCapacity);
		if (expectation->bufferCapacity <= 0.0) {
			fprintf(stderr, "Expected buffer capacity is not positive. This is probably a mistake.\n");
			*err = JLDR_E_BAD_INPUT;
			return expectations;
		}

		if (!assign_object(&jIonicStrength, value, "ionicStrength")) {
			*err = JLDR_E_BAD_INPUT;
			return expectations;
		}
		if (!json_is_real(jIonicStrength)) {
			fprintf(stderr, "Expected real for ionicStrength\n");
			*err = JLDR_E_BAD_INPUT;
			return expectations;
		}
		expectation->ionicStrength = json_real_value(jIonicStrength);
		if (expectation->ionicStrength <= 0.0) {
			fprintf(stderr, "Expected ionic strength is not positive. This is probably a mistake.\n");
			*err = JLDR_E_BAD_INPUT;
			return expectations;
		}

		if (!assign_object(&jConcentrations, value, "ions")) {
			*err = JLDR_E_BAD_INPUT;
			return expectations;
		}
		if (!json_is_array(jConcentrations)) {
			fprintf(stderr, "Expected array for ions\n");
			*err = JLDR_E_BAD_INPUT;
			return expectations;
		}

		expectation->ions = parse_expected_ions(jConcentrations, err);
		if (*err != JLDR_OK) {
			fprintf(stderr, "Expected ions are invalid\n");
			return expectations;
		}

		expectations->count++;
	}

	*err = JLDR_OK;

	return expectations;
}

static
expectation_array_t * parse_expected_array(const json_t *node, enum LoaderErrorCode *err)
{
	json_t *jEx;
	expectation_array_t *ex;

	if (!assign_object(&jEx, node, "expected")) {
		*err = JLDR_E_NOT_FOUND;
		return NULL;
	}

	ex = parse_expected(jEx, err);
	if (*err != JLDR_OK) {
		ldr_destroy_expectation_array(ex);
		return NULL;
	}

	return ex;
}

void ldr_destroy_expectation(expectation_t *expectation)
{
	for (size_t idx = 0; idx < expectation->ions.count; idx++)
		free(expectation->ions.ions[idx].name);
	free(expectation->ions.ions);
}

void ldr_destroy_expectation_array(expectation_array_t *array)
{
	for (size_t idx = 0; idx < array->count; idx++)
		ldr_destroy_expectation(&array->expectations[idx]);
	free(array);
}

expectation_array_t * ldr_load_expectations(const char *fileName, enum LoaderErrorCode *err)
{
	FILE *f;
	json_t *root;
	json_error_t jsonError;
	expectation_array_t *ex;

	root = json_load_file(fileName, JSON_REJECT_DUPLICATES, &jsonError);
	if (root == NULL) {
		print_json_error(&jsonError);
		json_decref(root);
		*err = JLDR_E_MALFORMED;
		return NULL;
	}

	ex = parse_expected_array(root, err);
	json_decref(root);

	return ex;
}
