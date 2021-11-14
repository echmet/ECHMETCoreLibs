#include <stdio.h>
#include "json_ldr_util.h"
#include "expected_json_ldr.h"

static int assign_object(json_t **item, json_t *root, const char *key)
{
	*item = json_object_get(root, key);
	if (*item == NULL) {
		fprintf(stderr, "Could not find key \"%s\" in JSON object\n", key);
		return 0;
	}
	return 1;
}

static size_t parse_expected_array(const json_t *array, expectation_t *expectations, size_t max, enum LoaderErrorCode *err)
{
	size_t idx;
	json_t *item;
	size_t ctr = 0;

	json_array_foreach(array, idx, item) {
		/* JSON objects */
		json_t *jDebyeHuckel;
		json_t *jOnsagerFuoss;
		json_t *jViscosity;
		json_t *jBufferCapacity;
		json_t *jIonicStrength;
		/* Values */
		int debyeHuckel;
		int onsagerFuoss;
		int viscosity;
		double bufferCapacity;
		double ionicStrength;

		if (idx > max)
			return ctr;

		if (!assign_object(&jDebyeHuckel, item, "debyeHuckel")) {
			*err = JLDR_E_BAD_INPUT;
			return ctr;
		}
		if (!json_is_boolean(jDebyeHuckel)) {
			fprintf(stderr, "Expected boolean for debyeHuckel\n");
			*err = JLDR_E_BAD_INPUT;
			return ctr;
		}
		debyeHuckel = json_boolean_value(jDebyeHuckel);

		if (!assign_object(&jOnsagerFuoss, item, "onsagerFuoss")) {
			*err = JLDR_E_BAD_INPUT;
			return ctr;
		}
		if (!json_is_boolean(jOnsagerFuoss)) {
			fprintf(stderr, "Expected boolean for onsagerFuoss\n");
			*err = JLDR_E_BAD_INPUT;
			return ctr;
		}
		onsagerFuoss = json_boolean_value(jOnsagerFuoss);

		if (!assign_object(&jViscosity, item, "viscosity")) {
			*err = JLDR_E_BAD_INPUT;
			return ctr;
		}
		if (!json_is_boolean(jViscosity)) {
			fprintf(stderr, "Expected boolean for viscosity\n");
			*err = JLDR_E_BAD_INPUT;
			return ctr;
		}
		viscosity = json_boolean_value(jViscosity);

		if (!assign_object(&jBufferCapacity, item, "bufferCapacity")) {
			*err = JLDR_E_BAD_INPUT;
			return ctr;
		}
		if (!json_is_real(jBufferCapacity)) {
			fprintf(stderr, "Expected real for bufferCapacity\n");
			*err = JLDR_E_BAD_INPUT;
			return ctr;
		}
		bufferCapacity = json_real_value(jBufferCapacity);
		if (bufferCapacity <= 0.0) {
			fprintf(stderr, "Expected buffer capacity is not positive. This is probably a mistake.\n");
			*err = JLDR_E_BAD_INPUT;
			return ctr;
		}

		if (!assign_object(&jIonicStrength, item, "ionicStrength")) {
			*err = JLDR_E_BAD_INPUT;
			return ctr;
		}
		if (!json_is_real(jIonicStrength)) {
			fprintf(stderr, "Expected real for ionicStrength\n");
			*err = JLDR_E_BAD_INPUT;
			return ctr;
		}
		ionicStrength = json_real_value(jIonicStrength);
		if (ionicStrength <= 0.0) {
			fprintf(stderr, "Expected ionic strength is not positive. This is probably a mistake.\n");
			*err = JLDR_E_BAD_INPUT;
			return ctr;
		}

		expectations[idx].debyeHuckel = debyeHuckel;
		expectations[idx].onsagerFuoss = onsagerFuoss;
		expectations[idx].viscosity = viscosity;
		expectations[idx].bufferCapacity = bufferCapacity;
		expectations[idx].ionicStrength = ionicStrength;

		ctr++;
	}

	return ctr;
}

static expectation_array_t parse_expected(const json_t *node, enum LoaderErrorCode *err)
{
	json_t *jExArray;
	size_t count;
	size_t realCount;
	expectation_t *expectations;
	expectation_array_t exArr = { NULL, 0 };

	jExArray = json_object_get(node, "expected");
	if (jExArray == NULL) {
		printf("No expected values\n");
		*err = JLDR_E_NOT_FOUND;
		return exArr;
	}

	count = json_array_size(jExArray);
	if (count == 0) {
		fprintf(stderr, "Expected values array is empty. This is not allowed\n");
		*err = JLDR_E_MALFORMED;
		return exArr;
	}

	expectations = (expectation_t *)calloc(count, sizeof(expectation_t));
	if (expectations == NULL) {
		fprintf(stderr, "Insufficient memory to store all expectations\n");
		*err = JLDR_E_NO_MEM;
		return exArr;
	}

	realCount = parse_expected_array(jExArray, expectations, count, err);

	exArr.expectations = expectations;
	exArr.count = realCount;

	return exArr;
}

void ldr_destroy_expectation(expectation_t *expectation)
{
	// Currently noop
}

void ldr_destroy_expectation_array(expectation_array_t *expectations)
{
	for (size_t idx = 0; idx < expectations->count; idx++)
	{
		ldr_destroy_expectation(&expectations->expectations[idx]);
	}
}

expectation_array_t ldr_load_expectations(const char *fileName, enum LoaderErrorCode *err)
{
	FILE *f;
	json_t *root;
	json_error_t jsonError;
	expectation_array_t exArr = { NULL, 0 };

	root = json_load_file(fileName, JSON_REJECT_DUPLICATES, &jsonError);
	if (root == NULL) {
		print_json_error(&jsonError);
		json_decref(root);
		*err = JLDR_E_MALFORMED;
		return exArr;
	}

	exArr = parse_expected(root, err);

	return exArr;
}
