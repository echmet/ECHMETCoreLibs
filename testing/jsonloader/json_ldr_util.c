#include "json_ldr_util.h"

void print_json_error(const json_error_t *error)
{
	fprintf(stderr, "Error has occured when parsing JSON structure:\n"
			"Error: %s\n"
			"Source: %s\n"
			"Line: %d\n"
			"Column: %d\n",
			error->text,
			error->source,
			error->line,
			error->column);
}
