#include <errno.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <jansson.h>
#include "constituents_json_ldr.h"

static void print_json_error(const json_error_t *error)
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

static const char * type_to_string(const enum ConstituentType t)
{
	switch (t) {
	case LIGAND:
		return "Ligand";
	case NUCLEUS:
		return "Nucleus";
	default:
		return "INVALID";
	}
}

static const char * role_to_string(const enum ConstituentRole r)
{
	switch (r) {
	case BACKGROUND:
		return "Background";
	case ANALYTE:
		return "Analyte";
	default:
		return "INVALID";
	}
}

static void print_double_array(const double *array, const int len, const char *prefix, const char *fmt)
{
	int idx;

	printf("%s{ ", prefix);
	for (idx = 0; idx < len; idx++) {
		printf(fmt, array[idx]);
		printf("; ");
	}
	printf("}\n");
}

static void print_all_constituents(const constituent_t *allCts, const size_t len)
{
	size_t idx;

	for (idx = 0; idx < len; idx++) {
		const constituent_t *ctuent = &allCts[idx];

		printf("Name: %s\n"
		       "Type: %s\n"
		       "Role: %s\n"
		       "low charge = %d, high charge = %d\n",
		       ctuent->name,
		       type_to_string(ctuent->ctype),
		       role_to_string(ctuent->crole),
		       ctuent->chargeLow, ctuent->chargeHigh);
		printf("concentration in BGE = %f\n", ctuent->concentrationBGE);
		printf("concentration in Sample = %f\n", ctuent->concentrationSample);
		printf("pKa = ");
		printf("viscosity coefficient = %f\n", ctuent->viscosityCoefficient);
		print_double_array(ctuent->pKas, ctuent->chargeHigh - ctuent->chargeLow, "", "%lf");
		printf("mobilities = ");
		print_double_array(ctuent->mobilities, ctuent->chargeHigh - ctuent->chargeLow + 1, "", "%lg");

		if (ctuent->ctype == NUCLEUS) {
			size_t iidx;

			for (iidx = 0; iidx < ctuent->complexFormsCount; iidx++) {
				size_t jidx;
				const complex_form_t *cForm = &(ctuent->complexForms[iidx]);

				printf("\tFor charge %d: \n", cForm->nucleusCharge);
				if (cForm->ligandGroups == NULL) {
					printf("\t\tThis constituent forms no complexes with charge %d\n", cForm->nucleusCharge);
					continue;
				}

				for (jidx = 0; jidx < cForm->count; jidx++) {
					const ligand_group_t *lGroup = &(cForm->ligandGroups[jidx]);

					printf("\t\tLigand group %zu\n", jidx);

					size_t kidx;
					for (kidx = 0; kidx < lGroup->count; kidx++) {
						const ligand_form_t *lForm = &(lGroup->ligandForms[kidx]);

						printf("\t\t\tLigand name: %s\n"
						       "\t\t\tLigand charge: %d\n"
						       "\t\t\tMax ligands: %d\n",
						       lForm->name,
						       lForm->charge,
						       lForm->maxCount);
						print_double_array(lForm->pBs, lForm->maxCount, "\t\t", "%lf");
						print_double_array(lForm->mobilities, lForm->maxCount, "\t\t", "%lg");
					}
				}
			}
		}
	}
}

static void init_constituent(constituent_t *ctuent,
			     const char *ctype, const char *crole, const char *name,
			     const int chargeLow, const int chargeHigh,
			     const double concentrationBGE, const double concentrationSample,
			     double *const pKas, double *const mobilities,
			     const double viscosityCoefficient)
{
	ctuent->chargeLow = chargeLow;
	ctuent->chargeHigh = chargeHigh;
	ctuent->concentrationBGE = concentrationBGE;
	ctuent->concentrationSample = concentrationSample;
	ctuent->pKas = pKas;
	ctuent->mobilities = mobilities;
	ctuent->viscosityCoefficient = viscosityCoefficient;

	ctuent->name = calloc(strlen(name) + 1, sizeof(char));
	strcpy(ctuent->name, name);

	fprintf(stderr, "Constituent type: %s\n", ctype);
	if (strcmp(ctype, "L") == 0) {
		ctuent->ctype = LIGAND;
		fprintf(stderr, "--- LIGAND\n");
	} else if (strcmp(ctype, "N") == 0) {
		ctuent->ctype = NUCLEUS;
		fprintf(stderr, "--- NUCLEUS\n");
	} else {
		fprintf(stderr, "--- INVALID\n");
		ctuent->ctype = INVALID_TYPE;
	}

	if (strcmp(crole, "B") == 0)
		ctuent->crole = BACKGROUND;
	else if (strcmp(crole, "A") == 0)
		ctuent->crole = ANALYTE;
	else
		ctuent->crole = INVALID_ROLE;
}

static bool is_json_number(const json_t *item)
{
	return json_is_integer(item) || json_is_real(item);
}

static double * parse_array_reals(const json_t *jArray, const size_t length)
{
	size_t idx;
	double *array;

	if (length < 1)
		return NULL;

	array = calloc(length, sizeof(double));
	if (array == NULL)
		return NULL;

	for (idx = 0; idx < length; idx++) {
		const json_t *item = json_array_get(jArray, idx);

		if (item == NULL) {
			free(array);
			return NULL;
		}

		if (is_json_number(item) == 0) {
			fprintf(stderr, "Item is not a real number\n");
			free(array);
			return NULL;
		}

		array[idx] = json_number_value(item);
	}

	return array;
}

static int parse_ligands(const json_t *jLigands, ligand_group_t *lGroup)
{
	size_t idx;
	json_t *item;
	size_t ligandCount;

	lGroup->count = json_array_size(jLigands);
	if (lGroup->count < 1) {
		lGroup->ligandForms = NULL;
		fprintf(stderr, "No ligands in \"ligands\" array\n");
		return 0;
	}
	lGroup->ligandForms = calloc(lGroup->count, sizeof(ligand_form_t));

	json_array_foreach(jLigands, idx, item) {
		json_t *jElement;
		const char *name;
		int charge;
		size_t maxCount;
		size_t iidx;
		double *pBs;
		double *mobilities;
		ligand_form_t *lForm = &lGroup->ligandForms[idx];

		lForm->pBs = NULL;
		lForm->mobilities = NULL;

		if (json_is_object(item) == 0) {
			fprintf(stderr, "Item is not an object\n");
			free(lGroup->ligandForms);
			return 0;
		}

		jElement = json_object_get(item, "name");
		if (jElement == NULL) {
			fprintf(stderr, "No key \"name\" in the ligandForm element\n");
			free(lGroup->ligandForms);
			return 0;
		}
		if (json_is_string(jElement) == 0) {
			fprintf(stderr, "Item \"name\" is not a string\n");
			free(lGroup->ligandForms);
			return 0;
		}
		name = json_string_value(jElement);

		jElement = json_object_get(item, "charge");
		if (jElement == NULL) {
			fprintf(stderr, "No key \"charge\" in the ligandForm element\n");
			free(lGroup->ligandForms);
			return 0;
		}
		if (json_is_integer(jElement) == 0) {
			fprintf(stderr, "Item \"charge\" is not an integer\n");
			free(lGroup->ligandForms);
			return 0;
		}
		charge = json_integer_value(jElement);

		jElement = json_object_get(item, "maxCount");
		if (jElement == NULL) {
			fprintf(stderr, "No key \"maxCount\" in the ligandForm element\n");
			free(lGroup->ligandForms);
			return 0;
		}
		if (json_is_integer(jElement) == 0) {
			fprintf(stderr, "Item \"maxCount\" is not an integer\n");
			free(lGroup->ligandForms);
			return 0;
		}
		maxCount = json_integer_value(jElement);

		jElement = json_object_get(item, "pBs");
		if (jElement == NULL) {
			fprintf(stderr, "No key \"pBs\" in the ligandForm element\n");
			free(lGroup->ligandForms);
			return 0;
		}
		if (json_is_array(jElement) == 0) {
			fprintf(stderr, "Item \"pBs\" is not an array\n");
			free(lGroup->ligandForms);
			return 0;
		}
		if (json_array_size(jElement) < maxCount) {
			fprintf(stderr, "Sizes of pBs and maxCount do not match\n");
			free(lGroup->ligandForms);
			return 0;
		}

		pBs = parse_array_reals(jElement, maxCount);
		if (pBs == NULL) {
			fprintf(stderr, "Cannot parse pBs array\n");
			free(lGroup->ligandForms);
			return 0;
		}

		jElement = json_object_get(item, "mobilities");
		if (jElement == NULL) {
			fprintf(stderr, "No key \"mobilities\" in the ligandForm element\n");
			free(lGroup->ligandForms);
			free(pBs);
			return 0;
		}
		if (json_is_array(jElement) == 0) {
			fprintf(stderr, "Item \"mobilities\" is not an array\n");
			free(lGroup->ligandForms);
			free(pBs);
			return 0;
		}
		if (json_array_size(jElement) != maxCount) {
			fprintf(stderr, "Sizes of \"mobilities\" and \"jMaxCount\" do not match\n");
			free(lGroup->ligandForms);
			free(pBs);
			return 0;
		}

		mobilities = parse_array_reals(jElement, maxCount);
		if (mobilities == NULL) {
			fprintf(stderr, "Cannot parse \"mobilities\" array\n");
			free(lGroup->ligandForms);
			free(pBs);
			return 0;
		}

		lForm->name = calloc(strlen(name) + 1, sizeof(char));
		strcpy(lForm->name, name);
		lForm->charge = charge;
		lForm->maxCount = maxCount;
		lForm->pBs = pBs;
		lForm->mobilities = mobilities;
	}

	return 1;
}

static int parse_ligand_groups(const json_t *jLgArray, ligand_group_t **lGroups, size_t *count)
{
	size_t idx;
	json_t *item;

	*count = json_array_size(jLgArray);

	if (*count == 0) {
		fprintf(stderr, "No ligands in \"ligandGroup\"\n");
		return -1;
	}

	*lGroups = calloc(*count, sizeof(ligand_group_t));

	json_array_foreach(jLgArray, idx, item) {
		json_t *jElement;
		ligand_group_t *lGrp = &(*lGroups)[idx];

		jElement = json_object_get(item, "ligands");
		if (jElement == NULL) {
			fprintf(stderr, "No key \"ligands\" in the ligandGroup element\n");
			free(*lGroups);
			return 0;
		}
		if (json_is_array(jElement) == 0) {
			fprintf(stderr, "Item \"ligands\" is not an array\n");
			free(*lGroups);
			return 0;
		}
		lGrp->count = json_array_size(jElement);

		if (parse_ligands(jElement, lGrp) == 0) {
			free(*lGroups);
			return 0;
		}
	}

	return 1;
}

static int parse_complex_forms(const json_t *node, constituent_t *ctuent)
{
	size_t idx;
	json_t *item;
	json_t *jComplexForms;
	complex_form_t *cpxForms;
	int cpxCount;

	jComplexForms = json_object_get(node, "complexForms");
	if (jComplexForms == NULL) {
		fprintf(stderr, "No key \"complexForms\" in the constituent element\n");
		return -1;
	}
	if (json_is_array(jComplexForms) == 0) {
		fprintf(stderr, "Item \"complexForms\" is not an array\n");
		return -1;
	}
	cpxCount = json_array_size(jComplexForms);
	if (cpxCount < 1) {
		ctuent->complexForms = NULL;
		ctuent->complexFormsCount = 0;
		return 1;
	}

	cpxForms = calloc(cpxCount, sizeof(complex_form_t));

	json_array_foreach(jComplexForms, idx, item) {
		int nucleusCharge;
		json_t *jElement;
		complex_form_t *cForm = &cpxForms[idx];

		jElement = json_object_get(item, "nucleusCharge");
		if (jElement == NULL) {
			fprintf(stderr, "No key \"nucleusCharge\" in the complexForm element\n");

			free(cpxForms);
			return -1;
		}
		if (json_is_integer(jElement) == 0) {
			fprintf(stderr, "Item \"nucleusCharge\" is not an integer\n");

			free(cpxForms);
			return -1;
		}
		nucleusCharge = json_integer_value(jElement);
		cForm->nucleusCharge = nucleusCharge;

		jElement = json_object_get(item, "ligandGroups");
		if (jElement == NULL) {
			fprintf(stderr, "No key \"ligandGroups\" in the complexForm element\n");

			free(cpxForms);
			return -1;
		}
		if (json_is_array(jElement) == 0) {
			fprintf(stderr, "Item \"ligandGroups\" is not an array\n");

			free(cpxForms);
			return -1;
		}

		if (parse_ligand_groups(jElement, &cForm->ligandGroups, &cForm->count) == 0) {
			fprintf(stderr, "Cannot parse array \"ligandGroups\"\n");

			free(cpxForms);
			return -1;
		}
	}

	ctuent->complexForms = cpxForms;
	ctuent->complexFormsCount = idx;
	return 1;
}

static size_t parse_cts_array(const json_t *array, constituent_t *allCts, size_t max, enum LoaderErrorCode *err)
{
	size_t idx;
	json_t *item;
	size_t ctr = 0;

	json_array_foreach(array, idx, item) {
		const char *typeId;
		const char *roleId;
		const char *name;
		json_int_t chargeLow;
		json_int_t chargeHigh;
		json_t *jTypeId;
		json_t *jRoleId;
		json_t *jName;
		json_t *jChargeLow;
		json_t *jChargeHigh;
		json_t *jConcentrationBGE;
		json_t *jConcentrationSample;
		json_t *jViscosityCoefficient;
		json_t *jPKaArray;
		json_t *jMobilities;
		double concentrationBGE;
		double concentrationSample;
		double viscosityCoefficient;
		double *pKas;
		double *mobilities;
		constituent_t *ctuent = &allCts[ctr];

		if (idx >= max)
			return ctr;

		jTypeId = json_object_get(item, "type");
		if (jTypeId == NULL) {
			fprintf(stderr, "No key \"type\" in the constituent element\n");
			*err = JLDR_E_BAD_INPUT;
			return ctr;
		}
		if (json_is_string(jTypeId) == 0) {
			fprintf(stderr, "Item \"type\" is not a string\n");
			*err = JLDR_E_BAD_INPUT;
			return ctr;
		}
		typeId = json_string_value(jTypeId);

		jRoleId = json_object_get(item, "role");
		if (jRoleId == NULL) {
			fprintf(stderr, "No key \"role\" in the constituent element\n");
			*err = JLDR_E_BAD_INPUT;
			return ctr;
		}
		if (json_is_string(jRoleId) == 0) {
			fprintf(stderr, "Item \"role\" is not a string\n");
			*err = JLDR_E_BAD_INPUT;
			return ctr;
		}
		roleId = json_string_value(jRoleId);


		jName = json_object_get(item, "name");
		if (jName == NULL) {
			fprintf(stderr, "No key \"name\" in the constituent element\n");
			*err = JLDR_E_BAD_INPUT;
			return ctr;
		}
		if (json_is_string(jName) == 0) {
			fprintf(stderr, "Item \"name\" is not a string\n");
			*err = JLDR_E_BAD_INPUT;
			return ctr;
		}
		name = json_string_value(jName);

		jChargeLow = json_object_get(item, "chargeLow");
		if (jChargeLow == NULL) {
			fprintf(stderr, "No key \"chargeLow\" in the constituent element\n");
			*err = JLDR_E_BAD_INPUT;
			return ctr;
		}
		if (json_is_integer(jChargeLow) == 0) {
			fprintf(stderr, "Item \"chargeLow\" is not an integer\n");
			*err = JLDR_E_BAD_INPUT;
			return ctr;
		}
		chargeLow = json_integer_value(jChargeLow);

		jChargeHigh = json_object_get(item, "chargeHigh");
		if (jChargeHigh == NULL) {
			fprintf(stderr, "No key \"chargeHigh\" in the constituent element\n");
			*err = JLDR_E_BAD_INPUT;
			return ctr;
		}
		if (json_is_integer(jChargeHigh) == 0) {
			fprintf(stderr, "Item \"chargeHigh\" is not an integer\n");
			*err = JLDR_E_BAD_INPUT;
			return ctr;
		}
		chargeHigh = json_integer_value(jChargeHigh);

		if (chargeHigh - chargeLow < 0) {
			fprintf(stderr, "Invalid values of \"n\" or \"p\"\n");
			*err = JLDR_E_BAD_INPUT;
			return ctr;
		}

		jConcentrationBGE = json_object_get(item, "concentrationBGE");
		if (jConcentrationBGE == NULL) {
			fprintf(stderr, "No key \"concentrationBGE\" in constituent element\n");
			*err = JLDR_E_BAD_INPUT;
			return ctr;
		}
		if (is_json_number(jConcentrationBGE) == 0) {
			fprintf(stderr, "Item \"concentrationBGE\" is not a real number\n");
			*err = JLDR_E_BAD_INPUT;
			return ctr;
		}
		concentrationBGE = json_number_value(jConcentrationBGE);

		jConcentrationSample = json_object_get(item, "concentrationSample");
		if (jConcentrationSample == NULL) {
			fprintf(stderr, "No key \"concentrationSample\" in constituent element\n");
			*err = JLDR_E_BAD_INPUT;
			return ctr;
		}
		if (is_json_number(jConcentrationBGE) == 0) {
			fprintf(stderr, "Item \"concentrationSample\" is not a real number\n");
			*err = JLDR_E_BAD_INPUT;
			return ctr;
		}
		concentrationSample = json_number_value(jConcentrationSample);

		jViscosityCoefficient = json_object_get(item, "viscosityCoefficient");
		if (jConcentrationSample == NULL) {
			fprintf(stderr, "No key \"viscosityCoefficient\" in constituent element\n");
			*err = JLDR_E_BAD_INPUT;
			return ctr;
		}
		if (is_json_number(jViscosityCoefficient) == 0) {
			fprintf(stderr, "Item \"viscosityCoefficient\" is not a real number\n");
			*err = JLDR_E_BAD_INPUT;
			return ctr;
		}
		viscosityCoefficient = json_number_value(jViscosityCoefficient);

		jPKaArray = json_object_get(item, "pKas");
		if (jPKaArray == NULL) {
			fprintf(stderr, "No key \"pKas\" in the constituent element\n");
			*err = JLDR_E_BAD_INPUT;
			return ctr;
		}
		if (json_is_array(jPKaArray) == 0) {
			fprintf(stderr, "Item \"pKas\" is not an array\n");
			*err = JLDR_E_BAD_INPUT;
			return ctr;
		}
		if (json_array_size(jPKaArray) < (chargeHigh - chargeLow)) {
			fprintf(stderr, "Array \"pKas\" is too short to cover all ionic forms\n");
			*err = JLDR_E_BAD_INPUT;
			return ctr;
		}

		pKas = parse_array_reals(jPKaArray, chargeHigh - chargeLow);
		if (pKas == NULL && (chargeHigh - chargeLow > 1)) {
			fprintf(stderr, "Cannot parse pKas values\n");
			*err = JLDR_E_BAD_INPUT;
			return ctr;
		}

		jMobilities = json_object_get(item, "mobilities");
		if (jMobilities == NULL) {
			fprintf(stderr, "No key \"mobilities\" in the constituent element\n");
			*err = JLDR_E_BAD_INPUT;
			free(pKas);
			return ctr;
		}
		if (json_is_array(jMobilities) == 0) {
			fprintf(stderr, "Item \"mobilities\" is not an array\n");
			*err = JLDR_E_BAD_INPUT;
			free(pKas);
			return ctr;
		}
		if (json_array_size(jMobilities) < (chargeHigh - chargeLow + 1)) {
			fprintf(stderr, "Array \"mobilities\" is too small to cover all ionic forms\n");
			*err = JLDR_E_BAD_INPUT;
			free(pKas);
			return ctr;
		}

		mobilities = parse_array_reals(jMobilities, chargeHigh - chargeLow + 1);
		if (mobilities == NULL) {
			fprintf(stderr, "Canot parse mobilities values\n");
			*err = JLDR_E_BAD_INPUT;
			free(pKas);
			return ctr;
		}

		init_constituent(ctuent, typeId, roleId, name, chargeLow, chargeHigh, concentrationBGE, concentrationSample, pKas, mobilities, viscosityCoefficient);

		if (ctuent->ctype == NUCLEUS) {
			if (parse_complex_forms(item, ctuent) <= 0) {
				*err = JLDR_E_BAD_INPUT;
				return ctr;
			}
		} else if (ctuent->ctype == LIGAND) {
			ctuent->complexForms = NULL;
			ctuent->complexFormsCount = 0;
		} else {
			*err = JLDR_E_BAD_INPUT;
			return 0;
		}

		ctr++;
	}

	return ctr;
}

static const constituent_array_t * parse_json(const json_t *node, enum LoaderErrorCode *err)
{
	json_t *jCtsArray;
	constituent_t *allCts;
	constituent_array_t *ctArr;
	size_t count;
	size_t realCount;

	jCtsArray = json_object_get(node, "constituents");
	if (jCtsArray == NULL) {
		fprintf(stderr, "No key \"constituents\" was found in the input\n");
		*err = JLDR_E_BAD_INPUT;
		return NULL;
	}
	if (json_is_array(jCtsArray) == 0) {
		fprintf(stderr, "Item \"constituents\" is not an array\n");
		*err = JLDR_E_BAD_INPUT;
		return NULL;
	}
	count = json_array_size(jCtsArray);

	ctArr = malloc(sizeof(constituent_array_t));
	if (ctArr == NULL) {
		*err = JLDR_E_NO_MEM;
		return NULL;
	}

	allCts = calloc(count, sizeof(constituent_t));
	if (allCts == NULL) {
		fprintf(stderr, "Insufficient memory to store all constituents\n");
		free(ctArr);
		*err = JLDR_E_NO_MEM;
		return NULL;
	}

	realCount = parse_cts_array(jCtsArray, allCts, count, err);
	if (*err != JLDR_OK) {
		size_t idx;

		for (idx = 0; idx < realCount; idx++)
			ldr_destroy_constituent(&allCts[idx]);

		free(ctArr);
		return NULL;
	}

	print_all_constituents(allCts, realCount);

	ctArr->constituents = allCts;
	ctArr->count = realCount;

	return ctArr;
}

const constituent_array_t * ldr_loadFromFile(const char *fileName, enum LoaderErrorCode *err)
{
	FILE *f;
	json_t *root;
	json_error_t jsonError;
	const constituent_array_t *ctArr;

	f = fopen(fileName, "r");
	if (f == NULL) {
		fprintf(stderr, "Cannot open file \"%s\" for reading: %s\n", fileName, strerror(errno));
		*err = JLDR_E_CANT_READ;
		return NULL;
	}

	root = json_loadf(f, JSON_REJECT_DUPLICATES, &jsonError);
	if (root == NULL) {
		print_json_error(&jsonError);
		json_decref(root);
		*err = JLDR_E_MALFORMED;
		return NULL;
	}

	ctArr = parse_json(root, err);
	json_decref(root);

	return ctArr;
}

void ldr_destroy_array(const constituent_array_t *array)
{
	free(array->constituents);
	free((void *)array);
}

void ldr_destroy_constituent(constituent_t *ctuent)
{
	size_t idx;

	if (ctuent->complexForms != NULL) {
		for (idx = 0; idx < ctuent->complexFormsCount; idx++) {
			size_t jdx;
			complex_form_t *cForm = &(ctuent->complexForms[idx]);

			if (cForm->ligandGroups == NULL)
				continue;

			for (jdx = 0; jdx < cForm->count; jdx++) {
				size_t kdx;
				ligand_group_t *lGroup = &cForm->ligandGroups[jdx];

				if (lGroup->ligandForms == NULL)
					continue;

				for (kdx = 0; kdx < lGroup->count; kdx++) {
					ligand_form_t *lForm = &lGroup->ligandForms[kdx];

					free(lForm->name);
					free(lForm->pBs);
					free(lForm->mobilities);
				}
				free(lGroup->ligandForms);
			}
			free(cForm->ligandGroups);
		}
		free(ctuent->complexForms);
	}

	free(ctuent->name);
	free(ctuent->pKas);
	free(ctuent->mobilities);
}
