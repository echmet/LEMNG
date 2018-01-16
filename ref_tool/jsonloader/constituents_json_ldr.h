#ifndef CONSTITUENTS_LDR_H
#define CONSTITUENTS_LDR_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>

enum ConstituentType {
	LIGAND,
	NUCLEUS,
	INVALID_TYPE
};

enum ConstituentRole {
	BACKGROUND,
	ANALYTE,
	INVALID_ROLE
};

enum LoaderErrorCode {
	JLDR_OK,
	JLDR_E_BAD_INPUT,
	JLDR_E_CANT_READ,
	JLDR_E_MALFORMED,
	JLDR_E_NO_MEM
};

struct ligand_form {
	char *name;
	int charge;
	int maxCount;
	double *pBs;
	double *mobilities;
};
typedef struct ligand_form ligand_form_t;

struct ligand_group {
	size_t count;
	ligand_form_t *ligandForms;
};
typedef struct ligand_group ligand_group_t;

struct complex_form {
	int nucleusCharge;
	size_t count;
	ligand_group_t *ligandGroups;
};
typedef struct complex_form complex_form_t;

struct constituent {
	enum ConstituentType ctype;
	enum ConstituentRole crole;
	char *name;
	int chargeLow;
	int chargeHigh;
	double concentrationBGE;
	double concentrationSample;
	double *pKas;
	double *mobilities;
	double viscosityCoefficient;

	size_t complexFormsCount;
	complex_form_t *complexForms;

};
typedef struct constituent constituent_t;

struct constituent_array {
	constituent_t *constituents;
	size_t count;
};
typedef struct constituent_array constituent_array_t;

const constituent_array_t * ldr_loadFromFile(const char *fileName, enum LoaderErrorCode *err);
void ldr_destroy_array(const constituent_array_t *array);
void ldr_destroy_constituent(constituent_t *ctuent);

#ifdef __cplusplus
}
#endif

#endif // CONSTITUENTS_LDR_H
