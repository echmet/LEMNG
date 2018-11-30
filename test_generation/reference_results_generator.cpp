#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <lemng.h>
#include <jsonloader/inputreader.h>
#include <json_input_processor.h>

void printResults(const ECHMET::LEMNG::Results &results, const char *outputFile)
{
	FILE *fh = fopen(outputFile, "w");
	if (fh == NULL)
		std::exit(EXIT_FAILURE);

	fprintf(fh, "%.11g\n", results.BGEProperties.pH);
	fprintf(fh, "%.11g\n", results.BGEProperties.conductivity);
	fprintf(fh, "%.11g\n", results.BGEProperties.ionicStrength);
	fprintf(fh, "%.11g\n", results.BGEProperties.bufferCapacity);

	fprintf(fh, "\n");

	for (size_t idx = 0; idx < results.eigenzones->size(); idx++) {
		const ECHMET::LEMNG::REigenzone &ez = results.eigenzones->at(idx);

		fprintf(fh, "%.11g\n", ez.mobility);
		fprintf(fh, "%.11g\n", ez.uEMD);
		fprintf(fh, "%.11g\n", ez.solutionProperties.pH);
		fprintf(fh, "%.11g\n", ez.solutionProperties.conductivity);
		fprintf(fh, "\n");
	}

	fclose(fh);
}

const char * isEnabledAns(const bool enabled)
{
	if (enabled)
		return "yes";
	return "no";
}

void applyConcentrations(ECHMET::LEMNG::InAnalyticalConcentrationsMap *acMap, const std::map<std::string, double> &rdAcMap)
{
	for (std::map<std::string, double>::const_iterator cit = rdAcMap.begin(); cit != rdAcMap.end(); cit++)
		acMap->item(cit->first.c_str()) = cit->second;
}

int launch(int argc, char **argv)
{
	const char *inputDataFile;
	const char *outputFile;
	bool correctForDH;
	bool correctForOF;
	bool correctForVS;
	ECHMET::LEMNG::JsonInputProcessor::InputDescription inputDesc;
	ECHMET::LEMNG::RetCode ltRet;
	ECHMET::LEMNG::CZESystem *czeSystem = NULL;

	if (argc < 6) {
		std::cout << "Usage: inputFile outputFile DH_CORRECTION(number) OF_CORRECTION(number) VS_CORRECTION(number)\n";
		return EXIT_FAILURE;
	}

	inputDataFile = argv[1];
	outputFile = argv[2];

	std::cout << inputDataFile << "\n" << outputFile << std::endl;

	try {
		correctForDH = std::atoi(argv[3]) >= 1;
		correctForOF = std::atoi(argv[4]) >= 1;
		correctForVS = std::atoi(argv[5]) >= 1;
	} catch (...) { /* THIS DOES NOT WORK WITH PRE-C++11 CONVERSION FUNCS. */
		std::cerr << "ERROR: Invalid input parameters\n";

		return EXIT_FAILURE;
	}

	try {
		ECHMET::LEMNG::JsonInputProcessor inputProc;
		InputReader reader;
		const constituent_array_t *ctarray = reader.read(inputDataFile);
		inputDesc = inputProc.process(ctarray);
	} catch (InputReader::MalformedInputException &ex) {
		std::cerr << ex.what();
		return EXIT_FAILURE;
	}

	ltRet = ECHMET::LEMNG::makeCZESystem(inputDesc.BGEComposition, inputDesc.SampleComposition, czeSystem);
	if (ltRet != ECHMET::LEMNG::OK) {
		std::cerr << "Cannot create CZESystem: " << ECHMET::LEMNG::LEMNGerrorToString(ltRet) << std::endl;
		ECHMET::SysComp::releaseInputData(inputDesc.BGEComposition);
		ECHMET::SysComp::releaseInputData(inputDesc.SampleComposition);

		return EXIT_FAILURE;
	}

	ECHMET::LEMNG::InAnalyticalConcentrationsMap *acBGEMap;
	ECHMET::LEMNG::InAnalyticalConcentrationsMap *acFullMap;

	if (czeSystem->makeAnalyticalConcentrationsMaps(acBGEMap, acFullMap) != ECHMET::LEMNG::OK) {
		std::cerr << "Failed to get analytical concentration maps" << std::endl;
		return EXIT_FAILURE;
	}

	applyConcentrations(acBGEMap, inputDesc.BGEConcentrations);
	applyConcentrations(acFullMap, inputDesc.SampleConcentrations);

	ECHMET::NonidealityCorrections corrections = ECHMET::defaultNonidealityCorrections();
	if (correctForDH)
		ECHMET::nonidealityCorrectionSet(corrections, ECHMET::CORR_DEBYE_HUCKEL);
	if (correctForOF)
		ECHMET::nonidealityCorrectionSet(corrections, ECHMET::CORR_ONSAGER_FUOSS);
	if (correctForVS)
		ECHMET::nonidealityCorrectionSet(corrections, ECHMET::CORR_VISCOSITY);

	ECHMET::LEMNG::Results results;
	ECHMET::LEMNG::RetCode tRet = czeSystem->evaluate(acBGEMap, acFullMap, corrections, results);

	if (tRet != ECHMET::LEMNG::OK)
		std::cout << "Failed to solve the system: " << czeSystem->lastErrorString() << "\n";
	else
		printResults(results, outputFile);

	acBGEMap->destroy();
	acFullMap->destroy();
	ECHMET::LEMNG::releaseCZESystem(czeSystem);
	ECHMET::LEMNG::releaseResults(results);

	ECHMET::SysComp::releaseInputData(inputDesc.BGEComposition);
	ECHMET::SysComp::releaseInputData(inputDesc.SampleComposition);

	return EXIT_SUCCESS;
}

int main(int argc, char **argv)
{
	return launch(argc, argv);
}

