#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <lemng.h>
#include "jsonloader/inputreader.h"
#include "json_input_processor.h"

void printComposition(ECHMET::LEMNG::RConstituentMap *composition) {
	ECHMET::LEMNG::RConstituentMap::Iterator *it = composition->begin();
	if (it == NULL) {
		std::cerr << "Cannot get iterator";
		return;
	}

	while (it->hasNext()) {
		const ECHMET::LEMNG::RConstituent &ctuent = it->value();

		std::cout << "- " << ctuent.name->c_str() << " " << ctuent.concentration << "\n";
		ECHMET::LEMNG::RFormMap::Iterator *fit = ctuent.forms->begin(); /* Here should be a nullptr check */
		while (fit->hasNext()) {
			const ECHMET::LEMNG::RForm &rForm = fit->value();

			std::cout << "\t";

			for (size_t jdx = 0; jdx < rForm.ions->size(); jdx++) {
				const ECHMET::LEMNG::RIon &ion = rForm.ions->at(jdx);

				std::cout << ion.name->c_str() << "(" << ion.charge << ")" << "[" << ion.count << "]";
			}
			std::cout << ": " << rForm.concentration << "\n";

			fit->next();
		}
		fit->destroy();
		it->next();

	}
	it->destroy();
	std::cout << "\n";
}

void printProperties(const ECHMET::LEMNG::RSolutionProperties &properties)
{
	std::cout << "pH = " << properties.pH << std::endl;
	std::cout << "conductivity = " << properties.conductivity << std::endl;
	std::cout << "ionic strength = " << properties.ionicStrength << std::endl;

	printComposition(properties.composition);
}

const char * ezType(const ECHMET::LEMNG::EigenzoneType ezType)
{
	if (ezType == ECHMET::LEMNG::ANALYTE)
		return "(ANALYTE)";
	return "(SYSTEM)";
}

void printResults(const ECHMET::LEMNG::Results &results, const double drivingVoltage, const double totalLength, const double effectiveLength, const double uEOF, const std::map<std::string, double> &aMap)
{
	std::cout << "*** BGE PROPERTIES ***\n";
	printProperties(results.BGEProperties);

	for (size_t idx = 0; idx < results.eigenzones->size(); idx++) {
		const ECHMET::LEMNG::REigenzone &ez = results.eigenzones->at(idx);

		std::cout << "*** EIGENZONE " << idx << " "
			<< ezType(ez.ztype)
			<< (ez.tainted ? " (TAINTED)" : "") << " ***\n";
		std::cout << "mobility: " << ez.mobility << "\n";
		std::cout << "a2t: " << ez.a2t << "\n";
		std::cout << "uEMD: " << ez.uEMD << "\n";

		printProperties(ez.solutionProperties);
	}

	std::cout << "--- Plotting EFG ---\n";

	std::vector<double> times;
	std::vector<std::vector<double> > signals;

	std::cout << "* Conductivity *\n";
	ECHMET::LEMNG::EFGPairVec *electrophoregram;
	ECHMET::LEMNG::RetCode tRet = plotElectrophoregram(electrophoregram, results, drivingVoltage, totalLength, effectiveLength, uEOF, 0.001, ECHMET::LEMNG::RESP_CONDUCTIVITY);
	if (tRet != ECHMET::LEMNG::OK) {
		std::cout << " Cannot plot EFG " << LEMNGerrorToString(tRet) << std::endl;
		return;
	}

	for (size_t idx = 0; idx < electrophoregram->size(); idx++) {
		const ECHMET::LEMNG::EFGPair &p = electrophoregram->at(idx);
		//std::cout << p.time << "; " << p.value << "\n";

		times.push_back(p.time);
		signals.push_back(std::vector<double>(p.value));
	}
	electrophoregram->destroy();

	for (std::map<std::string, double>::const_iterator cit = aMap.begin(); cit != aMap.end(); cit++) {
		const char *key = cit->first.c_str();

		std::cout << "* " << key << " *\n";

		plotElectrophoregram(electrophoregram, results, drivingVoltage, totalLength, effectiveLength, uEOF, 0.001, ECHMET::LEMNG::RESP_CONCENTRATION, key);

		for (size_t idx = 0; idx < electrophoregram->size(); idx++) {
			const ECHMET::LEMNG::EFGPair &p = electrophoregram->at(idx);
			signals[idx].push_back(p.value);
			//std::cout << p.time << "; " << p.value << "\n";
		}
		electrophoregram->destroy();
	}


	std::ofstream efgPlots("efgplots.csv");
	for (size_t idx = 0; idx < times.size(); idx++) {
		efgPlots << times.at(idx) << "; ";

		const std::vector<double> &sig = signals.at(idx);
		for (size_t jdx = 0; jdx < sig.size(); jdx++) {
			efgPlots << sig.at(jdx) << "; ";
		}
		efgPlots << "\n";
	}
}

const char * isEnabledAns(const bool enabled)
{
	if (enabled)
		return "yes";
	return "no";
}

void printTracepointInfo()
{
	ECHMET::LEMNG::TracepointInfoVec *tpiVec = ECHMET::LEMNG::tracepointInfo();
	if (tpiVec == NULL) {
		std::cout << "No tracepoints\n";
		return;
	}

	for (size_t idx = 0; idx < tpiVec->size(); idx++) {
		const ECHMET::LEMNG::TracepointInfo &tpi = tpiVec->at(idx);
		std::cout << "TRACEPOINT " << tpi.id << " " << tpi.description->c_str()
			  << ", enabled: " << isEnabledAns(ECHMET::LEMNG::tracepointState(tpi.id)) << "\n";
	}

	for (size_t idx = 0; idx < tpiVec->size(); idx++)
		tpiVec->at(idx).description->destroy();
	tpiVec->destroy();
}

void printTrace()
{
	ECHMET::FixedString *trace = ECHMET::LEMNG::trace();
	if (trace != NULL) {
		std::cout << trace->c_str();
		trace->destroy();
	} else
		std::cout << "No trace\n";
}

void applyConcentrations(ECHMET::LEMNG::InAnalyticalConcentrationsMap *acMap, const std::map<std::string, double> &rdAcMap)
{
	for (std::map<std::string, double>::const_iterator cit = rdAcMap.begin(); cit != rdAcMap.end(); cit++)
		acMap->item(cit->first.c_str()) = cit->second;
}

int launch(int argc, char **argv)
{
	const char *inputDataFile;
	bool correctForDH;
	bool correctForOF;
	bool correctForVS;
	double drivingVoltage;
	double totalLength;
	double effectiveLength;
	double uEOF;
	ECHMET::LEMNG::JsonInputProcessor::InputDescription inputDesc;
	ECHMET::LEMNG::RetCode ltRet;
	ECHMET::LEMNG::CZESystem *czeSystem = NULL;

	if (argc < 9) {
		std::cout << "Usage: inputFile DH_CORRECTION(number) OF_CORRECTION(number) VS_CORRECTION(number), DrivingVoltage(kV) TotalLength(cm) EffectiveLength(cm) uEOF(U)\n";
		return EXIT_FAILURE;
	}

	inputDataFile = argv[1];

	try {
		correctForDH = std::atoi(argv[2]) >= 1;
		correctForOF = std::atoi(argv[3]) >= 1;
		correctForVS = std::atoi(argv[4]) >= 1;
		drivingVoltage = strtod(argv[5], NULL) * 1000.0;
		totalLength = strtod(argv[6], NULL) / 100.0;
		effectiveLength = strtod(argv[7], NULL) / 100.0;
		uEOF = strtod(argv[8], NULL);
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
		std::cerr << "Cannot create CZESystem" << std::endl;
		ECHMET::SysComp::releaseInputData(inputDesc.BGEComposition);
		ECHMET::SysComp::releaseInputData(inputDesc.SampleComposition);

		return EXIT_FAILURE;
	}

	{
		ECHMET::LEMNG::TracepointInfoVec *tpVec = ECHMET::LEMNG::tracepointInfo();
		ECHMET::LEMNG::toggleAllTracepoints(false);
		ECHMET::LEMNG::toggleTracepoint(14, true);
		ECHMET::LEMNG::toggleTracepoint(15, true);
		ECHMET::LEMNG::toggleTracepoint(16, true);
		ECHMET::LEMNG::toggleTracepoint(tpVec->back().id, true);
		printTracepointInfo();

		for (size_t idx = 0; idx < tpVec->size(); idx++)
			tpVec->at(idx).description->destroy();
		tpVec->destroy();
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
		printResults(results, drivingVoltage, totalLength, effectiveLength, uEOF, inputDesc.SampleConcentrations);

	printTrace();

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

