/*
 * main.cpp
 *
 *  Created on: Jun 10, 2012
 *      Author: ctusche
 */

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <algorithm>
#include <iomanip>
#include "AntiPatch.h"
#include "Printer.h"

using namespace std;

struct Parameter {

	std::string pfdFilePath, projektName, pdbFilePath, pictureName,
			antigenicWeightsFilePath, surfaceInfoPath, chain;

	int minPatchSize;
	int antigenicWeightsThresh;
	int detail;

	enum kindOfAverage {
		mean = 3, median
	} avgKind;

	bool printSummary, onlySurface, printPymol, printOnKonsole,
			useAntigeniWeights, printPythonString, printAValue;

	float betaBreak, maxBeta, useDelta; // d$x[which.max(d$y)]*0.95

};

Parameter getParameterFromProgramCall(int, char**);

Parameter readConfigFile(std::string filename, char* prog, Parameter par);

void switchParameter(Parameter &par, char* prog, std::string p, char s);

string getHelp(std::string progname, Parameter par);

int main(int argc, char* argv[]) {
	Parameter par;

	if (argc > 1)
		par = getParameterFromProgramCall(argc, argv);

	try {

		AntiPatcher AP;

		// read in structure data
		// this will also set delta
		AP.readStructureData(par.pdbFilePath, par.chain);

		// read in info saying which residue is buried and which is exposed
		AP.readSurfaceData(par.surfaceInfoPath, !par.onlySurface);

		// read antigenic data
		AP.readAntigenicWeights(par.antigenicWeightsFilePath);

		// transform agw's to p-values
		AP.setPvalues();

		// get beta
		AP.setBeta(AP.makeBetaEvaluation());


		// make graph cut to get selection
		AP.findSelection(!par.onlySurface);

		// merging cluster
		SubPatcher SP = SubPatcher(AP.getSelection(), AP.getResidues());

		// set merging to half of AP patch size
		SP.setDelta(0.5 * AP.getDelta());
		PatchResultInfo pri = SP.patchWithInfo();

		Printer p(SP.getPatches(), AP.getResidues(), par.pdbFilePath,
				par.detail);
		p.printSimple();

		cout << "done." << endl;


	} catch (Exception &e) {
		cerr << e.error << endl;
		return -1;
	} catch (const char* msg) {
		cerr << msg << endl;
		return -1;
	}

	return 0;
}

Parameter getParameterFromProgramCall(int argc, char* argv[]) {
	Parameter par;
	std::string test;
	std::string::iterator it;

	for (int i = 1; i < argc; i++) {
		test = argv[i];
		if (*(it = test.begin()) == '-')
			if (i + 1 == argc)
				switchParameter(par, argv[0], "fail", *(++it)); //help?
			else
				switchParameter(par, argv[0], argv[++i], *(++it));
		else
			return readConfigFile(test, argv[0], par);
	}

	return par;
}

void switchParameter(Parameter &par, char* prog, std::string p, char s) {
	std::string test;

	switch (s) {
	case 'h':
		std::cout << getHelp(prog, par);
		exit(1);

	case 'a':
		test = p;
		if (test.compare("mean") == 0)
			par.avgKind = Parameter::mean;
		else if (test.compare("median") == 0)
			par.avgKind = Parameter::median;
		else
			throw Exception(getHelp(prog, par));
		break;

	case 'i':
		par.pfdFilePath = p;
		break;

	case 'n':
		par.projektName = p;
		break;

	case 'c':
		par.chain = p;
		break;

	case 'g':
		par.antigenicWeightsFilePath = p;
		break;

	case 's':
		par.surfaceInfoPath = p;
		break;

	case 'p':
		par.minPatchSize = atoi(p.c_str());
		break;

	case 't':
		par.antigenicWeightsThresh = atoi(p.c_str());
		break;

	case 'x':
		test = p;
		if (test.compare("false") == 0)
			par.printPythonString = false;
		else if (test.compare("true") == 0)
			par.printPythonString = true;
		else
			throw Exception(getHelp(prog, par));
		break;

	case 'y':
		test = p;
		if (test.compare("false") == 0)
			par.printAValue = false;
		else if (test.compare("true") == 0)
			par.printAValue = true;
		else
			throw Exception(getHelp(prog, par));
		break;

	case 'z':
		par.betaBreak = atof(p.c_str());
		break;

	case 'd':
		test = p;
		if (test.compare("false") == 0)
			par.printSummary = false;
		else if (test.compare("true") == 0)
			par.printSummary = true;
		else
			throw Exception(getHelp(prog, par));
		break;

	case 'u':
		test = p;
		if (test.compare("false") == 0)
			par.onlySurface = false;
		else if (test.compare("true") == 0)
			par.onlySurface = true;
		else
			throw Exception(getHelp(prog, par));
		break;

	case 'm':
		test = p;
		if (test.compare("false") == 0)
			par.printPymol = false;
		else if (test.compare("true") == 0)
			par.printPymol = true;
		else
			throw Exception(getHelp(prog, par));
		break;

	case 'b':
		par.pdbFilePath = p;
		break;

	case 'r':
		par.pictureName = p;
		break;

	case 'v':
		par.maxBeta = atof(p.c_str());
		break;

	case 'e':
		par.detail = atoi(p.c_str());
		break;

	case 'f':
		par.useDelta = atof(p.c_str());
		break;

	case 'k':
		test = p;
		if (test.compare("false") == 0)
			par.printOnKonsole = false;
		else if (test.compare("true") == 0)
			par.printOnKonsole = true;
		else
			throw Exception(getHelp(prog, par));
		break;

	default:
		Exception(getHelp(prog, par));
	}
}

Parameter readConfigFile(std::string filename, char* prog, Parameter par) {
	ifstream in(filename.c_str(), ios::in);
	if (in.fail()) {
		cout << "could not open File: '" << filename << "'!" << endl;
		exit(1);
	}
	char tmp[256];
	std::string p1, p2;

	while (!in.eof()) {
		in.getline(tmp, sizeof(tmp));
		stringstream ss(tmp);
		ss >> p1 >> p2;
		if (p1[0] == '-')
			switchParameter(par, prog, p2, p1[1]);
	}

	return par;
}

string getHelp(std::string progname, Parameter par) {
	stringstream help;
	help << boolalpha << "Usage: " << progname << std::endl
			<< "-h < get this Help >" << std::endl << std::endl
			<< "-a < mean, median >" << " (default): "
			<< ((par.avgKind == Parameter::mean) ? "mean" : "median")
			<< std::endl

//          << "-i < path to pfd-file >" << " (default): " << par.pfdFilePath << std::endl
			<< "-n < project name >" << " (default): " << par.projektName
			<< std::endl << "-c < used chain >" << " (default): " << par.chain
			<< std::endl << "-f < used delta >" << " (default): "
			<< par.useDelta << std::endl
			<< "-g < path to antigenic weights-file >" << " (default): "
			<< par.antigenicWeightsFilePath << std::endl
			<< "-s < path to structure file >" << " (default): "
			<< par.surfaceInfoPath << std::endl << "-p < minimum patch-size >"
			<< " (default): " << par.minPatchSize << std::endl
			<< "-x < print Values to use python Skript >" << " (default): "
			<< par.printPythonString << std::endl
			<< "-y < print aValues for agw Plot >" << " (default): "
			<< par.printAValue << std::endl
			<< "-z < threshhold for betasearch >" << " (default): "
			<< par.betaBreak << std::endl
			<< "-t < used thresh for antigenic weights >" << " (default): "
			<< par.antigenicWeightsThresh << std::endl
//          << "-d < print Summary in txt-File? >" << " (default): " << par.printSummary << std::endl
			<< "-u < use only surface data? >" << " (default): "
			<< par.onlySurface << std::endl
			<< "-m < print pymol skript (both)? >" << " (default): "
			<< par.printPymol << std::endl << "-b < path to pdb file >"
			<< " (default): " << par.pdbFilePath << std::endl
			<< "-r < picture name (will be extended) >" << " (default): "
			<< par.pictureName << std::endl
//          << "-e < size of detail >" << " (default): " << par.detail << std::endl
			<< "-k < print on Console? >" << " (default): "
			<< par.printOnKonsole << std::endl
//          << "-w < use antigenic weights? >" << " (default): " << par.useAntigeniWeights << std::endl
//          << "-v < maximum beta value >" << " (default): " << par.maxBeta << std::endl
			<< endl << endl;

	return help.str();
}

