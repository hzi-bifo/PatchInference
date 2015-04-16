#ifndef ANTIPATCH_H_
#define ANTIPATCH_H_

#include <iostream>
#include <vector>
#include <string>
#include <bitset>
#include <limits>
#include <set>
#include "maxflow/graph.h"
#include "SubPatcher.h"
#include "Exception.h" // throws Exception!

// ############################################################################

/**
 * Graph is a type defined in the maxflow package by Yuri Boykov.
 */
typedef Graph<float, float, float> GraphType;

// ############################################################################

/**
 * AdaPatcher is the class used to read in residue information and assign it
 * to source (positive selection) and sink (negative selection) nodes in a
 * graph cut approach. This is done with help of the method findSelection().
 * Afterwards, all selected nodes are grouped to patches based on their distance
 * to each other with the help of the SubPatcher (what an ugly name...).
 *
 * @ todo the greatest weakness of AdaPatch is its parameters. In the graph
 * cut, the p-value (measure for selective pressure) has to be weighted against
 * the distance to neighboring residues on the structure. This parameter has
 * to be set by hand. Also, the SubPatcher groups residues based on their
 * distance. The threshold for merging has to be set, too. Finally, the
 * SubPatcher sorts out patches of low quality (so far: of small size). A
 * measure for the quality of a patch has to be found and a threshold needs to
 * be defined.
 */
class AntiPatcher {
private:
	std::vector<Residue*> _residues; //< input information
	std::set<Residue*> _selection; //< selected residues after calling findSelection()

	float beta;
	float delta;
	bool allowBuried;

	// ---------------------------------------------------------------------------

	void readInResidueInformation(std::string pathToResidueInformation,
			std::string chain);


	std::set<std::string> readListOfResiduesOnSurface(std::string path,
			char separator, std::string token);

	// ---------------------------------------------------------------------------

	float getBestFlow(GraphType* g, Patch R, bool verbose);

	void makeSimpleNetwork(GraphType* g, std::vector<Residue*> nodes,
				float beta);

	size_t findSelection(float beta, bool allowBuried);

	// ---------------------------------------------------------------------------

	void scaleResidueCoords();

	// ---------------------------------------------------------------------------

	float getKernelDensity(std::vector<float> x);
	float getMean(std::vector<float>);
	float getKernel(float f);

	// ---------------------------------------------------------------------------

	std::vector<float> getPairwiseDistances(bool withWeight = false);

public:

	AntiPatcher();
	~AntiPatcher();

	// ---------------------------------------------------------------------------

	void readStructureData(std::string pdbPath, std::string chain = "A");
	void readSurfaceData(std::string surfaceInfoPath, bool allowBuried);
	void readAntigenicWeights(std::string antigenicWeightsPatch);

	// ---------------------------------------------------------------------------

	void setPvalues();

	// ---------------------------------------------------------------------------

	size_t findSelection(bool allowBuried);

	// ---------------------------------------------------------------------------

	std::vector<Residue*> getResidues();
	std::set<Residue*> getSelection();
	std::vector<float> getPairwiseDistance();
	float getBeta();
	float getDelta();
	Patch getResiduesWithSurfaceInformation(bool surface);

	// ---------------------------------------------------------------------------

	void setBeta(float beta);
	void setDelta(float delta);

	// ---------------------------------------------------------------------------

	float makeBetaEvaluation(bool verbose=false);


};

template<typename T>
std::string intToBin(const T & value) {
	std::bitset<std::numeric_limits<T>::digits + 1> bs(value);
	std::string s(bs.to_string());
	return s.substr(s.size() - 4);
}

template<typename T>
std::string NumberToString(T Number) {
	std::ostringstream ss;
	ss << Number;
	return ss.str();
}
// ############################################################################

#endif /* ANTIPATCH_H_ */

