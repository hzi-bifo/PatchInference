/*
 * SubPatcher.h
 *
 *  Created on: Jun 10, 2012
 *      Author: ctusche
 */
#ifndef SUBPATCHER_H_
#define SUBPATCHER_H_

#include <map>
#include <set>
#include <vector>
#include <string>
#include "Residue.h"

typedef std::vector<Residue*> Patch;
typedef std::vector<Patch> PatchList;
typedef std::map<float, Patch> PatchBook;

struct PatchResultInfo {
	size_t totalSize;
	float avgWeightPatches;
	float avgWeightNoPatch;
	float maxDistance;
	float averageDistance;
	size_t numberOfPatches;
};

class SubPatcher {
public:
	SubPatcher(std::set<Residue*> selection, Patch resid);

	SubPatcher() {}

	/**
	 * runs through the list of all patches (initialized: the list of all residues)
	 * and merges them if their distance is smaller than thresh
	 */
	void patch();

	PatchResultInfo patchWithInfo();

	void setDelta(float);

	PatchList getPatches() {
		return this->patches;
	}
	;

	Patch getIntermediates() {
		return this->intermediate;
	}

	Patch getSelectedResidues();

	/**
	 * find intermediate Residues
	 */
	void patchIntermediate();

	std::vector<int> getSortedPatchedResidue();

private:
	std::set<Residue*> S;
	Patch R;
	Patch allResidues;
	float delta;

	PatchList patches;
	PatchList ignore;
	PatchBook sortedPatches;
	Patch intermediate;

	/**
	 * run through all elements of list a and b, if only one pair of elements is
	 * closer than the threshold, insert all elements of b to a and
	 * delete b
	 */
	void merge(PatchList::iterator a, PatchList::iterator b);

	/** computes the average euklid distance between all pairs of residues
	 * from patches a and b
	 */
	float getEuklidMeanDist(PatchList::iterator a, PatchList::iterator b);

	/** computes the minimum euklidian distance between all pairs of residues
	 * from patches a and b
	 */
	float getEuklidMinDist(PatchList::iterator a, PatchList::iterator b);

	/**
	 * given a list of patches (stored in SubPatcher instance), returns the
	 * pair of patches that has the smallest residue distance to each other
	 * @see getMeanDist(PatchList::iterator a, PatchList::iterator b)
	 */
	std::pair<PatchList::iterator, PatchList::iterator> getClosestPair();

	/**
	 *
	 * @param withWeight set True if the weight has to be regard
	 * @return
	 */
	//vector<float> getPairwiseDistances(bool withWeight = false);
	/**
	 * sort Patches with avg
	 */
	void sortAvg(bool scored = false, bool useAgu = false);

};

#endif /* SUBPATCHER_H_ */
