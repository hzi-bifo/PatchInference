/*
 * SubPatcher.cpp
 *
 *  Created on: Jun 10, 2012
 *      Author: ctusche
 */

#include "SubPatcher.h"
#include <math.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <cfloat>
#define _USE_MATH_DEFINES

SubPatcher::SubPatcher(std::set<Residue*> selection, Patch resid) {
	this->delta = 0.0;
	this->R = resid;
	this->allResidues = resid;
	this->S = selection;
	for (size_t i = 0; i < resid.size(); i++) {
		Patch b;
		if (S.find(resid.at(i)) != S.end()) {
			//std::cerr << "  selected " << resid.at(i)->getAlias() << "\t" << resid.at(i)->getAgW() << std::endl;

			b.push_back(resid[i]);
			this->patches.push_back(b);
			R[i]->setPatched(true);
		}
	}
}

void SubPatcher::patch() {

	PatchList::iterator i, j;
	std::pair<PatchList::iterator, PatchList::iterator> r;

	i = patches.begin();

	while (true) {

		// get the pair of patches that have the smallest distance
		// to each other
		r = getClosestPair(); // returns <end, end> if > thresh

		i = r.first;
		j = r.second;

		// stop if nothing was found --> nothing to merge anymore
		if (i == patches.end())
			break;

		// else merge the patches == copy everything from patch j to
		// patch i, then delete j
		merge(i, j); // does nothing in case of <end, end>

		if (j != patches.end()) {
			patches.erase(j);
		}
	}

	for (i = patches.begin(); i != patches.end(); i++) {
		if (i->size() == 1) {
			patches.erase(i);
		}
	}

}

PatchResultInfo makePri() {
	PatchResultInfo pri;
	pri.avgWeightNoPatch = 0;
	pri.avgWeightPatches = 99999;
	pri.maxDistance = 0;
	pri.totalSize = 0;
	pri.averageDistance = 0;
}

float getSizeOfPatch(Patch p) {
	Patch::iterator a, b;
	float largestDistance = 0;
	float curDistance;
	for (a = p.begin(); a != p.end(); a++) {
		for (b = p.begin(); b != p.end(); b++) {
			curDistance = sqrt((*a)->getDistance(*b));
			if (curDistance > largestDistance) {
				largestDistance = curDistance;
			}
		}
	}
	return largestDistance;
}

#include <set>
PatchResultInfo SubPatcher::patchWithInfo() {

	PatchResultInfo pri = makePri();

	PatchList::iterator i, j;
	std::pair<PatchList::iterator, PatchList::iterator> r;

	i = patches.begin();

	while (true) {

		// get the pair of patches that have the smallest distance
		// to each other
		r = getClosestPair(); // returns <end, end> if > thresh

		i = r.first;
		j = r.second;

		// stop if nothing was found --> nothing to merge anymore
		if (i == patches.end())
			break;

		// else merge the patches == copy everything from patch j to
		// patch i, then delete j
		merge(i, j); // does nothing in case of <end, end>

		if (j != patches.end()) {
			patches.erase(j);
		}
	}

	size_t MIN_SIZE = 0;
	for (i = patches.begin(); i != patches.end();) {
		if (i->size() <= MIN_SIZE) {
			patches.erase(i);
		} else {
			++i;
		}
	}

	PatchList::iterator it;

	float weight = 0;
	size_t totalRes = 0;
	std::set<Residue*> resInPatch;

	float largestDistance = 0;
	float curDistance;
	float averageDistance = 0;

	for (it = patches.begin(); it != patches.end(); it++) {
		for (size_t r = 0; r < it->size(); r++) {
			resInPatch.insert(it->at(r));
			weight += it->at(r)->getAgW();
			totalRes++;
		}
		curDistance = getSizeOfPatch(*it);
		averageDistance += curDistance;
		if (curDistance > largestDistance) {
			largestDistance = curDistance;
		}
	}

	pri.maxDistance = largestDistance;
	pri.numberOfPatches = patches.size();

	if (patches.size() > 0) {
		pri.averageDistance = averageDistance / float(patches.size());
	}

	pri.totalSize = totalRes;
	if (totalRes > 0) {
		pri.avgWeightPatches = weight / float(totalRes);
	} else {
		pri.avgWeightPatches = 9999;
	}

	weight = 0;
	totalRes = 0;
	std::vector<Residue*>::iterator iter;

	for (iter = this->allResidues.begin(); iter != this->allResidues.end();
			iter++) {
		// if the residue is NOT in a patch
		if (resInPatch.find(*iter) == resInPatch.end()) {
			//weight += (*iter)->getAgW();
			//this is the 'inverse' antigenic weight, i.e. max(a(n)) - a(n)
			weight += (*iter)->getPValue();
			totalRes++;
		}
	}

	if (totalRes > 0) {
		pri.avgWeightNoPatch = weight / float(totalRes);
	}

	return pri;

}

void SubPatcher::patchIntermediate() {
	// init iter
	Patch::iterator it, kt;
	PatchList::iterator jt = this->patches.begin();
	bool inter = false;

	// all residues
	for (it = this->R.begin(); it != this->R.end(); ++it) {
		// except those which are patched already
		if (!(*it)->isPatched()) {
			for (; jt != this->patches.end(); ++jt) {
				// omly if there are big enough to give space for others
				if (jt->size() > 2) {
					inter = true;
					//all residues from this patch
					for (kt = jt->begin(); kt != jt->end(); ++kt) {
						// only those, which are close anough to all other
						if ((*it)->getDistance(*kt) > this->delta) {
							inter = false;
						}
					}
					if (inter) {
						jt->push_back(*it);
						(*it)->setPatched(true);
						intermediate.push_back(*it);
					}
				}
			}
		}
	}
}

bool vergleich(Residue *a, Residue *b) {
	return a->getPos() < b->getPos();
}

Patch SubPatcher::getSelectedResidues() {
	Patch res;
	for (PatchList::iterator it = this->patches.begin();
			it != this->patches.end(); ++it)
		for (Patch::iterator jt = it->begin(); jt != it->end(); ++jt)
			res.push_back(*jt);

	std::sort(res.begin(), res.end(), vergleich);
	return res;
}


void SubPatcher::merge(PatchList::iterator a, PatchList::iterator b) {

	// do this only if a and b are reasonable patches
	if ((a != patches.end()) & (b != patches.end())) {
		for (size_t k = 0; k < b->size(); k++) {
			a->push_back(b->at(k));
		}
	}
}

float SubPatcher::getEuklidMeanDist(PatchList::iterator a,
		PatchList::iterator b) {
	float meanDist = 0;

	for (size_t i = 0; i < a->size(); i++)
		for (size_t j = 0; j < b->size(); j++)
			meanDist += sqrt(a->at(i)->getDistance(b->at(j)));

	return (meanDist / float(a->size() * b->size()));
}

float SubPatcher::getEuklidMinDist(PatchList::iterator a,
		PatchList::iterator b) {

	float tmp, minDist = FLT_MAX;

	for (size_t i = 0; i < a->size(); i++) {
		for (size_t j = 0; j < b->size(); j++) {
			tmp = sqrt(a->at(i)->getDistance(b->at(j)));
			if (tmp < minDist) {
				minDist = tmp;
			}
		}
	}

	return minDist;
}

std::pair<PatchList::iterator, PatchList::iterator> SubPatcher::getClosestPair() {

	PatchList::iterator i, j;
	PatchList::iterator bestI, bestJ;
	float tmp, bestDist = 99999999.9f;

	for (i = patches.begin(); i != patches.end(); i++)
		for (j = patches.begin(); j != i; j++)
			if ((tmp = getEuklidMinDist(i, j)) < bestDist) {
				bestDist = tmp;
				bestI = i;
				bestJ = j;
			}

	std::pair<PatchList::iterator, PatchList::iterator> res;
	if (bestDist < delta) {
		res.first = bestI;
		res.second = bestJ;
	} else {
		res.first = patches.end();
		res.second = patches.end();
	}

	return res;
}

void SubPatcher::setDelta(float d) {
	this->delta = d;
}

void SubPatcher::sortAvg(bool scored, bool useAgu) {
	float avg;
	for (PatchList::iterator it = patches.begin(); it != patches.end(); ++it) {
		avg = 0.0;
		for (Patch::iterator jt = it->begin(); jt != it->end(); ++jt)
			if (!useAgu)
				avg += (*jt)->getPValue();
			else
				avg += (*jt)->getAgW();
		if (scored)
			avg = -(1 / (0.5 + avg) + (float(it->size())));
		this->sortedPatches[avg] = *it;
	}
}
