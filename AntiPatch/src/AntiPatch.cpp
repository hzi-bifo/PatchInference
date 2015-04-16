#include <sstream>
#include <iostream>
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <set>

#include "AntiPatch.h"
#include "Table.h"

AntiPatcher::AntiPatcher() {

	delta = 0;
}

std::vector<std::vector<std::string> > readPDB(std::string path,
		std::string chain) {
	std::ifstream in(path.c_str());
	std::vector<std::string> parts;
	std::vector<std::vector<std::string> > res;
	std::string testAtom, resInf, testCarbon, testChain, nextChain(1,
			((char) ((int) chain[0] + 1)));
	bool stop = false;
	int N = 512;
	char l[N];

	if (in.fail()) {
		std::string msg("readPDB: unable to open file \n\t");
		msg.append(path);
		throw Exception(msg);
	}

	while (in.good() || stop) {
		in.getline(l, N);
		std::stringstream ss(l);
		ss >> testAtom;
		if (testAtom.compare("ATOM") == 0) {
			ss >> testCarbon;
			ss >> testCarbon;
			if (testCarbon.compare("CA") == 0) {
				ss >> testChain;
				ss >> testChain;
				if (testChain.compare(chain) == 0) {
					ss >> resInf;
					parts.push_back(resInf); // residue
					for (int i = 0; i < 3; ++i) {
						ss >> resInf;
						parts.push_back(resInf); // x,y,z
					}
				} else if (testChain.compare(chain) == 0)
					stop = true;
			}
		}
		if (!parts.empty())
			res.push_back(parts);
		parts.clear();
	}

	return res;
}

void AntiPatcher::readStructureData(std::string pdbPath, std::string chain) {

	std::vector<std::vector<std::string> > pdbTable = readPDB(pdbPath, chain);

	std::vector<std::string> row;
	std::string pos, alias, epi, oldPos = "";
	Point coord;

	for (size_t i = 0; i < pdbTable.size(); i++) {

		row = pdbTable.at(i);

		if (row.size() < 4) {
			throw "ERROR in AntiPatcher::readStructureData: invalid table row!\n";
		}

		epi = "";
		alias = "";
		pos = row.at(0);

		coord.setPoint(atof(row[1].c_str()), atof(row[2].c_str()),
				atof(row[3].c_str()));
		alias = pos;
		epi = pos;

		//pValue becomes -1 cause value is not known until now
		if (pos != oldPos) {
			_residues.push_back(new Residue(pos, coord, -1.0, alias, epi));
		}
		oldPos = pos;
	}

	scaleResidueCoords();
}

void AntiPatcher::scaleResidueCoords() {
	Patch::iterator it = _residues.begin();

	float maxX = 0, maxY = 0, maxZ = 0;
	float minX = 1000, minY = 1000, minZ = 1000;
	for (it = _residues.begin(); it != _residues.end(); ++it) {
		minX = ((*it)->getCoord().getX() < minX) ?
				(*it)->getCoord().getX() : minX;
		minY = ((*it)->getCoord().getY() < minY) ?
				(*it)->getCoord().getY() : minY;
		minZ = ((*it)->getCoord().getZ() < minZ) ?
				(*it)->getCoord().getZ() : minZ;
		maxX = ((*it)->getCoord().getX() > maxX) ?
				(*it)->getCoord().getX() : maxX;
		maxY = ((*it)->getCoord().getY() > maxY) ?
				(*it)->getCoord().getY() : maxY;
		maxZ = ((*it)->getCoord().getZ() > maxZ) ?
				(*it)->getCoord().getZ() : maxZ;
	}

	//float max = ((maxX - minX) > (maxY - minY)) ? (maxX - minX) : (maxY - minY);
	//max = (max > (maxZ - minZ)) ? max : (maxZ - minZ);

	float max = (maxX > maxY) ? maxX : maxY;
	max = (max > maxZ) ? max : maxZ;

	for (it = _residues.begin(); it != _residues.end(); ++it) {
		(*it)->setCoord((*it)->getCoord().getX() / max,
				(*it)->getCoord().getY() / max, (*it)->getCoord().getZ() / max);
	}

	// set delta to half of largest epitope site on HA of H3N2
	this->delta = 21.95 / max;

}

std::vector<float> AntiPatcher::getPairwiseDistance() {
	std::vector<float> dist;
	Patch::iterator it, jt;

	for (it = this->_residues.begin(); it != this->_residues.end(); ++it)
		for (jt = this->_residues.begin(); jt != it; ++jt)
			dist.push_back((*it)->getEuklideanDistance(*jt));

	return dist;
}

// -----------------------------------------------------------------------------

void AntiPatcher::setPvalues() {
	Patch::iterator it;

	float max = 0;
	for (Patch::iterator it = _residues.begin(); it != _residues.end(); ++it) {
		if ((*it)->getAgW() > max) {
			max = (*it)->getAgW();
		}
	}

	float aw;
	for (it = this->_residues.begin(); it != this->_residues.end(); it++) {
		aw = (*it)->getAgW();
		(*it)->setPValue((max - aw)/max);
	}
}

// -----------------------------------------------------------------------------

void AntiPatcher::readAntigenicWeights(std::string antigenicWeightsPath) {

	std::cerr << antigenicWeightsPath << std::endl;

	Table antigenicInfo(antigenicWeightsPath, "\t");
	std::vector<std::string> row;
	std::string pos;
	float mean, median;

	std::map<std::string, float> weightMap;

	for (size_t r = 0; r < antigenicInfo.rows(); r++) {
		row = antigenicInfo.row(r);
		if (row.size() != 5) {
			continue;
		}
		pos = row.at(0);
		mean = atof(row.at(3).c_str());
		median = atof(row.at(4).c_str());

		weightMap[pos] = mean;
	}

	std::map<std::string, float>::iterator wInf;
	// run through all residues and update their agw
	for (Patch::iterator it = _residues.begin(); it != _residues.end(); ++it) {
		// do we have an weightMap entry for this residue?
		wInf = weightMap.find((*it)->getPos());
		// if yes, assing the agw to this residue
		if (wInf != weightMap.end())
			(*it)->setAgW(wInf->second);
		// if not, set agw to zero
		else
			(*it)->setAgW(0.0f);
	}

}

void AntiPatcher::readSurfaceData(std::string surfaceInfoPath,
		bool allowBuried) {
	std::set<std::string> positionsOnSurface = readListOfResiduesOnSurface(
			surfaceInfoPath, '\t', "TRUE");
	for (Patch::iterator it = _residues.begin(); it != _residues.end(); ++it) {
		if (positionsOnSurface.find((*it)->getPos())
				!= positionsOnSurface.end()) // is there a position with this nr?
			(*it)->setSurfaceInfo(true);
		else
			(*it)->setSurfaceInfo(false);
	}

	this->allowBuried = allowBuried;

}

// -----------------------------------------------------------------------------

size_t AntiPatcher::findSelection(bool allowBuried) {

	return this->findSelection(this->beta, allowBuried);
}

size_t AntiPatcher::findSelection(float beta, bool allowBuried) {

	_selection.clear();

	// make sure we only use the nodes we need
	std::vector<Residue*> nodes;
	if (allowBuried) {
		nodes = _residues;
	} else {
		for (Patch::iterator it = _residues.begin(); it != _residues.end();
				++it) {
			if ((*it)->isOnSurface()) {
				nodes.push_back((*it));
			}
		}
	}

	size_t N = nodes.size();

	GraphType* g = new GraphType(N + 2, N + 50);
	makeSimpleNetwork(g, nodes, beta);
	g->maxflow();

	for (size_t i = 0; i < N; i++) {
		if (g->what_segment(i) != GraphType::SOURCE) {
			_selection.insert(nodes.at(i));
		}
	}

	delete g;
	return _selection.size();
}

// ---------------------------------------------------------------------------
void AntiPatcher::makeSimpleNetwork(GraphType* g, std::vector<Residue*> nodes,
		float beta) {

	int k = 0;

	// -------------- add t-links ---------------------------------
	for (Patch::iterator it = nodes.begin(); it != nodes.end(); ++it) {
		g->add_node();
		float a_bar_n = (*it)->getPValue();
		float a_n = (*it)->getAgW();

		g->add_tweights(k++, beta * a_bar_n, (beta * a_n));
	}
	// -------------- add n-links ---------------------------------
	float dist, weight;
	for (size_t i = 0; i < nodes.size(); i++) {
		for (size_t j = (i + 1); j < nodes.size(); j++) {
			if (i == j)
				continue;

			dist = nodes[i]->getEuklideanDistance(nodes[j]);
			weight = exp(-(dist));

			if (dist < this->delta)
				g->add_edge(i, j, weight, weight);
		}
	}
}

void AntiPatcher::setDelta(float delta) {
	this->delta = delta;
}

// ---------------------------------------------------------------------------

AntiPatcher::~AntiPatcher() {
	for (size_t i = 0; i < _residues.size(); ++i)
		delete _residues[i];
}

// ---------------------------------------------------------------------------

std::set<std::string> AntiPatcher::readListOfResiduesOnSurface(std::string path,
		char separator, std::string token) {
	size_t N = 256;
	char l[N];
	std::fstream in;
	std::string line;
	std::set<std::string> res;

	in.open(path.c_str());
	if (!in.is_open()) {
		std::string msg("getTrueResList: unable to open file \n\t");
		msg.append(path);
		throw Exception(msg);
	}

	while (in.good()) {
		in.getline(l, N);
		line = l;

		size_t sem = line.find(';'), sep = line.find_first_not_of(separator),
				tok = line.find(token);

		if (!(sem != std::string::npos && sem <= sep)
				&& tok != std::string::npos)
			res.insert(line.substr(0, line.find(separator)));
	}

	in.close();
	return res;
}

// ---------------------------------------------------------------------------

Patch AntiPatcher::getResiduesWithSurfaceInformation(bool surface) {
	Patch res;
	for (Patch::iterator it = _residues.begin(); it != _residues.end(); ++it)
		if (!surface || (*it)->isOnSurface())
			res.push_back(*it);
	return res;
}

// ---------------------------------------------------------------------------

Patch AntiPatcher::getResidues() {
	return this->_residues;
}

// ---------------------------------------------------------------------------

std::set<Residue*> AntiPatcher::getSelection() {
	return this->_selection;
}

// ---------------------------------------------------------------------------

float AntiPatcher::makeBetaEvaluation(bool verbose) {

	float stopBeta = 5000;
	float incr = 5;
	float beta = 2;

	PatchResultInfo pri;

	size_t stableCount = 0;
	float oldWeight = 0;
	float lastBeta = 0;


	if (verbose) {
		std::cerr
				<< "beta\tpatch_size\tavg_weight\tavg_out_weight\tmax_distance\taverageDistance\tnumberOfPatches\n";
	}

	while (beta < stopBeta) {

		findSelection(beta, this->allowBuried);

		SubPatcher SP = SubPatcher(this->getSelection(), this->getResidues());
		SP.setDelta(0.5 * this->getDelta());
		pri = SP.patchWithInfo();

		/*if (pri.avgWeightPatches * 0.9 < pri.avgWeightNoPatch) {
			//if (pri.avgWeightPatches - 0.2 < pri.avgWeightNoPatch) {
			if (verbose) {
				std::cerr << "weight in  patches:\t" << pri.avgWeightPatches
						<< std::endl;
				std::cerr << "weight w/o patches:\t" << pri.avgWeightNoPatch
						<< std::endl;
				std::cerr
						<< "Stopped beta search due to saturation criterion at beta = "
						<< beta << "\n";
			}
			break;
		}*/

		if (pri.avgWeightPatches == oldWeight) {
			// no change
			if (stableCount > 1500) {
				// set beta to the last value for which we observed a relevant change
				beta = lastBeta;
				if (verbose) {
					std::cerr
							<< "Stopped beta search due to max count criterion at beta = "
							<< beta << "\n";
				}
				break;
			}
			stableCount += incr;
		} else {
			if (verbose) {
				std::cerr << beta << '\t' << pri.totalSize << '\t'
						<< pri.avgWeightPatches << '\t';
				std::cerr << pri.avgWeightNoPatch << '\t' << pri.maxDistance
						<< '\t';
				std::cerr << pri.averageDistance << '\t' << pri.numberOfPatches
						<< '\t' << stableCount << '\n';
			}
			// change
			lastBeta = beta;
			stableCount = 0;
		}
		oldWeight = pri.avgWeightPatches;
		beta += incr;
	}

	if (verbose) {
		std::cerr << "Setting final beta estimate to " << beta / 2.0f << "\n";
	}

	return size_t(beta / 2.0f);
}

float AntiPatcher::getBeta() {
	return this->beta;
}

void AntiPatcher::setBeta(float beta) {
	this->beta = beta;
}

float AntiPatcher::getDelta() {
	return this->delta;
}

float AntiPatcher::getKernelDensity(std::vector<float> x) {
	std::map<int, float> xDens, yDens;
	std::sort(x.begin(), x.end(), std::less<float>());

	float dens = 0.0, mean = getMean(x);

	std::vector<float>::iterator it = x.begin();
	// Bandwith = vectorsize / 100
	float h = 100.;

	int i = 0;
	for (it = x.begin(); it != x.end(); ++it) {
		xDens[i] = *it;
		yDens[i] = getKernel((mean - *it) / h);
		dens += yDens[i++];
	}

	float max = 0.0;
	int index = 0;
	for (std::map<int, float>::iterator it = yDens.begin(); it != yDens.end();
			++it)
		if (it->second > max)
			index = it->first;

	index = (int) ((float) index * 0.95);

	return xDens[index];
}

float AntiPatcher::getMean(std::vector<float> vec) {
	float sum = 0.0;

	std::vector<float>::iterator it;
	for (it = vec.begin(); it != vec.end(); ++it)
		sum += *it;

	return sum / ((float) vec.size());
}

float AntiPatcher::getKernel(float f) {
	return 1 / sqrt(2 * M_PI) * exp(-0.5 * f * f);
}

std::vector<float> AntiPatcher::getPairwiseDistances(bool withWeight) {
	std::vector<float> dist;
	Patch::iterator it, jt;

	for (it = this->_residues.begin(); it != this->_residues.end(); ++it)
		for (jt = this->_residues.begin(); jt != it; ++jt) {
			if (withWeight)
				dist.push_back((*it)->getWeightedDistance(*jt));
			else
				dist.push_back((*it)->getDistance(*jt));
		}
	return dist;
}

