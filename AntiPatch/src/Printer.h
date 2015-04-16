/* 
 * File:   Printer.h
 * Author: lklesper
 *
 * Created on October 29, 2012, 9:47 AM
 */

#ifndef PRINTER_H
#define	PRINTER_H

#include "Residue.h"
#include <stdio.h>
#include <map>
#include <vector>
#include <string>

/** class Printer
 *
 *	These class is in first line to print pymol skripts. You can also print summarys with it.
 */

class Printer
{
public:
  Printer (std::vector< std::vector<Residue*> > patches, std::vector<Residue*> res, std::string pdbPath, int detail);
  Printer (const Printer& orig);
  virtual ~Printer ();

  void printSurface (std::string expName, std::string pictureName, std::string chain = "", int chainNo = 6 /*A - F*/);

  void printCluster (std::string expName, std::string pictureName, std::string chain = "", int clusterSize = 1);

  void printSummary (std::string expName, int clusterSize, std::vector<int> beta, float uBeta, std::vector<float> dist, float delta, std::vector<int> sortUsedRes, std::vector<Residue*> intermediate);

  void printStringForPythonSkript (std::string expName, int clusterSize);

  void printPValueVSResidueNr (std::string expName, std::vector<Residue*> res);

  void printSimple();

private:
  std::vector<std::string> surfaceColor;
  std::vector<std::string> clusterColor;
  std::vector< std::vector<Residue*> > patches;
  std::vector<Residue*> residues;

  std::string pdbPath;
  int detail;

  std::string getObject (size_t k, std::vector<std::string> color, std::string name);
  std::string getResi (std::string chain, size_t k, std::string name);
  std::string getNameWithoutPath (std::string);
  std::string getFirstTwoPicturePlots (std::string name);
  std::string getRotatedPictureName (std::string name, int n, int rotate);
  std::string getSortedResidueAlias ();
  /**
   * returns a list with all Residues from chain sorted by atoms
   * @param chain
   * @return 
   */
  std::string getSortedResidueValues (std::string chain);

  std::string printResAndB (std::string chain, int chainLength);

  void printPml (std::string expName, std::string pictureName, std::string display, int rotate = 0);
  void sortResidues (std::vector<Residue*> res);
  std::map<int, int> getMapHowOften (std::string chain);

  void printPatchForSummary (std::ofstream &of, int clusterSize);
  void printBetaForSummary (std::ofstream &of, std::vector<int> beta, float uBeta);
  void printDeltaForSummary (std::ofstream &of, std::vector<float> dist, float delta);
  void printResPValForSummary (std::ofstream& of);
  void printSortedUsedResidue (std::ofstream& of, std::vector<int> sortRes);
  void printIntermediates (std::ofstream& of, std::vector<Residue*> intermediate);

  std::string getFQname (std::string name);
};

inline bool
sortVectorLength (std::vector<Residue*> a, std::vector<Residue*> b){ return a.size () > b.size (); }

#endif	/* PRINTER_H */

