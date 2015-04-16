/* 
 * File:   Residue.h
 * Author: l.Klesper
 *
 * Created on 2. August 2012, 12:38
 */

#ifndef RESIDUE_H
#define	RESIDUE_H

#include <string>
#include <vector>
#include "Point.h"

/** class Residue
 *
 *	In this class all Information about a residue will be saved.
 * 	position in data file, 3 d point in structure, pvalue (first maybe not known)
 *	a alias (default = position), epi for label, bool value if the residue is on surface.
 *	also there are setter and getter for all these variables and
 *	functions to calculate distances between residues and 
 */
class Residue
{
private:
  std::string pos;
  size_t count;
  float value;
  std::string label;
  std::string alias;
  std::string epi;
  Point coord;
  float pValue;
  float agw;
  bool patched;
  bool surface;

public:
  Residue (std::string pos, Point c, float pval, std::string alias, std::string epi);

  Residue () { };

  void setPValue (float);
  void setAgW (float);

  Point getCoord ();
  void setCoord (float x, float y, float z);
  float getPValue ();
  float getAgW ();
  float getWeightedDistance (Residue* res);
  float getDistance (Residue* res);
  float getEuklideanDistance (Residue* res);
  std::string getLabel ();
  std::string getAlias ();
  std::string getEpi ();
  std::string getPos ();


  void
  setSurfaceInfo (bool surface) { this->surface = surface; };

  bool
  isOnSurface () { return surface; };

  bool isPatched ();
  void setPatched (bool j);
};

#endif	/* RESIDUE_H */

