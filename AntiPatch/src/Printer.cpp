/* 
 * File:   Printer.cpp
 * Author: lklesper
 * 
 * Created on October 29, 2012, 9:47 AM
 */

#include "Printer.h"
#include "AntiPatch.h"
#include "Exception.h"
#include <iomanip>
#include <stdlib.h>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>

#ifdef WINDOWS
#include <direct.h>
#define GETDIR _getcwd
#define PATHSEP "\"
#else
#include <unistd.h>
#define GETDIR getcwd
#define PATHSEP "/"
#endif

Printer::Printer (std::vector< std::vector<Residue*> > patches, std::vector<Residue*> res, std::string pdb, int detail)
: patches (patches), detail (detail)
{
  // get fully qualified Path Name
  pdbPath = getFQname (pdb);
  // sort Residues 
  sortResidues (res);
  // TODO: use boost or config File to do this:
  surfaceColor.push_back ("grey");
  surfaceColor.push_back ("grey");
  surfaceColor.push_back ("white");
  surfaceColor.push_back ("white");
  surfaceColor.push_back ("white");
  surfaceColor.push_back ("white");
  clusterColor.push_back ("green");
  clusterColor.push_back ("blue");
  clusterColor.push_back ("red");
  clusterColor.push_back ("yellow");
  clusterColor.push_back ("purple");
  clusterColor.push_back ("orange");
  clusterColor.push_back ("cyan");
  clusterColor.push_back ("dgrey");
  clusterColor.push_back ("salmon");
  clusterColor.push_back ("pink");
  clusterColor.push_back ("dblue");
  clusterColor.push_back ("lgreen");
  clusterColor.push_back ("magenta");
}

Printer::Printer (const Printer& orig) { }

Printer::~Printer () { }

void
Printer::printSurface (std::string expName, std::string pictureName, std::string chain, int chainLength)
{
  expName = getFQname (expName).append ("_surface");
  printPml (expName, pictureName.append ("_surface"), "surface", 3);
  std::ofstream of (expName.append (".pym").c_str (), std::ios::out);
  of << "select sp_" << chain << ", chain " << chain << std::endl;
  of << "color grey" << std::endl;
  of << printResAndB (chain, chainLength);

}

std::string
Printer::printResAndB (std::string chain, int chainLength)
{
  int n = this->residues.size () / 100 + ((this->residues.size () % 100 == 0) ? 0 : 1), N, nr;
  std::stringstream of, alias, res, currentChain;
  std::map<int, int> howOften = getMapHowOften (chain);

  for (int i = 0; i < n; ++i)
    {
      currentChain << "sp_" << chain << "r" << i;
      of << "select " << currentChain.str () << ", chain " << chain << " and resi ";
      N = i * 100 + 100;
      if (N > this->residues.size ()) N = this->residues.size ();
      for (int j = i * 100; j < N; ++j)
        {
          alias << residues.at (j)->getAlias () << (j < N - 1 ? "+" : "");
          nr = howOften[atoi (residues.at (j)->getAlias ().c_str ())];
          for (int k = 0; k < nr; ++k)
            res << residues.at (j)->getPValue () << (j < N - 1 ? "," : "");
        }
      of << alias.str () << std::endl << "di = iter ([" << res.str () << "])" << std::endl;
      of << "alter " << currentChain.str () << ", b = 0" << std::endl;
      of << "alter " << currentChain.str () << ", b = di.next()" << std::endl;
      of << "spectrum b, red_blue, " << currentChain.str () << ", minimum=0, maximum=1" << std::endl;

      alias.str ("");
      res.str ("");
      currentChain.str ("");
    }

  return of.str ();
}

void
Printer::printCluster (std::string expName, std::string pictureName, std::string chain, int clusterSize)
{
  // extends Path from name:
  std::string name = getNameWithoutPath (expName);
  // open stream
  expName = getFQname (expName).append ("_cluster");
  printPml (expName, pictureName, "cartoon", 0);
  std::ofstream of (expName.append (".pym").c_str (), std::ios::out);

  size_t k = 1; // k objects
  std::vector< std::vector<Residue*> >::iterator it;

  for (it = this->patches.begin (); it != this->patches.end (); it++)
    if (it->size () >= clusterSize)
      {
        of << getResi (chain, k, name);

        for (size_t j = 0; j < it->size (); j++)
          {
            of << (*it)[j]->getAlias ();
            if (j != it->size () - 1) of << "+";
          }

        k++;
        of << getObject (k, clusterColor, name);
      }

  of << std::endl;
  of.close ();
}

void
Printer::printSummary (std::string expName, int clusterSize, std::vector<int> beta, float uBeta, std::vector<float> dist, float delta, std::vector<int> sortUsedRes, std::vector<Residue*> intermediate)
{
  std::string path = expName.append (".summary");
  std::ofstream of (path.c_str (), std::ios::out);

  // print Patches:
  of << "-----------------------------------------------" << std::endl;
  printPatchForSummary (of, clusterSize);
  of << "-----------------------------------------------" << std::endl;
  printBetaForSummary (of, beta, uBeta);
  of << "-----------------------------------------------" << std::endl;
  printDeltaForSummary (of, dist, delta);
  of << "-----------------------------------------------" << std::endl;
  printResPValForSummary (of);
  of << "-----------------------------------------------" << std::endl;
  printSortedUsedResidue (of, sortUsedRes);
  of << "-----------------------------------------------" << std::endl;
  printIntermediates (of, intermediate);

  of.close ();
}

#include <algorithm>    // std::sort
void Printer::printSimple() {
	  std::cout << "------------ " << "patches: " << " ------------" << std::endl;
	  std::vector< std::vector<Residue*> >::iterator iter;
	  std::vector<Residue*> patch;
	  std::vector<Residue*>::iterator resIt;
	  size_t k=1;
	  for (iter=patches.begin(); iter!=patches.end(); iter++) {
		  patch = *iter;
		  if (patch.size()<2) continue;
		  std::cout << k << ".\t";
		  std::sort (patch.begin(), patch.end());
		  for (resIt=patch.begin(); resIt!=patch.end(); resIt++) {
			  std::cout << (*resIt)->getAlias () << "\t";
		  }
		  std::cout << "\n";
		  k++;
	  }



}


void
Printer::printIntermediates (std::ofstream& of, std::vector<Residue*> intermediate)
{
  of << "------------ " << "intermediate Residues: " << " ------------" << std::endl;
  for (Patch::iterator it = intermediate.begin (); it != intermediate.end (); ++it)
    of << (*it)->getAlias () << "\t";
  of << std::endl;
}

void
Printer::printSortedUsedResidue (std::ofstream& of, std::vector<int> sortRes)
{
  of << "------------ " << "used Residues: " << " ------------" << std::endl;
  for (std::vector<int>::iterator it = sortRes.begin (); it != sortRes.end (); ++it)
    of << *it << "\t";
  of << std::endl;
}

void
Printer::printResPValForSummary (std::ofstream& of)
{
  of << "------------ " << "Residues in Chain wit p-Value/normalized antigenic weight: " << " ------------" << std::endl;
  of << "sorted Residues   p-Value/normailzed antigenic weight: " << std::endl << "------------" << std::endl;
  for (std::vector<Residue*>::iterator it = residues.begin (); it != residues.end (); ++it)
    of << std::setw (8) << std::left << (*it)->getPValue () << std::setw (8) << (*it)->getAlias () << std::endl;
}

void
Printer::printDeltaForSummary (std::ofstream& of, std::vector<float> dist, float delta)
{
  of << "------------ " << "Deltasearch: " << " ------------" << std::endl;
  of << "used delta -> " << delta << std::endl;
  of << "sorted pairwise distances: " << std::endl << "------------" << std::endl;
  int count = 0;
  for (std::vector<float>::iterator it = dist.begin (); it != dist.end (); ++it)
    {
      of << *it << "\t";
      if (count++ % 10 == 0 && count != 1)
        of << std::endl;
    }
  of << std::endl;
}

void
Printer::printBetaForSummary (std::ofstream &of, std::vector<int> beta, float uBeta)
{
  of << "------------ " << "Betasearch: " << " ------------" << std::endl;
  of << "used Beta -> " << uBeta << std::endl;
  of << "vector of patchsize for different beta: "
          << std::endl << "------------" << std::endl;
  int count = 0;
  for (std::vector<int>::iterator it = beta.begin (); it != beta.end (); ++it)
    {
      of << *it << "\t";
      if (count++ % 10 == 0 && count != 1)
        of << std::endl;
    }
  of << std::endl;
}

void
Printer::printPatchForSummary (std::ofstream &of, int clusterSize)
{
  // comments are if you want to se mean values
  std::vector< std::vector<Residue*> >::iterator it;
  std::vector<int> res;
  std::sort (patches.begin (), patches.end (), sortVectorLength);
  of << "------------ " << "found Patches: " << " ------------" << std::endl;
  of << "patchnr\t[residues]\tcounted residues" << std::endl << "------------" << std::endl;
  int count = 1;
  //float avg, agwAvg;

  for (it = this->patches.begin (); it != this->patches.end (); it++)
    {
      if (it->size () >= clusterSize)
        {
          //avg = 0.0;
          //agwAvg = 0.0;
          of << count++ << "\t[";
          for (size_t j = 0; j < it->size (); j++)
            {
              //avg += (*it)[j]->getPValue ();
              //agwAvg += (*it)[j]->getAgW ();
              res.push_back (atoi ((*it)[j]->getAlias ().c_str ()));
            }
          std::sort (res.begin (), res.end ());
          for (size_t j = 0; j < it->size (); j++)
            {
              of << res.at (j);
              if (j != it->size () - 1) of << ", ";
              else
                {
                  of << "]\t" << res.size();
                  /*of << "\t(" << avg / it->size ();
                  if (agwAvg > 0.0) of << ", " << agwAvg / it->size ();
                  of << ")";*/
                }
              
            }
          of << std::endl;
          res.clear ();
        }
    }
}

void
Printer::printPml (std::string expName, std::string pictureName, std::string display, int rotate)
{
  std::string pml = expName, pym = expName;
  std::ofstream of (pml.append (".pml").c_str (), std::ios::out);

  of << "load " << pdbPath << std::endl << "hide everything, all" << std::endl
          << "show " << display << ", all" << std::endl << "bg_color white" << std::endl
          << "color white" << std::endl << "@" << pym.append (".pym") << std::endl
          << "reset" << std::endl;
  of << getFirstTwoPicturePlots (pictureName);
  if (rotate == 0)
    of << getRotatedPictureName (pictureName, 0, 0);
  for (int i = 0; i < rotate; ++i)
    of << getRotatedPictureName (pictureName, i, 360 / rotate);
  of << "save " << pictureName.append (".psw") << ", format=pse";
}

std::string
Printer::getFirstTwoPicturePlots (std::string name)
{
  std::string pdb = getNameWithoutPath (pdbPath);
  pdb = pdb.substr (0, pdb.find ('.'));
  std::string left = name, right = name, top = name;
  std::stringstream res;
  res << "turn z, -35" << std::endl
          << "turn x, 20" << std::endl
          << "turn y, -20" << std::endl;
  res << "translate [0,5,-15], object=" << pdb << std::endl;
  res << "ray " << detail << ", " << detail << std::endl;
  res << "png " << right.append ("_right.png") << ", dpi=300" << std::endl;

  res << "turn y, 70" << std::endl;
  res << "ray " << detail << ", " << detail << std::endl;
  res << "png " << left.append ("_left.png") << ", dpi=300" << std::endl;

  res << "turn y, 140" << std::endl;
  res << "turn x, 90" << std::endl;
  res << "ray " << detail << ", " << detail << std::endl;
  res << "png " << top.append ("_top.png") << ", dpi=300" << std::endl;
  res << "turn x, -90" << std::endl;

  return res.str ();
}

std::string
Printer::getRotatedPictureName (std::string name, int n, int rotate)
{
  std::stringstream end, res;

  if (n != 0) end << "_" << n << ".png";
  else end << ".png";

  if (rotate != 0 and n != 0)
    res << "turn y, " << rotate << std::endl;
  res << "ray " << detail << ", " << detail << std::endl;
  res << "png " << name.append (end.str ()) << ", dpi=300" << std::endl;

  return res.str ();
}

std::string
Printer::getNameWithoutPath (std::string expName)
{
  // extends Path from name:
  size_t found = expName.find_last_of ("/\\");
  std::string name;
  if (found != std::string::npos)
    name = expName.substr (found + 1);
  else
    name = expName;

  return name;
}

std::string
Printer::getObject (size_t k, std::vector<std::string> color, std::string name)
{
  std::stringstream res;
  res << std::endl << "create obj" << (k - 1) << ", sp" << (k - 1)
          << "_" << name << std::endl << "show spheres, sp" << (k - 1)
          << "_" << name << std::endl << "color "
          << color.at ((k - 2) % color.size ()) << ", sp"
          << (k - 1) << "_" << name << std::endl << "show surface, "
          << "obj" << (k - 1) << std::endl << "set transparency=0.75"
          << std::endl;
  return res.str ();
}

std::string
Printer::getResi (std::string chain, size_t k, std::string name)
{
  std::stringstream res, ch;
  ch << "chain " << chain << " and ";
  res << "select sp" << k << "_" << name << ", "
          << ((chain != "" & chain != "X") ? ch.str () : "") << "resi ";
  return res.str ();
}

std::string
Printer::getSortedResidueAlias ()
{
  std::stringstream res;

  for (int i = 0; i < residues.size (); ++i)
    res << residues.at (i)->getAlias () << (i < residues.size () - 1 ? "+" : "");
  res << std::endl;

  return res.str ();
}

std::string
Printer::getSortedResidueValues (std::string chain)
{
  std::map<int, int> howOften = getMapHowOften (chain);

  std::stringstream res;
  int n;

  for (int i = 0; i < residues.size (); ++i)
    {
      n = howOften[atoi (residues.at (i)->getAlias ().c_str ())];
      for (int j = 0; j < n; ++j)
        res << residues.at (i)->getPValue () << (i < residues.size () - 1 ? "," : "");
    }

  return res.str ();
}

std::map<int, int>
Printer::getMapHowOften (std::string chain)
{
  std::ifstream ifs (pdbPath.c_str (), std::ios::in);
  if (!ifs)
    throw Exception (pdbPath.c_str (), "could not open");

  std::string tmp;
  std::map<int, int> howOften;
  int current = 0, last = -1;

  while (getline (ifs, tmp, '\n'))
    {
      std::istringstream iss (tmp);
      iss >> tmp;
      if (tmp.compare ("ATOM") == 0)
        {
          for (int i = 0; i < 4; ++i)
            iss >> tmp;

          if (tmp.compare (chain) == 0)
            {
              iss >> current;

              if (current == last)
                howOften[current] = howOften[current] + 1;
              else
                {
                  last = current;
                  howOften[current] = 1;
                }
            }
        }
    }

  return howOften;
}

void
Printer::sortResidues (std::vector<Residue*> res)
{
  std::vector<int> vec;
  std::vector<int>::iterator itv;
  std::vector<Residue*>::iterator it;

  for (size_t i = 0; i < res.size (); ++i)
    vec.push_back (atoi (res[i]->getAlias ().c_str ()));

  std::sort (vec.begin (), vec.end ());

  for (itv = vec.begin (); itv != vec.end (); ++itv)
    for (it = res.begin (); it != res.end (); ++it)
      if (atoi ((*it)->getAlias ().c_str ()) - *itv == 0)
        {
          residues.push_back (*it);
          break;
        }
}

void
Printer::printStringForPythonSkript (std::string expName, int clusterSize)
{
  std::string path = expName.append (".pythonList");
  std::ofstream of (path.c_str (), std::ios::out);

  std::vector<int> res;
  std::sort (patches.begin (), patches.end (), sortVectorLength);
  int currentPatchSize;

  int n = this->patches.size (); // letztes Patch das noch ueber der Patchsize liegt (Danach kein '-' mehr
  for (int i = 0; i < this->patches.size (); ++i)
    if (this->patches.at (i).size () < clusterSize)
      {
        n = i;
        break;
      }

  for (int i = 0; i < n; ++i)
    {
      currentPatchSize = this->patches.at (i).size ();
      if (currentPatchSize >= clusterSize)
        {
          for (size_t j = 0; j < currentPatchSize; j++)
              res.push_back (atoi (this->patches.at (i).at (j)->getAlias ().c_str ()));
            
          std::sort (res.begin (), res.end ());
          for (size_t j = 0; j < currentPatchSize; j++)
              of << res.at (j) << ((j != currentPatchSize - 1)? ",":"");
            
          res.clear ();
        }
      of << ((i == n - 1) ? "" : "-");
    }
  of.close ();
}

void
Printer::printPValueVSResidueNr (std::string expName, Patch res)
{
  std::string path = expName.append (".pValueRes");
  std::ofstream of (path.c_str (), std::ios::out);

  of << "\tpValue\tResidue" << std::endl;
  for (Patch::iterator it = res.begin (); it != res.end (); ++it)
    of << (*it)->getPValue () << "\t" << (*it)->getAlias () << std::endl;

}

std::string
Printer::getFQname (std::string name)
{
  char path[FILENAME_MAX];
  GETDIR (path, sizeof (path));
  std::stringstream ss;
  ss << path << PATHSEP << name;
  return ss.str ();
}

