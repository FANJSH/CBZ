#ifndef ENDFIDTRANSLATOR
#define ENDFIDTRANSLATOR

#include <vector>
#include <string>
#include <iostream>
#include <stdlib.h>

#include "Numeric.h"

using namespace std;

class ENDFIDTranslator{
 private:
  int nuclnum;
  vector<string> nucln;
  vector<string> nucln_srac;
  vector<int> mass_b;
 public:
  ENDFIDTranslator();
  int GetENDFID(string name);
  string GetNuclideName(int id);
  int ID(string name){return GetENDFID(name);};
  int NewID(int ia,int iz,int il);
  string Name(int id){return GetNuclideName(id);};
  string NewName(int id);
  string Name(int iz,int ia,int il,bool srac=false){return SearchNameFromParameter(iz,ia,il,srac);};
  void GetParameter(string name,int &iz,int &ia,int &il);
  void GetParameter(int id,int &iz,int &ia,int &il){GetParameter(Name(id),iz,ia,il);};
  void GetParameterNew(int id,int &iz,int &ia,int &il);
  string SearchNameFromParameter(int iz,int ia,int il,bool srac=false);
  int SearchMatIDFromParameter(int iz,int ia,int il)
  {return ID(SearchNameFromParameter(iz,ia,il));};
  string NameFromAtomicNumber(int an);
  int AtomicNumberFromName(string name);
  string ExtractIsotopeName(string name);
  void check();
};

#endif
