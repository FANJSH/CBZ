#ifndef MATIDTRANSLATOR
#define MATIDTRANSLATOR

#include <vector>
#include <string>
#include <iostream>
#include <stdlib.h>

#include "Numeric.h"

using namespace std;

class MATIDTranslator{
 private:
  int nuclnum;
  vector<string> nucln;
  vector<int> mass_b;
 public:
  MATIDTranslator();
  string ExtractIsotopeName(string name);
  string NameFromAtomicNumber(int an);
  int AtomicNumberFromName(string name);
  void GetParameter(string name,int &iz,int &ia,int &il);
  void GetParameter(int id,int &iz,int &ia,int &il);
  int GetMATID(int iz,int ia,int il);
  int GetMATID(string name);
  string GetName(int id);
  string GetName(int iz,int ia,int il);
  int ID(string name){return GetMATID(name);};
  int ID(int iz,int ia,int il){return GetMATID(iz,ia,il);};
  string Name(int id){return GetName(id);};
  string Name(int iz,int ia,int il){return GetName(iz,ia,il);};
  //
  int GetMATIDFromENDFID(int mat);
};

#endif
