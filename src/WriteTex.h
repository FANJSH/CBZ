#ifndef WRITETEX
#define WRITETEX

#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include "CartCore.h"

using namespace std;

class WriteTex{
 protected:
  FILE *outp;
 public:
  WriteTex(){outp=fopen("map.tex","w");};
  void Begin(int xwid,int ywid);
  void End();
  void WriteCircle(int x,int y,int r);
  void WriteBCircle(int x,int y,int r);
  void WriteCross(int x,int y,int l);
  void WriteKagis(int x,int y,int l);
  void WriteKagi1(int x,int y,int l);
  void WriteKagi2(int x,int y,int l);
  void WriteKagi3(int x,int y,int l);
  void WriteSquare(int x,int y,int xl,int yl);
  void WriteMinus(int x,int y,int l);
  void WriteOMinus(int x,int y,int l);
  //void WriteFigure3D(CartCore &inp,float unit); //unit:size of assembly
  // void WriteFigure2D(CartCore &inp,float unit,float z);
  void WriteAssembly(int id,int x,int y,int size);
  void WriteMaterial(int id,int x,int y);
  //void CalXY(CartCore &inp,int &x,int &y,float unit);
};

#endif
