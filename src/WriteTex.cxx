#include <iostream>
#include "WriteTex.h"

using namespace std;

void WriteTex::Begin(int xwid,int ywid)
{
  fprintf(outp,"\\documentclass[10pt]{article}\n");
  fprintf(outp,"\\begin{document}\n");
  fprintf(outp,"\\begin{figure}[H]\n");
  fprintf(outp,"\\centering\n");
  fprintf(outp,"\\begin{picture}(%d,%d)\n",xwid,ywid);
};

void WriteTex::End()
{

  fprintf(outp,"\\end{picture}\n");
  fprintf(outp,"\\end{figure}\n");
  fprintf(outp,"\\end{document}\n");
};

void WriteTex::WriteCircle(int x,int y,int r)
{
  fprintf(outp,"\\put(%d,%d){\\circle{%d}}\n",x,y,r);
};

void WriteTex::WriteBCircle(int x,int y,int r)
{
  fprintf(outp,"\\put(%d,%d){\\circle*{%d}}\n",x,y,r);
};

void WriteTex::WriteCross(int x,int y,int l)
{
  fprintf(outp,"\\put(%d,%d){\\line(1,0){%d}}\n",x-l,y,l*2);
  fprintf(outp,"\\put(%d,%d){\\line(0,1){%d}}\n",x,y-l,l*2);
};

void WriteTex::WriteKagis(int x,int y,int l)
{
  fprintf(outp,"\\put(%d,%d){\\line(1,0){%d}}\n",x,y,l);
  fprintf(outp,"\\put(%d,%d){\\line(0,1){%d}}\n",x,y,l);
};
void WriteTex::WriteKagi1(int x,int y,int l)
{
  fprintf(outp,"\\put(%d,%d){\\line(-1,0){%d}}\n",x,y,l);
  fprintf(outp,"\\put(%d,%d){\\line(0,1){%d}}\n",x,y,l);
};
void WriteTex::WriteKagi2(int x,int y,int l)
{
  fprintf(outp,"\\put(%d,%d){\\line(1,0){%d}}\n",x,y,l);
  fprintf(outp,"\\put(%d,%d){\\line(0,-1){%d}}\n",x,y,l);
};
void WriteTex::WriteKagi3(int x,int y,int l)
{
  fprintf(outp,"\\put(%d,%d){\\line(-1,0){%d}}\n",x,y,l);
  fprintf(outp,"\\put(%d,%d){\\line(0,-1){%d}}\n",x,y,l);
};

void WriteTex::WriteMinus(int x,int y,int l)
{
  fprintf(outp,"\\put(%d,%d){\\line(1,0){%d}}\n",x-l,y,l*2);
};

void WriteTex::WriteOMinus(int x,int y,int l)
{
  fprintf(outp,"\\put(%d,%d){\\line(0,1){%d}}\n",x,y-l,l*2);
};

void WriteTex::WriteAssembly(int id,int x,int y,int size)
{
  if(id==1){WriteCircle(x,y,size);};
  if(id==2){WriteBCircle(x,y,size);};
  if(id==3){WriteCross(x,y,size/2);};
  if(id==4){WriteMinus(x,y,size/2);};
  if(id==5){WriteSquare(x-size/2,y-size/2,size,size);};
  if(id==6){WriteOMinus(x,y,size/2);};
  if(id==7){WriteKagis(x,y,size/2);};
};

void WriteTex::WriteMaterial(int mat,int x,int y)
{
  if(mat==1){WriteCircle(x,y,3);};
  if(mat==2){WriteBCircle(x,y,3);};
  if(mat==3){
    WriteCircle(x-2,y,3);
    WriteCircle(x+2,y,3);
  };
  if(mat==4){
    WriteBCircle(x-2,y,3);
    WriteBCircle(x+2,y,3);
  };
  if(mat==5){WriteMinus(x,y,3);};
  if(mat==6){WriteOMinus(x,y,3);};
  if(mat==7){WriteKagis(x,y,3);};
  if(mat==8){WriteSquare(x-2,y-2,4,4);};
  if(mat==9){WriteKagi1(x,y,3);};
  if(mat==10){WriteKagi2(x,y,3);};
  if(mat==11){WriteKagi3(x,y,3);};
};

void WriteTex::WriteSquare(int x,int y,int xl,int yl)
{
  fprintf(outp,"\\put(%d,%d){\\line(1,0){%d}}\n",x,y,xl);
  fprintf(outp,"\\put(%d,%d){\\line(0,1){%d}}\n",x+xl,y,yl);
  fprintf(outp,"\\put(%d,%d){\\line(-1,0){%d}}\n",x+xl,y+yl,xl);
  fprintf(outp,"\\put(%d,%d){\\line(0,-1){%d}}\n",x,y+yl,yl);
};

/*
void WriteTex::CalXY(CartCore &core,int &x,int &y,float unit)
{
  y=0;
  for(int i=0;i<core.GetYr();i++){
    y+=int(12*core.GetYwid(i)/unit);
  };

  x=0;
  for(int i=0;i<core.GetXr();i++){
    x+=int(12*core.GetXwid(i)/unit);
  };
};
*/
/*
void WriteTex::WriteFigure3D(CartCore &core,float unit)
{
  int xx,yy;
  int xr=core.GetXr();
  int yr=core.GetYr();

  CalXY(core,xx,yy,unit);

  int xwid=xx+150;
  if(xwid<300)xwid=300;
  int ywid=yy+260;

  Begin();

  fprintf(outp,"\\begin{picture}(%d,%d)\n",xwid,ywid);
  WriteSquare(0,0,xwid,ywid);

  int x,y;
  int ind=0;
  y=ywid-30;

  fprintf(outp,"\\put(10,%d){(Radial map)}\n",y+10);
  for(int i=0;i<yr;i++){
    int yl=int(12*core.GetYwid(i)/unit);
    y-=yl;
    x=10;
    for(int j=0;j<xr;j++){
      int xl=int(12*core.GetXwid(j)/unit);
      int sml=yl/2;
      if(xl/2<sml)sml=xl/2;
      WriteSquare(x,y,xl,yl);
      WriteAssembly(core.GetAsmMap(ind),x+xl/2,y+yl/2,sml);
      x+=xl;
      ind++;
    };
  };

  y=ywid-25;
  x+=20;
  for(int i=0;i<core.GetAset()->GetAsmNum();i++){
    y-=15;
    WriteSquare(x,y,12,12);
    WriteAssembly(i,x+6,y+6,6);
    fprintf(outp,"\\put(%d,%d){%s}\n",x+15,y+2,core.GetAset()->GetAsm(i).GetName());
  };

  y=ywid-30-yy-15;
  fprintf(outp,"\\put(20,%d){(Boundary condition)}\n",y);
  if(core.GetRightBC()==0){fprintf(outp,"\\put(30,%d){E: Vacuum}\n",y-12);}
  else{fprintf(outp,"\\put(30,%d){E: Reflective}\n",y-12);};
  if(core.GetLeftBC()==0){fprintf(outp,"\\put(100,%d){W: Vacuum}\n ",y-12);}
  else{fprintf(outp,"\\put(100,%d){W: Reflective}\n",y-12);};
  if(core.GetBackBC()==0){fprintf(outp,"\\put(30,%d){N: Vacuum}\n",y-12*2);}
  else{fprintf(outp,"\\put(30,%d){N: Reflective}\n",y-12*2);};
  if(core.GetFrontBC()==0){fprintf(outp,"\\put(100,%d){S: Vacuum}\n",y-12*2);}
  else{fprintf(outp,"\\put(100,%d){S: Reflective}\n",y-12*2);};

  x=10;
  y-=50;
  fprintf(outp,"\\put(10,%d){(Axial map)}\n",y+10);
  y-=10;
  float w;

  for(int i=0;i<core.GetAset()->GetAsmNum();i++){
    WriteSquare(x,y,12,12);
    WriteAssembly(i,x+6,y+6,6);
    WriteSquare(x,y-105,12,100);
    int ww=0;
    int wwold;

    int zmax=core.GetAset()->GetAsm(i).GetZdiv();
    for(int j=0;j<zmax;j++){
      w=core.GetAset()->GetAsm(i).GetZdist(zmax-1);
      wwold=ww;
      ww=int(100*core.GetAset()->GetAsm(i).GetZdist(j)/w);
      if(j!=zmax-1)fprintf(outp,"\\put(%d,%d){\\line(1,0){12}}\n",x,y-5-ww);
      int mat=core.GetAset()->GetAsm(i).GetZmap(j);
      int za=(ww-wwold)/2+wwold;
      WriteMaterial(mat,x+6,y-5-za);
    };
    x+=15;
  };

  x=10+12*xr+10;
  //  y-=10;
  for(int i=0;i<core.GetAset()->GetMatNum();i++){
    y-=15;
    WriteSquare(x,y,12,12);
    WriteMaterial(i,x+6,y+6);
    fprintf(outp,"\\put(%d,%d){%s}\n",x+15,y+2,core.GetAset()->GetMaterialName(i));
  };

  y=ywid-yy-220;
  fprintf(outp,"\\put(20,%d){(Boundary condition)}\n",y);
  if(core.GetUpperBC()==0){fprintf(outp,"\\put(40,%d){Top: Vacuum} ",y-12);}
  else{fprintf(outp,"\\put(40,%d){Top: Reflective} ",y-12);};
  if(core.GetBottomBC()==0){fprintf(outp,"\\put(120,%d){Bottom: Vacuum} ",y-12);}
  else{fprintf(outp,"\\put(120,%d){Bottom: Reflective} ",y-12);};

  fprintf(outp,"\\put(5,5){Size of Assembly : %7.4f*%7.4f*%8.3f(cm)}",unit,unit,w);

  End();
};
*/

/*

void WriteTex::WriteFigure2D(CartCore &core,float unit,float z)
{
  Begin();

  int xx,yy;
  int xr=core.GetXr();
  int yr=core.GetYr();

  CalXY(core,xx,yy,unit);

  int xwid=xx+150;
  if(xwid<300)xwid=300;
  int ywid=yy+100;

  fprintf(outp,"\\begin{picture}(%d,%d)\n",xwid,ywid);
  WriteSquare(0,0,xwid,ywid);

  int x,y;
  int ind=0;
  y=ywid-20;

  int xymap[xr*yr];
  core.GetMaterialMapXY(z,xymap);

  fprintf(outp,"\\put(10,%d){(Radial map)}\n",y+10);
  for(int i=0;i<yr;i++){
    int yl=int(12*core.GetYwid(i)/unit);
    y-=yl;
    x=10;
    for(int j=0;j<xr;j++){
      int xl=int(12*core.GetXwid(j)/unit);
      int sml=yl/2;
      if(xl/2<sml)sml=xl/2;
      WriteSquare(x,y,xl,yl);
      WriteAssembly(xymap[ind],x+xl/2,y+yl/2,sml);
      x+=xl;
      ind++;
    };
  };

  y=ywid-30;
  x=10+12*xr+30;
  for(int i=0;i<core.GetAset()->GetMatNum();i++){
    bool exist=false;
    for(int j=0;j<xr*yr;j++){
      if(xymap[j]==i)exist=true;
    };
    if(exist){
      y-=15;
      WriteSquare(x,y,12,12);
      WriteAssembly(i,x+6,y+6,6);
      fprintf(outp,"\\put(%d,%d){%s}\n",x+15,y+2,core.GetAset()->GetMaterialName(i));
    };
  };

  y=ywid-30-yy;
  fprintf(outp,"\\put(20,%d){(Boundary condition)}\n",y);
  if(core.GetRightBC()==0){fprintf(outp,"\\put(40,%d){E: Vacuum}\n",y-12);}
  else{fprintf(outp,"\\put(40,%d){E: Reflective}\n",y-12);};
  if(core.GetLeftBC()==0){fprintf(outp,"\\put(40,%d){W: Vacuum}\n",y-12*2);}
  else{fprintf(outp,"\\put(40,%d){W: Reflective}\n",y-12*2);};
  if(core.GetBackBC()==0){fprintf(outp,"\\put(40,%d){N: Vacuum}\n",y-12*3);}
  else{fprintf(outp,"\\put(40,%d){N: Reflective}\n",y-12*3);};
  if(core.GetFrontBC()==0){fprintf(outp,"\\put(40,%d){S: Vacuum}\n",y-12*4);}
  else{fprintf(outp,"\\put(40,%d){S: Reflective}\n",y-12*4);};

  fprintf(outp,"\\put(30,5){Size of Assembly : %7.4f*%7.4f(cm)}",unit,unit);

  End();
};
*/
