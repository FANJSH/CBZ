#ifndef SPHERICALHARMONICS
#define SPHERICALHARMONICS

#include "Numeric.h"

real Kaijo(int m);
real AssociatedLegendre(int l,int m,real mu);
real Legendre(int l,real mu);
real Ylm(int l,int m,real mu);
real Harmonics(int l,int m,real mu,real et,real xi);
void RootOfLegendre(int pl,real *root);
real WeightOfGaussian(int pl, real mu);

real Chebyshev(int l,real x);
real ChebyshevSecondKind(int l,real x);

#endif
