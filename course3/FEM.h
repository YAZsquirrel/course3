#pragma once
#include "Filtration.h"
#include <cmath>
#include <iostream>
#include <functional>
#include <list>
#include "GridMaker.h"
typedef double real;

using namespace filtration;
using namespace mesh_comps;

class FEM
{
   private:

   struct bound {
      int knots_num[4];
      real value;
   };
   std::list<bound*> bounds1;
   
   //struct bound2 {
   //   face* bound_face;
   //   real value; // theta
   //};
   std::list<bound*> bounds2;

   struct Matrix {
      real *l, *u, *d;
      int *ig, *jg;
   };
   hexahedron* hexas;
   //real lambda(real knot[2], real t);
   real f(knot *knot_);
   real ug(knot *knot_);
   //real gamma = 1, theta = 1, lambda = 1;
   int num_of_knots, num_of_FE, un;

   real localM2d[4][4];
   real localM[8][8]; // 8*8
   real localG[8][8];
   real localA[8][8];
   real reversed_J[3][3];
   real J[3][3];
   real Jgrad_i[3];
   real gradi[3];
   real Jgrad_j[3];
   real gradj[3];

   void WriteMatrix(Matrix* A);
   void MakeSparseFormat();
   void AddElement(Matrix *A, int knot_num[8], int i, int j, real elem);
   void AddLocal(Matrix* A, int knot_num[8], real localA[8][8], real coeff);
   void AddFirstBounds();
   void AddSecondBounds(/*hexahedron* hexa*/);
   void CreateA(hexahedron* hexa);
   void CreateSLAE();
   void CreateM(hexahedron *hexa); // можно посылать параметр, на кот. умножается матрица
   void CreateG(hexahedron *hexa); // лямбду надо при формировании
   void Createb(hexahedron* hexa);

   Matrix *A;
   real *b, *qk, *temp, *t;
   real *q; 
   void copy(real *x, real *y);
   real scalar(real* v, real* u, int size);
   void MatxVec(real* v, Matrix* A, real* b);
   real* z, *r, *p, *ff, *x;
   Filtration* filtr;
   Mesh* mesh;
   void SolveSLAE();

   real Integrate(const std::function<real(real, real, real, int, int, int[8])> f, int i, int j, int knot_num[8]);
   real Integrate2D(const std::function<real(real, real, int, int, int[4])> f, int i, int j, int knot_num[4]);
   inline int mu(int index);
   inline int v(int index);
   inline int nu(int index);
   real W(int index, real alpha);
   real d_phi(int index, int what, real ksi, real etta, real tetha);
   real prime_by_var(int what, int varOnFE, int knot_num[8], real ksi, real etta, real tetha);
   inline real phi(int index, real ksi, real etta, real tetha);
   inline real det_J();
   void calc_grad(int ij, int index, real ksi, real etta, real tetha);
   std::function<real(real, real, real, int, int, int[8])> Gij;
   std::function<real(real, real, real, int, int, int[8])> Mij;


   public:
   FEM();
   void SolveElliptic();
   void GetSolutionOnPlane(real z);
   void Output(std::ofstream &out);
};