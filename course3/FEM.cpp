#include "FEM.h"
#include <fstream>
using namespace mesh_comps;

FEM::FEM()
{
#pragma region input
   filtr = new Filtration();
   filtr->Start();
   mesh = filtr->mesh;

   std::ifstream fknots("Knots.txt");
   std::ifstream fhexas("Hexahedrons.txt");
   std::ifstream fbounds1("FirstBounds.txt");
   std::ifstream fbounds2("SecondBounds.txt");
   std::ifstream fparams("EnvParams.txt");
   fparams >> un;
   fparams.close();

   mesh->hexas.clear();
   mesh->knots.clear();

   fknots >> num_of_knots;
   mesh->knots.reserve(num_of_knots);
   
   for (int i = 0; i < num_of_knots; i++)
   {
      knot* k = new knot();
      fknots >> k->x >> k->y >> k->z;
      mesh->knots.push_back(k);
   }
   fknots.close();

   hexahedron* hexa;
   fhexas >> num_of_FE;
   mesh->hexas.reserve(num_of_FE);
   int num;
   for (int i = 0; i < num_of_FE; i++)
   {
      mesh->hexas.push_back(hexa = new hexahedron());
      for (int k = 0; k < 8; k++)
      {
         fhexas >> num;
         hexa->knots_num[k] = num;
      }
      for (int k = 0; k < 8; k++)
      {

         //real lam;
         int minp = -1;
         real ph = 1e10;
         for (int p = 0; p < filtr->phases.size(); p++)
         {
            if (filtr->phases[p].h >= mesh->knots[hexa->knots_num[4]]->z)
            {
               ph = std::min(ph, filtr->phases[p].h);
               if (ph == filtr->phases[p].h || minp != p)
                  minp = p;
            }
         }
         hexa->lam = filtr->por.K * filtr->phases[minp].penetrability / filtr->phases[minp].viscosity;
      }
   }
   fhexas.close();

   int numOfBounds;
   fbounds1 >> numOfBounds;
   for (int i = 0; i < numOfBounds; i++)
   {
      bound* cond = new bound;
      //   fbounds1 >> cond->bound_param;
      for (int j = 0, number; j < 4; j++)
      {
         fbounds1 >> number;								
         cond->knots_num[j] = number;				
         //cond->knots[i] = &knots[number - 1];	
      }
      fbounds1 >> cond->value;
      bounds1.push_back(cond);
   }
   fbounds1.close();

   fbounds2 >> numOfBounds;
   for (int i = 0; i < numOfBounds; i++)
   {
      bound* cond = new bound;

      for (int j = 0, number; j < 4; j++)
      {
         fbounds2 >> number;
         cond->knots_num[j] = number;
      }
      fbounds2 >> cond->value;
      bounds2.push_back(cond);
   }
   fbounds2.close();

#pragma endregion

   MakeSparseFormat();
   q = new real[num_of_knots]{};
   b = new real[num_of_knots]{};
   temp = new real[num_of_knots]{};
   ff = new real[num_of_knots]{};
   z = new real[num_of_knots]{};
   r = new real[num_of_knots]{};
   p = new real[num_of_knots]{};

   Mij = [this](real ksi, real etta, real theta, int i, int j, int knot_num[8])
   {
      for (int i = 0; i < 3; i++)                                  
         for (int j = 0; j < 3; j++)                               
            J[i][j] = prime_by_var(i, j, knot_num, ksi, etta, theta);

      return phi(i, ksi, etta, theta) * phi(j, ksi, etta, theta) * det_J();
   };

   Gij = [this](real ksi, real etta, real theta, int i, int j, int knot_num[8])
   {
      for (int ip = 0; ip < 3; ip++)
         Jgrad_i[ip] = Jgrad_j[ip] = 0.;
      
      for (int ip = 0; ip < 3; ip++)                                       // | dx/d(ksi)   | dy/d(ksi)   | dz/d(ksi)   |
         for (int jp = 0; jp < 3; jp++)                                    // | dx/d(etta)  | dy/d(etta)  | dz/d(etta)  |
            J[ip][jp] = prime_by_var(ip, jp, knot_num, ksi, etta, theta);  // | dx/d(theta) | dy/d(theta) | dz/d(theta) |
      
      // J^-1
      for (int ip = 0; ip < 3; ip++)
         for (int jp = 0; jp < 3; jp++)
         {
            real min[4]{};
            int k = 0;
            for (int im = 0; im < 3; im++)
               for (int jm = 0; jm < 3; jm++)
               {
                  if (im != ip && jm != jp)
                     min[k++] = J[im][jm];
               }
            reversed_J[jp][ip] = ((ip + jp + 2) % 2 ? -1 : 1) * (min[0] * min[3] - min[1] * min[2]);
         }
      
      // grad(phi(ksi, etta, theta))
      calc_grad(1, i, ksi, etta, theta);
      calc_grad(2, j, ksi, etta, theta); 

      // J^-1 * grad(phi)
      for (int ip = 0; ip < 3; ip++)
         for (int jp = 0; jp < 3; jp++)
         {
            Jgrad_i[ip] += reversed_J[ip][jp] * gradi[jp];
            Jgrad_j[ip] += reversed_J[ip][jp] * gradj[jp];
         }
      //TODO: diffuse component (lambda) nado dobavit' ???
      // 
      // 
 
      // Jgrad_i^T * Jgrad_j
      real res = 0;
      for (int ip = 0; ip < 3; ip++)
         res += Jgrad_i[ip] * Jgrad_j[ip];
     
      return res / det_J();
   };

}

void FEM::MakeSparseFormat()
{
   const int N = 8;
   int* list1, * list2;
   int* listbeg = new int[num_of_knots];

   for (int i = 0; i < num_of_knots; i++)
      listbeg[i] = -1;

   list1 = new int[num_of_knots * num_of_knots]{};
   list2 = new int[num_of_knots * num_of_knots]{};
   int listsize = -1, iaddr, ind1, ind2, k;

   for (int iel = 0; iel < num_of_FE; iel++) // ï
   {
      for (int i = 0; i < N; i++) // 
      {
         k = mesh->hexas[iel]->knots_num[i]; //
         for (int j = i + 1; j < N; j++) // need to set N = ?
         {
            ind1 = k;
            ind2 = mesh->hexas[iel]->knots_num[j];  //
            if (ind2 < ind1) //
            {
               ind1 = ind2;
               ind2 = k;
            }
            iaddr = listbeg[ind2];
            if (iaddr == -1) // 
            {
               listsize++;
               listbeg[ind2] = listsize;
               list1[listsize] = ind1;
               list2[listsize] = -1;
            }
            else // 
            {
               while (list1[iaddr] < ind1 && list2[iaddr] >= 0)
                  iaddr = list2[iaddr];
               if (list1[iaddr] > ind1)  // 
               {                         // 
                  listsize++;
                  list1[listsize] = list1[iaddr];
                  list2[listsize] = list2[iaddr];
                  list1[iaddr] = ind1;
                  list2[iaddr] = listsize;
               }
               else if (list1[iaddr] < ind1) // 
               {
                  listsize++;
                  list2[iaddr] = listsize;
                  list1[listsize] = ind1;
                  list2[listsize] = -1;
               }
            }
         }
      }
   }

   A = new Matrix;
   A->ig = new int[num_of_knots + 1]{};
   A->jg = new int[listsize + 1]{};  // +1???

   for (int i = 0; i < num_of_knots; i++)
   {
      A->ig[i + 1] = A->ig[i];

      for (iaddr = listbeg[i]; iaddr != -1; )
      {
         A->jg[A->ig[i + 1]] = list1[iaddr];
         A->ig[i + 1]++;
         iaddr = list2[iaddr];
      }
   }

   //for (int i = 0; i < num_of_knots + 1; i++)
   //   std::cout << A->ig[i] << " ";
   //std::cout << '\n';
   //for (int i = 0; i < listsize + 1; i++)
   //   std::cout << A->jg[i] << " ";
   //std::cout << '\n';

   delete[] listbeg;
   delete[] list1;
   delete[] list2;

   A->l = new real[listsize + 1]{};
   A->u = new real[listsize + 1]{};
   A->d = new real[num_of_knots]{};

}

void FEM::AddElement(Matrix* A, int knot_num[8], int i, int j, real elem)
{
   bool found = false;
   if (i == j)
      A->d[i] += elem;
   else if (i < j)
   {
      int m;
      for (m = A->ig[j]; m < A->ig[j + 1]; m++)
         if (A->jg[m] == i) { found = true; break; }
      if (found)
         A->u[m] += elem; // i-1?
   }
   else
   {  
      int n;
      for (n = A->ig[i]; n < A->ig[i + 1]; n++)
         if (A->jg[n] == j) { found = true; break; }
      if (found)
         A->l[n] += elem; // i-1??
   }

}

void FEM::AddLocal(Matrix* A, int knot_num[8], real localA[8][8], real coeff)
{
   int ibeg, iend, ind;
   for (int i = 0; i < 8; i++)
      A->d[knot_num[i]] += localA[i][i];
   for (int i = 0; i < 8; i++)
   {
      ibeg = A->ig[knot_num[i]];
      for (int j = 0; j < i; j++) // i - 1?
      {
         iend = A->ig[knot_num[i] + 1];  // -1 ?
         while (A->jg[ibeg] != knot_num[j])
         {
            ind = (ibeg + iend) / 2;
            if (A->jg[ind] <= knot_num[j])
               ibeg = ind;
            else
               iend = ind;
         }
         A->l[ibeg] += localA[i][j];
         A->u[ibeg] += localA[j][i];
         ibeg++;
      }
   }

}

void FEM::copy(real* from, real* to)
{
   for (int i = 0; i < num_of_knots; i++)
      to[i] = from[i];
}

void FEM::SolveElliptic()
{
   std::ofstream out("Result.txt");
   CreateSLAE();

   SolveSLAE();
   Output(out);
   out.close();

   //real* check = new real[num_of_knots];
   ////real* check1 = new real[num_of_knots];
   ////for (int i = 0; i < num_of_knots; i++)
   ////   check1[i] = 1.0;
   //MatxVec(check, A, q);
   //for (int i = 0; i < num_of_knots; i++)
   //   std::cout << b[i] << " ";
   //std::cout << '\n';
}

void FEM::GetSolutionOnPlane(real z)
{
   std::ofstream zout("ResultZ.txt");

   for (int i = 0; i < num_of_knots; i++)
   {
      if (abs(mesh->knots[i]->z - z) < 1e-10)
         zout << mesh->knots[i]->x << " " 
               << mesh->knots[i]->y << " "
               //<< mesh->knots[i]->z << " "
               << q[i] / 101325.0 << '\n';
   }

}

void FEM::Output(std::ofstream& out)
{
   //out.scientific;
   //out.precision(15);
   //std::cout.scientific;
   //std::cout.precision(15);
   out.setf(std::ios::right);
   out.width(15);
   out << "\n| x" << std::fixed;
   out.width(15);
   out << "| y";
   out.width(15);
   out << "| z";
   out.width(15);
   out << "| q";
   out.width(15);
   out << "|\n";
   //std::cout << title;

   for (int i = 0; i < num_of_knots; i++)
   {
      out << "|" << mesh->knots[i]->x;
      out.width(15);
      out << "|" << mesh->knots[i]->y;
      out.width(15);
      out << "|" << mesh->knots[i]->z;
      out.width(15);
      out << "|" << q[i];
      out.width(15);
      out << "|\n";
      //out  << "| " << "\t| " << mesh->knots[i]->x << "\t| " << mesh->knots[i]->y << mesh->knots[i]->z << "\t| "
      //   << std::scientific << ug(&mesh->knots[i]) << "\t| " << q[i] << "\t| "
      //   << abs(q[i] - ug(&mesh->knots[i])) << "\t|\n";
      // std::cout << std::scientific << "| " << ug(&knots[i]) << "\t| " << q[i] << "\t| "
      //    << abs(q[i] - ug(&knots[i])) << "\t|\n";
   }

}

void FEM::AddFirstBounds()
{
   for (auto cond : bounds1)
   {
      for (int i = 0; i < 4; i++)
      {
         A->d[cond->knots_num[i]] = 1;
         for (int j = A->ig[cond->knots_num[i]]; j < A->ig[cond->knots_num[i] + 1]; j++)
            A->l[j] = 0.;
         for (int j = 0; j < A->ig[num_of_knots]; j++)
            if (A->jg[j] == cond->knots_num[i])
               A->u[j] = 0.;
      
         b[cond->knots_num[i]] = cond->value;// = ug(mesh->knots[cond->knots_num[i]]);  /// надо будет поменять, наверно, для неизвестных функций из таблицы/по функции
      }
   }
}

void FEM::AddSecondBounds()
{  
   for (auto bound : bounds2)
   {
      
      //real x1 = 0.0, x2 = 0.0, x3 = 0.0, x4 = 0.0,
      //     y1 = 0.0, y2 = 0.0, y3 = 0.0, y4 = 0.0;
      //real a0 = (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1), 
      //     a1 = (x2 - x1) * (y4 - y3) - (y2 - y1) * (x4 - x3), 
      //     a2 = (x3 - x1) * (y4 - y2) - (y3 - y1) * (x4 - x2);
      //real sg = -1; 
      //   if (a0 > 0) sg = 1;
      //localM2d[0][0] = sg * (2 * a0 + a1 + a2) / 72.; localM2d[0][1] = sg * (2 * a0 + a1 + a2) / 72.; localM2d[0][2] = sg * (2 * a0 + a1 + a2) / 72.; localM2d[0][3] = sg * (2 * a0 + a1 + a2) / 72.;
      //localM2d[1][0] = sg * (2 * a0 + a1 + a2) / 72.; localM2d[1][1] = sg * (2 * a0 + a1 + a2) / 72.; localM2d[1][2] = sg * (2 * a0 + a1 + a2) / 72.; localM2d[1][3] = sg * (2 * a0 + a1 + a2) / 72.;
      //localM2d[2][0] = sg * (2 * a0 + a1 + a2) / 72.; localM2d[2][1] = sg * (2 * a0 + a1 + a2) / 72.; localM2d[2][2] = sg * (2 * a0 + a1 + a2) / 72.; localM2d[2][3] = sg * (2 * a0 + a1 + a2) / 72.;
      //localM2d[3][0] = sg * (2 * a0 + a1 + a2) / 72.; localM2d[3][1] = sg * (2 * a0 + a1 + a2) / 72.; localM2d[3][2] = sg * (2 * a0 + a1 + a2) / 72.; localM2d[3][3] = sg * (2 * a0 + a1 + a2) / 72.;
      real x1 = mesh->knots[bound->knots_num[1]]->x - mesh->knots[bound->knots_num[0]]->x,
         y1 = mesh->knots[bound->knots_num[1]]->y - mesh->knots[bound->knots_num[0]]->y,
         z1 = mesh->knots[bound->knots_num[1]]->z - mesh->knots[bound->knots_num[0]]->z;
      real x2 = mesh->knots[bound->knots_num[2]]->x - mesh->knots[bound->knots_num[0]]->x,
         y2 = mesh->knots[bound->knots_num[2]]->y - mesh->knots[bound->knots_num[0]]->y,
         z2 = mesh->knots[bound->knots_num[2]]->z - mesh->knots[bound->knots_num[0]]->z;
      real x3 = mesh->knots[bound->knots_num[1]]->x - mesh->knots[bound->knots_num[3]]->x,
         y3 = mesh->knots[bound->knots_num[1]]->y - mesh->knots[bound->knots_num[3]]->y,
         z3 = mesh->knots[bound->knots_num[1]]->z - mesh->knots[bound->knots_num[3]]->z;
      real x4 = mesh->knots[bound->knots_num[2]]->x - mesh->knots[bound->knots_num[3]]->x,
         y4 = mesh->knots[bound->knots_num[2]]->y - mesh->knots[bound->knots_num[3]]->y,
         z4 = mesh->knots[bound->knots_num[2]]->z - mesh->knots[bound->knots_num[3]]->z;

      real aa = sqrt(x1 * x1 + y1 * y1 + z1 * z1);   // k0 - k1
      real ba= sqrt(x2 * x2 + y2 * y2 + z2 * z2);   // k0 - k2
      real ca = sqrt(x3 * x3 + y3 * y3 + z3 * z3);   // k3 - k1
      real da = sqrt(x4 * x4 + y4 * y4 + z4 * z4);   // k3 - k2
      real area = .25 * sqrt((-aa + ba + ca + da) * (aa - ba + ca + da) * (aa + ba - ca + da) * (aa + ba + ca - da));
      
      for (int i = 0; i < 4; i++)
         b[bound->knots_num[i]] += bound->value * area / 4.;
   }
}

void FEM::CreateSLAE()
{
   hexahedron* hexa;
   for (int i = 0; i < num_of_FE; i++)
   {
      hexa = mesh->hexas[i];
      CreateG(hexa);
      CreateM(hexa);
      CreateA(hexa);
      Createb(hexa);
   }

   AddSecondBounds();
   //for (int i =0 ; i < num_of_knots; i++) std::cout << b[i] << "\n";
   AddFirstBounds();
   //WriteMatrix(A);
}

void FEM::CreateA(hexahedron* hexa)
{
   for (int i = 0; i < 8; i++ )
      for (int j = 0; j < 8; j++ )
         AddElement(A, hexa->knots_num, hexa->knots_num[i], hexa->knots_num[j], localG[i][j]);
         //localA[i][j] = localG[i][j] + localM[i][j];
   //AddLocal(A, hexa->knots_num, localA, 1);
}

void FEM::CreateM(hexahedron* hexa)
{
    for (int i = 0; i < 8; i++)
      for (int j = 0; j < 8; j++)
         localM[i][j] = Integrate(Mij, i, j, hexa->knots_num);
}

void FEM::CreateG(hexahedron* hexa)
{
   hexa->lam = 1.;
   for (int i = 0; i < 8; i++)
      for (int j = 0; j < 8; j++)
         localG[i][j] = hexa->lam * Integrate(Gij, i, j, hexa->knots_num);
}

void FEM::Createb(hexahedron* hexa) 
{
   real localb[8]{};
   real f_[8]{};
   for (int i = 0; i < 8; i++)
      f_[i] = f(mesh->knots[hexa->knots_num[i]]);
   
   for (int i = 0; i < 8; i++)
      for (int j = 0; j < 8; j++)
         localb[i] += localM[i][j] * f_[j];
   
   for (int i = 0; i < 8; i++)
      b[hexa->knots_num[i]] += localb[i];
}

real FEM::scalar(real* u, real* v, int size)
{
   real sum = 0;
   for (int i = 0; i < size; i++)
      sum += v[i] * u[i];

   return sum;
}

void FEM::MatxVec(real* v, Matrix* A, real* b)
{
   //real *out = new real[num_of_knots];

   for (int i = 0; i < num_of_knots; i++)
      v[i] = A->d[i] * b[i];

   for (int i = 0; i < num_of_knots; i++)
      for (int j = A->ig[i]; j < A->ig[i + 1]; j++) // -1?
      {
         v[i] += A->l[j] * b[A->jg[j]];
         v[A->jg[j]] += A->u[j] * b[i];
      }
}

void FEM::SolveSLAE( )
{
   real res, alpha, beta, skp, eps = 1e-17;
   int i, k;
   x = q;

   real lastres;
   MatxVec(ff, A, x);
   for (i = 0; i < num_of_knots; i++)
      z[i] = r[i] = b[i] - ff[i];
   MatxVec(p, A, z);
   res = sqrt(scalar(r, r, num_of_knots)) / sqrt(scalar(b, b, num_of_knots));

   for (k = 1; k < 100000 && res > eps; k++)
   {
      lastres = res;
      skp = scalar(p, p, num_of_knots);
      alpha = scalar(p, r, num_of_knots) / skp;
      for (i = 0; i < num_of_knots; i++)
      {
         x[i] += alpha * z[i];
         r[i] -= alpha * p[i];
      }
      MatxVec(ff, A, r);
      beta = -scalar(p, ff, num_of_knots) / skp;
      for (i = 0; i < num_of_knots; i++)
      {
         z[i] = r[i] + beta * z[i];
         p[i] = ff[i] + beta * p[i];
      }
      res = sqrt(scalar(r, r, num_of_knots)) / sqrt(scalar(b, b, num_of_knots));
   }
   std::cout << "iter: " << k << " Residual: " << res << std::endl;
}

int FEM::mu(int index)
{
   return index % 2;
}

int FEM::v(int index)
{
   return (index / 2) % 2;
}

int FEM::nu(int index)
{
   return index/ 4;
}

real FEM::W(int index, real alpha)
{
    if (!index) return 1. - alpha;
    return alpha;
}

real FEM::d_phi(int index, int var, real ksi, real etta, real tetha)
{
   real d_phi = 0.;
   switch (var)
   {
      case 0:    // ksi
      {
         d_phi = W(v(index), etta) * W(nu(index), tetha);
         if (!mu(index)) d_phi *= -1;
         break;
      }
      case 1:     // etha
      {
         d_phi = W(mu(index), ksi) * W(nu(index), tetha);
         if (!v(index)) d_phi *= -1;
         break;
      }
      case 2:     // theta
      {
         d_phi = W(mu(index), ksi) * W(v(index), etta);
         if (!nu(index)) d_phi *= -1;
         break;
      }
   }

   return d_phi;
}

real FEM::prime_by_var(int varOnCube, int varOnFE, int knot_num[8], real ksi, real etta, real tetha)
{
   real var = 0.;
   for (int i = 0; i < 8; i++)
   {
      switch (varOnFE)
      {
         case 0: var += mesh->knots[knot_num[i]]->x * d_phi(i, varOnCube, ksi, etta, tetha); break;
         case 1: var += mesh->knots[knot_num[i]]->y * d_phi(i, varOnCube, ksi, etta, tetha); break;
         case 2: var += mesh->knots[knot_num[i]]->z * d_phi(i, varOnCube, ksi, etta, tetha); break;
      }
   }
   return var;
}

real FEM::phi(int index, real ksi, real etta, real tetha)
{
   return  W(mu(index), ksi) * W(v(index), etta) * W(nu(index), tetha);
}

real FEM::det_J()
{
   return J[0][0] * J[1][1] * J[2][2] + J[2][0] * J[0][1] * J[1][2] + J[1][0] * J[2][1] * J[0][2] 
       - (J[2][0] * J[1][1] * J[0][2] + J[0][0] * J[2][1] * J[1][2] + J[1][0] * J[0][1] * J[2][2]);
}

real FEM::Integrate(const std::function<real(real, real, real, int, int, int[8])> f, int i, int j, int knot_nums[8])
{
   const int nKnot = 3;//5; // Knots num

   const real xj[nKnot] 
   = { .7745966692414833, 0., -.7745966692414833 }; // sqrt(0.6)
   //= { -sqrt(5. + 2. * (sqrt(10. / 7.))) / 3., -sqrt(5. - 2. * (sqrt(10. / 7.))) / 3.,	// Scales
   //           0. , sqrt(5. - 2. * (sqrt(10. / 7.))) / 3. , sqrt(5. + 2. * (sqrt(10. / 7.))) / 3. };

   const real qj[nKnot] 
   = { .55555555555555555, .8888888888888888, .55555555555555555 };
   //= { (322. - 13. * sqrt(70.)) / 900., (322. + 13. * sqrt(70.)) / 900., 128. / 225.,	// Weights
   //               (322. + 13. * sqrt(70.)) / 900., (322. - 13. * sqrt(70.)) / 900. };

   real result = 0.;
   for (int ix = 0; ix < nKnot; ix++)
      for (int iy = 0; iy < nKnot; iy++)
         for (int iz = 0; iz < nKnot; iz++)
            result += qj[ix] * qj[iy] * qj[iz] * (f(.5 + xj[ix] * .5, .5 + xj[iy] * .5, .5 + xj[iz] * .5, i, j, knot_nums));
   return result / 8.; 
}

real FEM::Integrate2D(const std::function<real(real, real, int, int, int[4])> f, int i, int j, int knot_nums[4])
{
   const int nKnot = 3;//5; // Knots num

   const real xj[nKnot]
      = { .7745966692414833, 0., -.7745966692414833 }; // sqrt(0.6)
      //= { -sqrt(5. + 2. * (sqrt(10. / 7.))) / 3., -sqrt(5. - 2. * (sqrt(10. / 7.))) / 3.,	// Scales
      //           0. , sqrt(5. - 2. * (sqrt(10. / 7.))) / 3. , sqrt(5. + 2. * (sqrt(10. / 7.))) / 3. };

   const real qj[nKnot]
      = { .55555555555555555, .8888888888888888, .55555555555555555 };
   //= { (322. - 13. * sqrt(70.)) / 900., (322. + 13. * sqrt(70.)) / 900., 128. / 225.,	// Weights
   //               (322. + 13. * sqrt(70.)) / 900., (322. - 13. * sqrt(70.)) / 900. };

   real result = 0.;
   for (int ix = 0; ix < nKnot; ix++)
      for (int iy = 0; iy < nKnot; iy++)
            result += qj[ix] * qj[iy] * (f(.5 + xj[ix] * .5, .5 + xj[iy] * .5, i + 1, j + 1, knot_nums));
   return result / 4.;
}

void FEM::calc_grad(int ij, int index, real ksi, real etta, real tetha)
{
   switch (ij)
   {
      case 1: 
         for (int i = 0; i < 3; i++)
            gradi[i] = d_phi(index, i, ksi, etta, tetha);
         break;
      case 2:
         for (int i = 0; i < 3; i++)
            gradj[i] = d_phi(index, i, ksi, etta, tetha);
         break;
   }
}

void FEM::WriteMatrix(Matrix* A)
{
   double** mat = new double* [num_of_knots] {};
   for (int i = 0; i < num_of_knots; i++)
   {
      mat[i] = new double[num_of_knots] {};
   }

   for (int i = 0; i < num_of_knots; i++)
   {
      mat[i][i] = A->d[i];
      for (int j = A->ig[i]; j < A->ig[i + 1]; j++)
      {
         mat[i][A->jg[j]] = A->l[j];
         mat[A->jg[j]][i] = A->u[j];
      }
   }

   std::ofstream out("matrix.txt");

   for (int i = 0; i < num_of_knots; i++)
   {
      for (int j = 0; j < num_of_knots; j++)
      {
         out.setf(std::ios::left);
         out.width(15);
         out << mat[i][j];
      }
      out << "\n";
   }
}