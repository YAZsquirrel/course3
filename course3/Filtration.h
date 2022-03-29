#pragma once
#include "GridMaker.h"
#include <fstream>
#include <vector>
typedef double real;
using namespace mesh_comps;

namespace filtration
{
   struct component
   {
      int num;          // num of a component
      real proportion;  // in a phase
   };

   struct phase
   {
      real h;
      real viscosity;      // etta
      real penetrability;  // k (kappa)
      real oil_over_water;     // ksi

      std::vector<component*> compsInPhase; 
      static int L;                   // size of comp array  
   };

   struct poroda
   {
      real K;
      real Fi;
   };

   class Filtration
   {
      public:
      //Filtration();
      poroda por;
      Mesh* mesh;
      real* flow_in_face;  // size face_size
      std::vector<phase> phases;       // size M
      std::vector<component*> comps;    // size 1..2..3....

      real* alpha;   // size 
      real* beta;    // size knot_num
      void Start();
   };
}
