#pragma once
#include <vector>
typedef double real;

namespace filtration
{
   struct component
   {
      int num;          // num of a component
      real proportion;  // in a phase
   };

   struct phase
   {
      real viscosity;      // etta
      real penetrability;  // k (kappa)
      real saturation;     // Sm

      component* compsInPhase; 
      static int L;                   // size of comp array  
   };

   real* flow_in_face;  // size face_size

   std::vector<phase> phases;       // size M
   std::vector<component*> comps;    // size 1..3...

   real* alpha;   // size 
   real* beta;    // size knot_num
}
