#pragma once
#include <vector>
#include <set>
#include <map>
#include <list>
//#include "Filtration.h"
typedef double real;

namespace mesh_comps
{

   struct knot
   {
      //unsigned int knot_num;
      knot(real _x, real _y, real _z) : x(_x), y(_y), z(_z) {}
      knot() : x(0.0), y(0.0), z(0.0) {}
      real x, y, z;
   };

   struct well
   {
      real x, y;
      real radius;
      real h1, h2; 
      real intake;

      well(real _x, real _y, real rad, real _h1, real _h2, real th) :
         x(_x), y(_y), radius(rad), h1(_h1), h2(_h2), intake(th) {}
   };

   struct hexahedron
   {
      int knots_num[8];
      real lam;
   };

   class Mesh
   {
      public: 
	    Mesh();

       std::vector<knot*> knots;
       std::vector<hexahedron*> hexas;

       void set_layers(real* _layers, int layers_num)
       {
            layers.reserve(layers_num);
            for (int i = 0; i < layers_num; i++)
               layers.push_back(_layers[i]);
            layers.push_back(w_info.wells[0].h1);
            layers.push_back(w_info.wells[0].h2);
       }

       void GenerateMesh(); //

       private: 
       int xn = 0, yn = 0, zn = 0;
       std::vector<real> layers;
       knot env_corner1, env_corner2, step;
       struct well_info
       {
          std::vector<well> wells;
          real conc_rad, rad_knots, conc;
       } w_info;

       void GenerateKnots();
       void FindAllHexasAndBounds(int plain_size, int* well_inds, std::vector<int>** inwell_indecies, real* zs);
       //void 
   };

}