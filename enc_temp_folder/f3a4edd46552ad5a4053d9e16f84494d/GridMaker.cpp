#include "GridMaker.h"
#include <fstream>
#include <cmath>
namespace mesh_comps
{
   Mesh::Mesh()
   {
      std::ifstream fparam("EnvParams.txt");
      // get size of an env and steps
      {
         real x, y, z;
         // TODO: delete un
         fparam >> x;

         fparam >> x >> y >> z;
         env_corner1 = knot(x, y, z);
         fparam >> x >> y >> z;
         env_corner2 = knot(x, y, z);
         fparam >> x >> y >> z;
         step = knot(x, y, z);
      }

      // get info on wells
      {
         int well_num;
         real x, y, v, well_rad, h1, h2; // v = tetha,  
         fparam >> well_num;
         w_info.wells.reserve(well_num);

         for (size_t i = 0; i < well_num; i++)
         {
            fparam >> x >> y >> h1 >> h2 >> v >> well_rad;
            w_info.wells.push_back( well( x, y, well_rad, h1, h2, v ) );
         }
         fparam >> w_info.conc_rad >> w_info.rad_knots >> w_info.conc;
         w_info.rad_knots = std::max(std::ceil(w_info.rad_knots / step.x) * step.x, std::ceil(w_info.rad_knots / step.y) * step.y);
      }

      fparam.close();
   }

   void Mesh::GenerateMesh()
   {
      GenerateKnots();

   }

   void Mesh::GenerateKnots()
   {
      // env_corner1, env_corner2, step
      std::ofstream bounds1("FirstBounds.txt");
      std::ofstream bounds2("SecondBounds.txt");
      std::ofstream fknots("Knots.txt");

      xn = (int)((env_corner2.x - env_corner1.x) / step.x) + 1;
      yn = (int)((env_corner2.y - env_corner1.y) / step.y) + 1;
      zn = (int)((env_corner2.z - env_corner1.z) / step.z) + 1;

      int plane_size = xn * yn;
      int total_knots = plane_size;// * zn; // 
      int* well_start_indecies = new int[w_info.wells.size()];
      std::vector<int>** inwell_indecies = new std::vector<int>*[w_info.wells.size()];
      knots.reserve(total_knots);
      real hx = step.x, hy = step.y, hz = step.z;
      bool is_near_well = false;

      real z = env_corner1.z;
      for (int k = 0; k < zn; k++)
      {
         real z = k * hz;
         for (int l = 0; l < layers.size(); l++)            
            if (z + hz > layers[l] && z < layers[l])
            {
               z = layers[l];
               /*k++;*/
               break;
            }

         if (k == 0)
         {
            for (int j = 0; j < yn; j++)
            {
               real y = j * hy;
               for (int i = 0; i < xn; i++)
               {
                  is_near_well = false;
                  real x = i * hx;
                  for (int m = 0; m < w_info.wells.size() && !is_near_well; m++)
                     is_near_well = abs(w_info.wells[m].x - x) < w_info.conc_rad && abs(w_info.wells[m].y - y) < w_info.conc_rad;
                  if (!is_near_well)
                     knots.push_back(new knot(x, y, z));
                  else 
                     plane_size--;
               }
            }

            for (int wk = 0; wk < w_info.wells.size(); wk++)
            {
               well *w = &w_info.wells[wk];
               real lx = w->x - w_info.conc_rad - step.x, 
                    rx = w->x + w_info.conc_rad + step.x, 
                    yu = w->y + w_info.conc_rad + step.y, 
                    yd = w->y - w_info.conc_rad - step.y;
            
               int well_area_x = ceil(2 * w_info.conc_rad / step.x),
                   well_area_y = ceil(2 * w_info.conc_rad / step.y); // +1
               bool x1 = true, y1 = true;
               int rays = 2 * (well_area_y + well_area_x);

               std::vector<int> inwell_inds;
               inwell_inds.reserve(rays + 4);

               for (int kn = 0; kn < knots.size() && (x1 || y1) && ( w->x + 2. * step.x > knots[kn]->x || w->y + 2. * step.y > knots[kn]->y); kn++)
               {
                  if (abs(knots[kn]->x - w->x) / abs(w->x) < 1e-7 && x1)
                  {
                     well_area_x++;
                     x1 = !x1;
                  }
                  if (abs(knots[kn]->y - w->y) / abs(w->y) < 1e-7 && y1)
                  {
                     well_area_y++;
                     y1 = !y1;
                  }
                  if (knots[kn]->x < rx && knots[kn]->x > lx && knots[kn]->y < yu && knots[kn]->y > yd)
                     inwell_inds.push_back(kn);
               }

               // sort clockwise
               std::vector<int> *sorted_inds = new std::vector<int>;
               sorted_inds->reserve(inwell_inds.size());
               for (int iw = 0; iw < well_area_x + 1; iw++)
                  sorted_inds->push_back(inwell_inds[iw]);
               for (int iw = well_area_x, iw1 = 1; iw1 < well_area_y;  iw1++)
                  sorted_inds->push_back(inwell_inds[iw + 2 * iw1]);
               for (int iw = 0; iw < well_area_x + 1; iw++)
                  sorted_inds->push_back(inwell_inds[inwell_inds.size() - 1 - iw]);
               for (int iw = well_area_x, iw1 = well_area_y -1 ; iw1 > 0; iw1--)
                  sorted_inds->push_back(inwell_inds[iw + 2 * iw1 - 1]);

               inwell_indecies[wk] = sorted_inds;

               rays = sorted_inds->size();
               if (wk == 0) well_start_indecies[0] = knots.size();
               else well_start_indecies[wk] = well_start_indecies[wk-1] + w_info.rad_knots * rays;

               plane_size += w_info.rad_knots * rays;
               //knots.resize(plane_size);

               
               for (int iw = 0; iw < sorted_inds->size(); iw++)
               {
                  real xiw = knots[sorted_inds->at(iw)]->x - w->x, yiw = knots[sorted_inds->at(iw)]->y - w->y;
                  real len = sqrt(xiw * xiw + yiw * yiw);
                  xiw /= len; yiw /= len;
                  //xiw += w->x;  yiw += w->y;

                  knots.push_back(new knot(xiw * w->radius + w->x, yiw * w->radius + w->y, z));

                  real length = sqrt((xiw * w->radius + w->x - knots[sorted_inds->at(iw)]->x) * (xiw * w->radius + w->x - knots[sorted_inds->at(iw)]->x)
                                   + (yiw * w->radius + w->y - knots[sorted_inds->at(iw)]->y) * (yiw * w->radius + w->y - knots[sorted_inds->at(iw)]->y));
                  real stepr = length / w_info.rad_knots;

                  //for (float r = stepr, j = 1; j < num_of_knots_on_rad; r += stepr, j++)
                  for (real r = w_info.conc != 1. ? length * (1. - w_info.conc) / (1. - pow(w_info.conc, w_info.rad_knots)) : stepr,
                     rr = 1; rr < w_info.rad_knots; rr++)
                  {
                     knots.push_back(new knot(xiw * (w->radius + r) + w->x, yiw * (w->radius + r) + w->y, z));
                     if (w_info.conc == 1.) r += stepr;
                     else r *= w_info.conc;
                  }
               }
               
            }

            //total_knots *= zn;
            //knots.resize(total_knots);
         }
         else
            for (int i = 0; i < plane_size; i++)
               knots.push_back(new knot(knots[i]->x, knots[i]->y, z));
      }

      fknots << knots.size() << '\n';
      for (int i = 0; i < knots.size(); i++)
         fknots << knots[i]->x << " " << knots[i]->y << " " << knots[i]->z << '\n';
      fknots.close();
      FindAllHexas(plane_size, well_start_indecies, inwell_indecies);

      delete[] well_start_indecies;
   }

   void Mesh::FindAllHexas(int plain_size, int* well_inds, std::vector<int>** inwell_inds)
   {
      std::vector<int*> lower_faces;
      std::vector<int*> upper_faces;
      std::ofstream fhexas("Hexahedrons.txt");

      lower_faces.reserve(2 * plain_size);
      upper_faces.reserve(2 * plain_size);

      bool inwell_prevx = false;
      for (int i = 0, y = 0, x = 0; i < well_inds[0] && y < yn; i++, x++)
      {
         if (x > xn - 1)
         {
            x = 0;
            y++;
         }

         int *face = new int[4]{};
         //bool inwell = false;
         bool found = false;
         //for (int iw = 0; iw < w_info.wells.size() && !inwell; iw++)
         //   for (int iwi = 0; iwi < inwell_inds[iw]->size() && !inwell; i++)
         //      if (inwell_inds[iw]->at(iwi) == i && inwell_inds[iw]->at(iwi) == i + 1)
         //         inwell_prevx = inwell = true;
         
         face[0] = i;
         if (abs(knots[i]->x - knots[i + 1]->x) <= (step.x + 1e-10) && abs(knots[i]->y - knots[i + 1]->y) < 1e-10 )  // o -> o
            face[1] = i + 1;
         for (int j = 0; j < well_inds[0]; j++)
            if (abs(knots[i]->y - knots[j]->y) <= (step.y + 1e-10) && abs(knots[i]->x - knots[j]->x) < 1e-10)
            {  
               face[2] = j;
               if (abs(knots[j]->x - knots[j + 1]->x) <= (step.x + 1e-10) && abs(knots[j]->y - knots[j + 1]->y) < 1e-10)
               {
                  face[3] = j + 1;
                  found = true;
                  break;
               }
            }

         if (found)
            lower_faces.push_back(face);
      }
      for (int w = 0; w < w_info.wells.size(); w++)
      {
         for (int i = well_inds[w], iw = 0; iw < inwell_inds[w]->size(); i += w_info.rad_knots, iw++) // r = inwell_inds[w]->size() - колво лучей, r * nr = колво в 
            for (int r = 0; r < w_info.rad_knots; r++)
            {
               int* face = new int[4]{};
               face[0] = i + r; face[1] = i + w_info.rad_knots + r;
               if (r = w_info.rad_knots - 1)
                  { face[2] = inwell_inds[w]->at(iw); face[3] = inwell_inds[w]->at(iw) + 1; }
               else
                  { face[2] = i + r + 1; face[3] = i + r + 1 + w_info.rad_knots; }
               lower_faces.push_back(face);
            }

         for (int r = 0; r < w_info.rad_knots; r++)
         {
            int* face = new int[4]{};
            face[0] = well_inds[w] + w_info.rad_knots * (inwell_inds[w]->size() - 1) + r; face[1] = well_inds[w] + r;
            if (r = w_info.rad_knots - 1)
            {
               face[2] = inwell_inds[w]->at(inwell_inds[w]->size() - 1); 
               face[3] = inwell_inds[w]->at(0);
            }
            else
            {
               face[2] = well_inds[w] + w_info.rad_knots * (inwell_inds[w]->size() - 1) + r + 1; 
               face[3] = well_inds[w] + w_info.rad_knots * inwell_inds[w]->size() + r + 1;
            }
            lower_faces.push_back(face);
         }
      }
      //for (int i = well_inds[w_info.wells.size() - 2]; i < well_inds[w_info.wells.size() - 1]; i++)
      fhexas << plain_size * (zn - 1) << '\n';

      for (int k = 0; k < zn - 1; k++)
      {
         real z = k * step.z;
         for (int l = 0; l < layers.size(); l++)
            if (z + step.z < layers[l] && z > layers[l])
            {
               z = layers[l];
               //k++;
               break;
            }

         for (int i = 0; i < lower_faces.size(); i++)
         {
            //delete[] upper_faces[i];
            upper_faces.clear();
            int *face = new int[4]{ lower_faces[i][0] + plain_size, lower_faces[i][1] + plain_size,
                                    lower_faces[i][2] + plain_size, lower_faces[i][3] + plain_size };
            upper_faces.push_back(face);
         }
         for (int i = 0; i < lower_faces.size(); i++)
         {
            fhexas << lower_faces[i][0] << " " << lower_faces[i][1] << " " << lower_faces[i][2] << " " << lower_faces[i][3]
                   << " " << upper_faces[i][0] << " " << upper_faces[i][1]
                   << " " << upper_faces[i][2] << " " << upper_faces[i][3] << '\n';

            int* face = new int[4]{ upper_faces[i][0], upper_faces[i][1],
                                    upper_faces[i][2], upper_faces[i][3]};

            delete[] lower_faces[i];
            lower_faces.clear();
            lower_faces.push_back(face);
         }
      }
      fhexas.close();
   }

}