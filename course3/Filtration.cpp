#include "Filtration.h"

using namespace mesh_comps;
//using namespace filtration;

namespace filtration {

	void Filtration::Start()
	{
		mesh = new Mesh();
		std::ifstream f_filtr("FiltrParams.txt");
		real trash;
		int comps_num;
		f_filtr >> trash;

		f_filtr >> por.K >> por.Fi;

		f_filtr >> comps_num;
		for (int i = 0; i < comps_num; i++)
		{
			component* comp = new component();
			f_filtr >> comp->num;

			comps.push_back(comp);
		}

		int phases_num;
		f_filtr >> phases_num;
		phases.reserve(phases_num);

		real* heights;
		heights = new real[phases_num];
		for (int i = 0; i < phases_num; i++)
		{
			phase ph;
			int comp_num;
			f_filtr >> heights[i] >> ph.oil_over_water >> ph.viscosity >> ph.penetrability;
			ph.h = heights[i];
			ph.compsInPhase.push_back(new component);
			//ph.compsInPhase.push_back(); //comps[comp_num];
			//for (int j = 0; j < comp_num; j++)
			//{
			//	ph.compsInPhase[j].n
			//}
			phases.push_back(ph);
		}

		mesh->set_layers(heights, phases_num);
		delete[] heights;

		mesh->GenerateMesh();
	}
	
}