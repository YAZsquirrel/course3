#include "Filtration.h"
#include "GridMaker.h"
#include <fstream>
using namespace mesh_comps;
using namespace filtration;

namespace filtration {

	Mesh* mesh = new Mesh();

	void StartMakingMesh()
	{
		std::ifstream f_filtr("FiltrParams.txt");
		int comps_num;
		f_filtr >> comps_num;
		for (int i = 0; i < comps_num; i++)
		{
			component* comp = new component();
			f_filtr >> comp->num >> comp->proportion;

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
			f_filtr >> heights[i] >> comp_num >> ph.saturation >> ph.viscosity >> ph.penetrability;
			ph.compsInPhase = comps[comp_num];
			phases.push_back(ph);
		}



		mesh->set_layers(heights, phases_num);
		delete[] heights;


		mesh->GenerateMesh();
	}
	
}