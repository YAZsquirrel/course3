#include "FEM.h"
using namespace filtration;
int main()
{
	FEM *fem = new FEM();
	fem->SolveElliptic();
	fem->GetSolutiononPlane(5.1);
	//Mesh* mesh = new Mesh();
	//Filtration *filtr = new Filtration;
	//filtr->Start();
}
