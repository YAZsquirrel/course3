#include "FEM.h"
#include <math.h>

//real FEM::lambda(real knot[2], real t)
//{
//	return real();
//}

real FEM::f(knot *point)
{
	real x = point->x;
	real y = point->y;
	real z = point->z;
	//return 0;
	switch (un)
	{
	case 1: return 1;
	case 2: return x + y + z;
	case 3: return 0;
	case 4: return -4;
	case 5: return x * x + y * y; //
	case 6: return -4 + x * x + y * y;
	default:
		return 0;
	}
}

real FEM::ug(knot *point)
{
	real x = point->x;
	real y = point->y;
	real z = point->z;

	switch (un)
	{
	case 1: return 1;
	case 2: return x + y + z; 
	case 3: return z;
	case 4: return x * x + y * y;
	case 5: return x * x + y * y; // x^3 + y^3 + xy + 1
	case 6: return x * x + y * y;
	default:
		return 0;
	}
}
