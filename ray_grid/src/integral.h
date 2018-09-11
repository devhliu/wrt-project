#ifndef _H_RAY_GRID_INTEGRALS
#define _H_RAY_GRID_INTEGRALS

double ray_transform(double*** values, int ngrid_tfunc, 
		        double zslice,
			double shift, double phi, int ngrid);

#endif