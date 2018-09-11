#ifndef _H_SLICING_REDUCTION
#define _H_SLICING_REDUCTION

double Rt_slicing(int i_shift, int i_phi, int i_theta,
		  double* shift, double* phi, double* theta,
		  int nshift, int nphi, int ntheta, 
		  double*** values, 
		  int ngrid
);

#endif