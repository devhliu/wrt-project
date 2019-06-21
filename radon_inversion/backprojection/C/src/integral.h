#ifndef __H_INTEGRALS
#define __H_INTEGRALS

double radon3d_adjoint(double*** radon_values, double x, double y, double z, 
	double* phi, double* theta, double* theta_weights, double* shift, int nphi, int ntheta, int nshift);

#endif