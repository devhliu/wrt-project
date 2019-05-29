#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

double ray_interpolate(double slice_z, double slice_shift, 
		       double* phi, int i_phi, int nphi,
		       double* shift, int nshift,
		       int ngrid,
		       double*** values) {

    //linear interpolation in variable z
    //spline interpolation in variable s
  
    if ( (fabs(slice_z) > 1.0) || (fabs(slice_shift) > 1.0) ) {
      return 0;
    }
    //----------------------------------------------- get z-cell ---------------------------------/
    const double delta_z = 2.0 / (ngrid - 1);
    
    int iz = (int)( (slice_z + 1.0) / delta_z);			  // number of the cell in z-slices
    if (iz == (ngrid - 1)) {
      iz -= 1;
    }
    //--------------------------------------------------------------------------------------------/
    
    //----------------------------------get s-cell -----------------------------------------------/
    const double delta_s = 2.0 / (nshift - 1);
    int is = (int)( (slice_shift + 1.0) / delta_s);	  // number of the cell in s
    if (is == (nshift - 1)) {
      is -= 1;
    }
    //--------------------------------------------------------------------------------------------//
    //--------------------------------interpolations in variables (s,z)---------------------------//
    gsl_interp_accel *acc;
    // spline interpolation in variable s (spline interpolation in shift argument)
    gsl_spline *spline;
	
    double interp_s_vals[3], z_args[3];
    //set interpolaztion z arguments (quadratic interpolation in z argument)
    z_args[0] = (-1.0) + delta_z * iz;
    z_args[1] = (-1.0) + delta_z * (iz + 1);
    z_args[2] = (-1.0) + delta_z * (iz + 2);
    
    //set interpolation s values at z-points
    acc = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(gsl_interp_cspline, nshift);

    // spline interpolations
    {
      gsl_spline_init(spline, shift, values[iz][i_phi], nshift);
      interp_s_vals[0] = gsl_spline_eval(spline, slice_shift, acc); 

      gsl_spline_init(spline, shift, values[iz + 1][i_phi], nshift);
      interp_s_vals[1] = gsl_spline_eval(spline, slice_shift, acc);

      if ( iz >= (ngrid - 2) ) {
        interp_s_vals[2] = 0;
      } else {
        gsl_spline_init(spline, shift, values[iz + 2][i_phi], nshift);
        interp_s_vals[2] = gsl_spline_eval(spline, slice_shift, acc);
      }
    }
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    
    
    //quadratic lagrange interpolation 
    acc = gsl_interp_accel_alloc();
    gsl_interp *quad_interp = gsl_interp_alloc(gsl_interp_polynomial, 3);
    gsl_interp_init(quad_interp, z_args, interp_s_vals, 3);
    
    double ray_transform_value = gsl_interp_eval(quad_interp, z_args, interp_s_vals, slice_z, acc);
    gsl_interp_free(quad_interp);
    gsl_interp_accel_free(acc);
    
    return ray_transform_value;
}


double radon_reduction(double*** ray_values, 
             double* shift, double* phi, double* theta,
             int nshift, int nphi, int ntheta, int ngrid,
             int i_shift, int i_phi, int i_theta) {

  // integration along the plane (shift, phi, theta) using values of ray integrals values[ngrid][nshift][nphi]
  
  double tau; 				    	    // parametrization of slicing of plane into rays
  const double dtau = (2.0 / (ngrid - 1))/(2.0);    // half of dshift for step across plane
  const int nsize = 4 * ngrid;  
  int i_tau;
  
  const double p_shift = shift[i_shift];	    // shift point
  const double p_phi = phi[i_phi];            // longitude angle phi
  const double p_theta = theta[i_theta];      // latitude angle theta
  
  double integral_value = 0,			// value of the Radon transform
	       add = 0;				          // term to add
  
  //TODO
  //optimize computation of the length of parameter tau

  for (i_tau = -(nsize / 2); i_tau < (nsize / 2) + 1; ++i_tau) {
      tau = (i_tau) * dtau;
      double slice_z = p_shift * cos(p_theta) + tau * sin(p_theta);	        // z=z(tau) section
      double slice_shift = p_shift * sin(p_theta) - tau * cos(p_theta);     // dist_z = dist_z(tau) distance to the ray in XY plane

      //cycle moves across the plane, tau -- distance between current line and 'starting line' on the plane
      //'starting line' -- line on the plane which is closest to the origin
      
      if ((fabs(tau *tau + p_shift *p_shift) > 1.0) || (fabs(slice_z) > 1.0) || (fabs(slice_shift) > 1.0)) {
	        add = 0;
      } else {
	        // get ray transform using interpolations; add term to the summ
	        add = ray_interpolate(slice_z, slice_shift, phi, i_phi, nphi, shift, nshift, ngrid, ray_values);
      }
      integral_value += (add * dtau);
  }
  return integral_value;
}
