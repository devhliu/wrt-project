#ifndef __H_RGRID_INTERPOLATION
#define __H_RGRID_INTERPOLATION

  double cube_trilinear_interp(double*** values, int ngrid,
    double x, double y, double z);

  double cube_zero_interp(double*** values, int ngrid,
    double x, double y, double z);
  
#endif