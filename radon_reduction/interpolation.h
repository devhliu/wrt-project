#ifndef _H_INTERPOLATION_REDUCTION
#define _H_INTERPOLATION_REDUCTION

double cell_linear_interp ( double x, double v1, double v2 );

double cell_bilinear_interp ( double x, double y, double v1, double v2, double v3, double v4 );

#endif
