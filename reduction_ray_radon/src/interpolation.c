double cell_linear_interp(double x, double v1, double v2) {
  return ((1-x)*v1 + x*v2); 
}
  
double cell_bilinear_interp(double x, double y, double v1, double v2, double v3, double v4) {
  double c1 = cell_linear_interp(x, v1, v2),
	 c2 = cell_linear_interp(x, v4, v3);
	 
  return cell_linear_interp(y, c1, c2);
}
