### Inversions in 2D/3D via backprojection algorithms

<p float="center"> 
  <img src="https://github.com/fedor-goncharov/Weighted-ray-Radon-transforms-in-3D/blob/master/pictures/backprojectionv1.png" width="360" />
  <img src="https://github.com/fedor-goncharov/Weighted-ray-Radon-transforms-in-3D/blob/master/pictures/backprojectionv0.png" width="360" />
  <br>
  <em> Adjoint of the Radon transform (left) and applied Laplacian in 3D (right) for Radon data corrupted with strong Poisson noise </em>
</p>


Octave/Matlab/C codes for computation of the adjoint Radon transforms in 2D/3D. In 3D the inversion of Radon transforms 
can be implemented as a backprojection algorithm when the adjoint transform is applied first and then the laplacian is applied after. In 2D the after the adjoint transform the inversion is not given by a laplacian (but a Hilbert transform), however 
for completeness reasons we put here code for laplacian in 2D. 


  * **adjrad3d.m** - Matlab/Octave scrpit for computation of the adjoint Radon transform in 3D
  * **laplace2dfft.m** - spectral laplacian in 2D
  * **laplace3dfft.m** - spectral laplacian in 3D
  * **hilbert1dfft.m** - spectral Hilbert transform in 2D

#### Usage

 * **adjrad3d.m**
 
       output_adjoint = adjrad3d(filename, ngrid, nphi, ntheta, nshift, interp_method=1, rsupp=1.0, exp_coeff=2.0)
       
       Script reads Radon transforms in 3D from file and performs 1D Fourier transforms 
       along shift variable. 
       
       Returned values
         nodes            : points in frequency domain where Fourier transforms is evaluated (size Nx3)
         values           : values for Fourier transform in nodes (size Nx1complex)
         jacobian_weights : volumes associated to each node in frequency domain

       Usage of the script
         filename          : file where the data is stored in CSV format
                             Data is expected in the following format : "[shift], [phi], [theta], [value]\n",
                             Variables 'shift, phi, theta' must vary in the following order : 
                 
                                 for (shift) 
                                    for (phi) 
                                       for (theta)
                                        ....
                                       end
                                    end
                                 end
                                      
         nphi              : number of azimuth nodes in [0, 2*pi)
         ntheta            : number of polar nodes in (0, pi)
        
                             Angles 'phi' are uniform on the circle and angles 'theta' 
                             correspond to Gaussian quadrature points, i.e., theta_j = arccos(t_j), 
                             (t_j, j = 1, ntheta) - Gauss-Lebato points on [-1, 1]. 

         nshift            : number of hyperplanes per one direction
                             Shifts are uniform along [-1,1]
         rsupp             : radius of the support of the test function
         padding_coeff (default=4) : parameter to padd Radon transforms with zeros along shift
 
 * **laplace2dfft.m**
 
       output_laplace = laplace2dfft(f, ngrid, period)
       
       Script reads Radon transforms in 2D from file and performs 1D Fourier transforms 
       along shift variable. 
       
       Returned values
         nodes            : points in frequency domain where Fourier transforms is evaluated (size Nx2)
         values           : values of Fourier transform in nodes (size Nx1complex)
         jacobian_weights : volumes associated to each node in frequency domain

       Usage of the script
         filename          : file where the data is stored in CSV format
                             Data modeled by Radon transforms is expected in the 
                             following format : "[shift], [phi], [value]\n",
                             Variables 'shift, phi' vary in the following order : 
                 
                                 for (shift) 
                                    for (phi) 
                                        ...
                                    end
                                 end
                                      
         nphi              : number of azimuthal angles [0, 2*pi)
                             Angles 'phi' are assumed to be uniform on [0, 2*pi).

         nshift            : number of hyperplanes per one direction
                             Shifts are assumed uniform on [-1,1].
         rsupp             : radius of the support of the test function
         padding_coeff(default=4) : parameter to padd Radon transforms with zeros along shift
 
 * **laplace3dfft.m**
 
       output_laplace = laplace3dfft(f, ngrid, period)
       
       Script performs reconstruction of a function from its 'values' of Fourier transforms 
       computed at 'nodes'.

       Returned values 
         test_function : real-valued matrix of size (ngrid x ngrid x ngrid)

       Usage of the script
        ngrid : number of points in [-1.0, 1.0); ngrid must be even;
        
        nodes            : matrix of size (NNodesx3); these are the points 
                           in 3d-space where the Fourier transform of the function is known. 
        values           : vector(NNodesx1(complex)) of values of Fourier transforms;
        jacobian_weights : volumes assigned nodes in the discretization of the Riemann integral for 
                           the Fourier transforms
## Future plans

  * Matlab/Octave code works is not parallelized and works sufficiently slow (it could be optimized, using the special 
  structure of the backprojection integral)
  
  * On the images border artifcats are present -- this is due to the property of a spectral laplacian (the argument is 
  not a periodic function) and can be fixed in a few ways -- either to compute adjoint on a bigger image or use sampling 
  by Chebyshev polynomials or apply periodization techniques (Lionel Moisan -- Periodic plus smooth image decomposition).
  
  * Backprojection for Radon transforms could be efficiently realized in C with help of GSL or other matrix libraries. 
