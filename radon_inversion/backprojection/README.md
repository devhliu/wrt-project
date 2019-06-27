### Inversions in 2D/3D via backprojection algorithms

<p float="center"> 
  <img src="https://github.com/fedor-goncharov/Weighted-ray-Radon-transforms-in-3D/blob/master/pictures/backprojectionv1.png" width="360" />
  <img src="https://github.com/fedor-goncharov/Weighted-ray-Radon-transforms-in-3D/blob/master/pictures/backprojectionv0.png" width="360" />
  <br>
  <em> Adjoint of the Radon transform (left) and applied Laplacian in 3D (right) for Radon data corrupted with strong Poisson noise </em>
</p>


Octave/Matlab/C codes for computation of the adjoint Radon transforms in 2D/3D. In 3D the inversion of Radon transforms 
can be implemented as a backprojection algorithm when the adjoint transform is applied first and then the laplacian is applied after. In 2D the after the adjoint transform the inversion is not given by the laplacian (but by the Hilbert transform), however 
for completeness reasons we put here code for laplacian in 2D. 

The implementation of spectral laplacians is based on the notes from MIT from [here](https://math.mit.edu/~stevenj/fft-deriv.pdf) by Steven J. Johnson. 


  * **adjrad3d.m** - Matlab/Octave scrpit for computation of the adjoint Radon transform in 3D
  * **laplace2dfft.m** - spectral laplacian in 2D
  * **laplace3dfft.m** - spectral laplacian in 3D
  * **hilbert1dfft.m** - spectral Hilbert transform in 2D

#### Usage

 * **adjrad3d.m**
 
       The Radon data are assumed to be given on grid uniform in shifts and uniform-gauss in directions on the sphere.

       Input parameters : 
            filename : file with Radon transforms
            ngrid    : number of pixels in one direction in the XYZ-domain
            nphi     : number of longitude angles [0, 2pi]
            ntheta   : number of latitude angles [0, pi]
            nshift   : 
       interp_method : interpolation method, because Radon transforms are given on a discrete grid
                      an interpolated data is needed for computation of the adjoint integral
                      possible methdods : 1 is "linear" (default), 0 is "nearest"
 
            rsupp    : length of the XYZ-domain in one direction (by default = 1.0)
       expand_factor : compute adjoint transform on the grid on the size [expand_factor * ngrid]^3

       Return value : 
            A three-dimensional matrix is returned of size (expand_factor x ngrid) per dimension, 
            where elements stand for the adjoint Radon transform at given pixel.
       
       A remark : 
  
       Radon data in file is expected to be ordered as follows : 
                                for (shift) 
                                    for (phi) 
                                       for (theta)
                                         fprintf("%f\n", radon(f, shift, phi, theta)
                                       end
                                    end
                                 end
                                 
 
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
