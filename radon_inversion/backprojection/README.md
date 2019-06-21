### Radon inversions in 2D/3D via backprojection algorithms

<p float="center"> 
  <img src="https://github.com/fedor-goncharov/Weighted-ray-Radon-transforms-in-3D/blob/master/pictures/backprojectionv1.png" width="360" />
  <img src="https://github.com/fedor-goncharov/Weighted-ray-Radon-transforms-in-3D/blob/master/pictures/backprojectionv0.png" width="360" />
</p>


Octave/Matlab scripts for inversion of classical Radon transforms in 2D/3D via backprojection.
Functions here take on input a CSV file containing Radon transforms and perform reconstructions by applying 
first a dual Radon transform (backprojection) and then a filter is applied (R* and then spectral derivative - 
square root of Laplacian in 2D, Laplacian in 3D).

  * **adjrad3d.m** - computation of adjoint Radon transform in 3D
  * **laplace2dfft.m** - spectral laplacian in 2D
  * **laplace3dfft.m** - spectral laplacian in 3D

#### Usage

 * **adjrad3d.m**
 
       [nodes, values, jacobian_weights] = rtft2d(filename, nphi, ntheta, nshift, rsupp, padding_coeff = 4)
       
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
 
       [nodes, values, jacobian_weights] = rtft2d(filename, nphi, nshift, rsupp, padding_coeff = 4)
       
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
 
       test_function = nfft_reconstruct_3d(ngrid, nodes, values, jacobian_weights)
       
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
