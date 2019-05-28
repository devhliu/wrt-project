### Radon inversion in 2D/3D

Octave/Matlab scripts for inversion of classical Radon transforms in 2D/3D via Slice Projection Theorem. 
Functions here take on input a CSV file containing Radon transforms and perform reconstructions by applying 
two times Fourier transforms (1D Fourier transform and then 3D (or 2D) inverse Fourier transform).

  * **rtft2D.m** - Fourier transform of Radon transforms in 2D with respect to shift
  * **rtft3D.m** - Fourier transforms of Radon transforms in 3D with respect to shift
  * **nfft_reconstruct_3d.m** - inverse 3D Fourier transform of the data returned by rtft3D.m
  * **nfft_reconstruct_2d.m** - inverse 2D Fourier transofmr of the data returned by rtft2D.m


#### Dependencies

   You must have Matlab/Octave interface for [NFFT](https://www-user.tu-chemnitz.de/~potts/nfft/) installed. 
   The link to this interface (mex-libraries) can be found [here](https://www-user.tu-chemnitz.de/~potts/nfft/download.php).
   
   Also I used Matlab script 'lgwt.m' for Gauss-Legendre quadrature rule taken from here [here](https://fr.mathworks.com/matlabcentral/fileexchange/4540-legendre-gauss-quadrature-weights-and-nodes?requestedDomain=) and a script 'cartprod.m' 
   taken from [here](https://fr.mathworks.com/matlabcentral/fileexchange/5475-cartprod-cartesian-product-of-multiple-sets).

#### Usage

 * **rtft3d.m**
 
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
 
 * **rtft2d.m**
 
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
 
 * **nfft_reconstruct_3d.m**
 
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
        
* **nfft_reconstruct_2d.m**
 
       test_function = nfft_reconstruct_2d(ngrid, nodes, values, jacobian_weights)
       
       Script performs reconstruction of a function from its 'values' of Fourier transforms at 'nodes'. 
       Result is given as a 3D matrix of size ngrid x ngrid pixels. 

       Returns values 
         test_function : real-valued matrix of size (ngrid x ngrid x ngrid)

       Usage of the script
        ngrid : number of points in [-1.0, 1.0); ngrid must be even;
                Note that it is number of points on non-closed interval. 
                NFFT assumes that your signal is periodic, so the values on the missing edge points
                of the grid is can reconstructed from periodicity.  
        
        nodes            : matrix of size (Nnodesx3); these are the points 
                           in space where Fourier transform of signal is known. 
        values           : vector of size (Nnodesx1(complex)); these are the values of 
                           Fourier transform of the signal at nodes;
        jacobian_weights : volumes of cells related to nodes in the Riemann summ of discretized Fourier integral

