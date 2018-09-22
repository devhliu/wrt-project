### Radon inversion in 2D/3D

Octave/Matlab scripts for inversion of classical Radon transforms in 2D/3D via projection theorem. 
Generally speaking, functions here take an input a CSV file containing Radon transforms and inverts them via 
projection theorem (1D Fourier transform + 3D inverse Fourier transform).

  * **RtFt_2d.m** - Fourier transform of Radon transforms in 2D along shift variable
  * **RtFt_3d.m** - Fourier transforms of Radon transforms in 3D along shift variable
  * **lgwt.m** - script for computation of Gaussian nodes in interval [-1, 1] (used in RtFt_3d.m)
  * **nfft_reconstruct_3d.m** - application of inverse Fourier transform to the data returned by RtFt_3d.m


#### Dependencies

   You must have Matlab/Octave interface for [NFFT](https://www-user.tu-chemnitz.de/~potts/nfft/) installed. 
   The link to this interface (mex-libraries) can be found [here](https://www-user.tu-chemnitz.de/~potts/nfft/download.php).

   Download Matlab/Octave binaries from the link given above and make path in your Matlab/Octave to the libraries 
   (so that NFFT    functions could be called locally).

#### Usage

 Load lgwt.m to your Matlab/Octave.

 * **RtFt_2d.m**
 
       [nodes, values, jacobian_weights] = RtFt_2d(filename, nphi, ntheta, nshift, rsupp, padding_coeff = 4)
       
       Script reads data given by Radon transforms in 3D from file and performs 
       1D Fourier transform along shift variable. Returns arguments : nodes, values, jacobian_weights. 
       
       nodes : points in 3D frequency space where Fourier transforms is evaluated (size Nx3)
       values : values of Fourier transform in nodes (size Nx1complex)
       jacobian_weights : volume associated to each node in frequency space

       Usage of the script
         filename          : file where the data is stored in CSV format (sigma, phi, theta, rt)
         nphi              : number of projections in azimuth angle [0, 2*pi)
         ntheta            : number of projections in polar angle (0, pi)
         nshift            : number of hyperplanes per one direction; shifts uniformly vary [-1, 1]
         rsupp             : radius of the support of the test function
         padding_coeff (4) : parameter to padd Radon transforms with zeros along shifts
       
       Angles 'phi' are uniform on the circle and angles 'theta' correspond to Gaussian quadrature points, 
       i.e. theta_j = arccos(t_j), (t_j, j = 1, ntheta) - Gauss-Lebato points on [-1, 1] 

 
 * **RtFt_3d.m**
 
       [nodes, values, jacobian_weights] = RtFt_3d(filename, nphi, ntheta, nshift, rsupp, padding_coeff = 4)
 
 * **nfft_reconstruct_3d.m**
 
       test_function = nfft_reconstruct_3d(ngrid, nodes, values, jacobian_weights)

