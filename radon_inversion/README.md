### Radon inversion in 2D/3D

Octave/Matlab scripts for inversion of classical Radon transforms in 2D/3D via projection theorem. 

  * **RtFt_2d.m** - Fourier transform of Radon transforms in 2D along shift variable
  * **RtFt_3d.m** - Fourier transforms of Radon transforms in 3D along shift variable
  * **lgwt.m** - script for computation of Gaussian nodes in interval [-1, 1] (used in RtFt_3d.m for computations
              of latitude angles on the two-dimensional sphere)
  * **nfft_reconstruct_3d.m** - application of inverse Fourier transform to the data returned by RtFt_3d.m


#### Dependencies

You must have Matlab/Octave interface for [NFFT](https://www-user.tu-chemnitz.de/~potts/nfft/) installed. 
The link to this interface (mex-libraries) can be found [here](https://www-user.tu-chemnitz.de/~potts/nfft/download.php).

Download Matlab/Octave binaries from the link given above and make path in your Matlab/Octave to the libraries (so that NFFT functions could be called locally).

#### Usage

 Load lgwt.m to your Matlab/Octave.

 * **RtFt_2d.m**
 
 * **RtFt_3d.m**
             [nodes, values, jacobian_weights] = RtFt_3d(filename, nphi, ntheta, nshift, rsupp, padding_coeff = 4)
