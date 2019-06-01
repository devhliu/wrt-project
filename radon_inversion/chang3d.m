% chang3d.m

function reconstruction = chang3d(ngrid, w0, input_filename, nphi, ntheta, nshift, rsupp=1.0, padd=4.0)

% depends on rtft3d.m, 
%            nfft_reconstruct_3d.m
 
% The function takes on input values for weighted Radon transforms along 
% planes in 3D % from file (input_filename) and applies Chang-type formula.
% The formula is based on classical inversion of Radon transforms and here 
% it is implemented via the Fourier Slice Projection Theorem. 

% NOTE: for usage of Chang-type formula the values for w0 (zero order harmonic 
% of the weight) must be provided.

%% Usage of the script:

%   ngrid             : size of the output image in pixels (ngrid x ngrid x ngrid)
%   input_filename    : file with Radon transforms (order [shift, phi, theta])
%   nshift            : number of shifts in [-1,1] per one direction
%   nphi              : number of projections in azimuth angle (0, 2pi)
%   ntheta            : number of latitude angles in (0, pi)
%   rsupp             : radius of the support of the test-function (usually taken 1.0)
%   padding_coeff (4) : multiplying factor for appending data with zeros for each projection

% Output: voxel image of size ngrid x ngrid x ngrid
 
  [nodes, values, jweights] = rtft3d(input_filename, nphi, ntheta, nshift, rsupp, padd);
  reconstruction = nfft_reconstruct_3d(ngrid, nodes, values, jweights) ./ w0;
  
endfunction