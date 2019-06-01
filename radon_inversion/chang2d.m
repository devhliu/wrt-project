% chang2d.m

function reconstruction = chang2d(ngrid, w0, input_filename, nphi, nshift, rsupp, padd)
  
% depends on rtft2d.m, 
%            nfft_reconstruct_2d.m
 
% The function takes on input values for weighted Radon transforms along 
% lines in 2D from file (input_filename) and applies Chang-type formula.
% The formula is based on the classical inversion of Radon transforms and here 
% it is implemented via the Fourier Slice Projection Theorem. 

% NOTE: for usage of Chang-type formula the values for w0 (zero order harmonic 
% of the weight) must be provided.

%% Usage of the script:

%   ngrid             : size of the output image in pixels (ngrid x ngrid)
%   input_filename    : file with Radon transforms (order [shift, phi])
%   nphi              : number of projections in azimuth angle (0, 2pi)
%   nshift            : number of shifts in [-1,1] per one direction
%   rsupp             : radius of the support of the test-function (usually taken 1.0)
%   padding_coeff (4) : multiplying factor for appending data with zeros for each projection

% Output: plane image of size ngrid x ngrid
  
  [nodes, values, jweights] = rtft2d(input_filename, nphi, nshift, rsupp, padd);
  reconstruction = nfft_reconstruct_2d(ngrid, nodes, values, jweights) ./ w0;
  
endfunction