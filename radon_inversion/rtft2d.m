% 1D Fourier transform of the data given by Radon transforms in 2D

function [nodes, values, jacobian_weights] = rtft2d(filename, nphi, nshift, rsupp, padding_coeff = 4)

% This script reads data given by ray transforms in 2D from file and performs 
% 1D Fourier transform along shift argument. Returns values (values) of the 
% 1D Fourier integral and respective nodes (nodes) with volumes (jacobian_weights). 

% Usage:

%   filename          : file with Radon transforms [sigma, phi]
%   nphi              : number of projections in azimuth angle
%   nshift            : number of rays per one direction, shift [-1, 1]
%   rsupp             : radius of the support of the test-function
%   padding_coeff (4) : number of frequencies per direction (padding_coeff / dshift)
 
% NOTE: It is assumed that angles 'phi' and 'shifts' are uniformly spaced 
% in their intervals.
  
  data = csvread(filename); 
  rt = data(: , size(data)(2));                 % values of the weighted ray transforms; always take last column
  data_matrix = reshape(rt, nphi, nshift);      % reshape data into matrix [theta, phi, shift]
  
  dphi = 2*pi/ nphi;
  dshift = 2 * rsupp / (nshift-1);
  ntotal = padding_coeff * floor(1 / dshift);   % total number of points with zero padding (should be odd)
  dfreq   = 1 / (ntotal * dshift);              % discretization step in frequency domain
  
  phi_vec = (0 : nphi-1)' * dphi;               % directions 'phi' on circle [0, 2*pi)
  frequencies = (0 : ntotal-1)' * dfreq;
  frequencies_centered = (-ntotal/2 : ntotal/2)' * dfreq;

  nodes = [];            % nodes in frequency domain
  values = [];           % values of Fourier integral in nodes 
  jacobian_weights = []; % jacobian multipliers for each value
  
  ft_shift_correction = exp(2 * pi * 1i * rsupp * frequencies);    % correction of FFT due to displacement [-1, 1] -> [0, 2]

  for i_phi = 1 : nphi
      % angle phi
      phi = phi_vec(i_phi);
      
      % append nodes
      direction = [cos(phi) sin(phi)];
      nodes = [nodes; frequencies_centered * direction];
      
      % append jacobian
      jacobian_weight = 0.5 * dfreq * dphi * abs(frequencies_centered);
      jacobian_weights = [jacobian_weights; jacobian_weight];
      
      % append 1D Fourier transform of RT using FFT along fixed direction
      rt_at_direction = vec(data_matrix(i_phi, :));           % data on fixed direction
      fft_vec = fft(rt_at_direction, ntotal);                          % 1D Fourier integral using 1D FFT
      ft1d_vec = dshift * fftshift(fft_vec .* ft_shift_correction);    % centralizing frequencies and jacobian correction
      ft1d_vec = [ft1d_vec; ft1d_vec(1)];                              % append to the end periodic frequency
      
      values = [values; vec(ft1d_vec)];       % 'values' of the Fourier transform of test function at 'nodes'
  end
end 
