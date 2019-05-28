% 1D Fourier transform of the data given by Radon transforms in 2D

function [nodes, values, jacobian_weights] = rtft2d(filename, nphi, nshift, rsupp, padding_coeff = 4)

% This script reads data given by ray transforms in 2D from file and performs 
% 1D Fourier transform along shift argument. Returns values (values) of the 
% 1D Fourier integral and respective nodes (nodes) with volumes (jacobian_weights). 

% Usage of the script:

%   filename          : file with Radon transforms [sigma, phi]
%   nphi              : number of projections in azimuth angle
%   nshift            : number of rays per one direction, shift [-1, 1]
%   rsupp             : radius of the support of the test-function
%   padding_coeff (4) : multiplying factor for appending data with zeros
 
% It is assumed that angles 'phi' and 'shifts' are uniformly spaced 
% in their intervals.
  
  data = csvread(filename); 
  rt = data(: , size(data)(2));                 % values of the weighted ray transforms; always take last column
  data_matrix = reshape(rt, nphi, nshift);      % reshape data into matrix [theta, phi, shift]
  
  dphi = 2*pi/ nphi;
  dshift = 2 * rsupp / (nshift-1);
  ntotal = padding_coeff * floor(1 / dshift);   % total number of points with zero padding (should be odd)
  dfreq   = 1 / (ntotal * dshift);              % discretization step in frequency domain
  
  phi_vec = vec((0 : nphi-1) * dphi);               % directions 'phi' on circle [0, 2*pi)
  frequencies = vec((0 : ntotal-1) * dfreq);
  frequencies_centered = vec(centered_frequencies(ntotal, dfreq));

  nodes = [];            % nodes in frequency domain
  values = [];           % values of Fourier integral in nodes 
  jacobian_weights = []; % jacobian multipliers for each value
  
  ft_shift_correction = exp(2 * pi * 1i * rsupp * frequencies);    % correction of FFT due to displacement [-1, 1] -> [0, 2]

  for i_phi = 1 : nphi
      % current angle phi
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
      
      values = [values; vec(ft1d_vec)];       % 'values' of the Fourier transform of test function at 'nodes'
  end
  
  % stabilization - append nodes outside of the ball of radius of Nyquist frequency
  append_nodes = cartprod(frequencies_centered, frequencies_centered); % all points in Fourier domain
  append_nodes = [append_nodes, sqrt(sum(append_nodes.^2, 2))];                            % append norms column to the right
  append_nodes = append_nodes(append_nodes(:, 3) > ((ntotal/2 + 1) * dfreq), :);           % choose nodes outside the ball
  append_nodes = append_nodes(:, [1 2]);                                                   % remove norm column
  nodes = [nodes; append_nodes];
  size_append = size(append_nodes, 1);  % remember the number of nodes that have been appended
  clear append_nodes;
  values = [values; zeros(size_append, 1)];                               % append zero values
  jacobian_weights = [jacobian_weights; ones(size_append, 1)*(dfreq^2)];  % append square volumes
    
end 
