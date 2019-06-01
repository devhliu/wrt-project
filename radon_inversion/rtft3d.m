% 1D Fourier transform of the data given by Radon transforms in 3D

function [nodes, values, jacobian_weights] = rtft3d(filename, nphi, ntheta, nshift, rsupp = 1.0, padding_coeff = 4.0)
% depends on lgwt.m, 
%            cartprod.m

% This script reads data given by Radon transforms in 3D from file and computes 
% 1D Fourier transforms along all projections. The returned result consists of nodes of the 
% aforementioned 1D Fourier integral, its values and 
% volume size (jacobian weights) of each node.

% Usage of the script

% filename : file where the data is stored in the form [sigma, phi, theta]
% nphi     : number of projections in azimuth angle [0, 2*pi)
% ntheta   : number of projections in polar angle (0, pi)
% nshift   : number of hyperplanes per one direction
% rsupp    : radius of the support, where the test function is lying
% padding_coeff (4) : number of frequencies per direction
 
% Geometry assumes using uniform distribution of angles 'phi', Gaussian polar 
% angles 'theta' (i.e. theta_j = arccos(t_j), (t_j, j = 1, ntheta) - Gaussian points on [-1, 1])
% and uniform shifts 
  
  % read input data
  rt = csvread(filename);                                        % values of the weighted Radon transforms
  rt = rt(:, size(rt)(2));                                       % always take last column
  rt_matrix = reshape(rt, ntheta, nphi, nshift);                 % reshape data into matrix [theta, phi, shift]
  clear rt;
  
  % set grid in Radon space
  phi_vec = vec(0 : nphi-1) * (2*pi / nphi);                     % directions 'phi' on circle [0, 2*pi)
  [gauss_nodes, theta_weights] = lgwt(ntheta, -1.0, 1.0);        % Gaussian nodes on [-1,1] and weights
  theta_vec = acos(gauss_nodes);                                 % directions 'theta' on (0, pi)
  
  % set variables for computing Fourier transform and output
  dphi   = 2 * pi / nphi;
  dshift = 2 * rsupp / (nshift-1);              % time step of the signal
  ntotal = padding_coeff * floor(1 / dshift);   % Nyguist frequency times total number of points with zero padding
  
  dfreq  = 1 / (ntotal * dshift);               % discretization step in frequency domain
  
  % array of frequencies per direction
  frequencies = vec(0 : ntotal-1) * dfreq;
  frequencies_centered = vec(centered_frequencies(ntotal, dfreq));

  nodes = [];               % array of nodes in frequency domain
  values = [];              % array of values of Fourier integral in nodes 
  jacobian_weights = [];    % array of jacobian multipliers for each node
  
  % correction of FFT due to displacement [-1, 1] -> [0, 2]
  ft_shift_correction = exp(2*pi * 1i * rsupp * frequencies);
  
  % array of nodes -------------------------------------------------------------
  printf("Generating nodes in Fourier space...\n");
  fflush(stdout);
  
  for i_theta = 1 : ntheta
    theta = theta_vec(i_theta);
    
    % long vector of directions : ( vec_f1.x, vec_f1.y, vec_f1.z, vec_f2.x, vec_f2.y, vec_f3.z ... ) 
    direction_theta = reshape([sin(theta)*cos(phi_vec'); sin(theta)*sin(phi_vec'); cos(theta)*ones(1,nphi)], 1, 3 * nphi);
    
    % matrix of frequencies (theta=const) : (frequencies x (direction_theta))
    nodes_theta = frequencies_centered * direction_theta;
    
    % each triple of columns are nodes for fixed 'phi' along shift : need to concatenate triples of columns vertically
    % reshape 'nodes_theta' to threed-dim matrice, where last index is a number of triple
    nodes_add_threedim = reshape(nodes_theta, size(frequencies_centered,1), 3, nphi);
    
    % trick to concatenate them vertically
    nodes_add = permute(nodes_add_threedim, [1 3 2]);
    nodes_add = reshape(nodes_add, [], size(nodes_add_threedim, 2), 1);
    
    nodes = [nodes; nodes_add];
  endfor
  
  printf("Primary nodes created. Appending nodes...");
  fflush(stdout);
  
  % stabilization - append zero nodes outside the ball of radius of Nyquist frequency / 2
  append_nodes = cartprod(frequencies_centered, frequencies_centered, frequencies_centered); % all points in Fourier domain
  append_nodes = [append_nodes, sqrt(sum(append_nodes.^2, 2))];                   % append norms column to the right
  append_nodes = append_nodes(append_nodes(:, 4) > ((ntotal/2 + 1) * dfreq), :);  % choose nodes outside the ball
  append_nodes = append_nodes(:, [1 2 3]);                                        % remove norm column
  nodes = [nodes; append_nodes];
  
  size_append = size(append_nodes, 1);  % remember the number of nodes that have been appended
  printf("%d primary nodes created. Appending nodes...", size(values, 1));
  fflush(stdout);
  clear append_nodes;
  
  printf("Computing jacobian weights for all nodes...");
  fflush(stdout);
  
  % jacobian weights ----------------------------------------------
  for i_theta = 1 : ntheta 
    jacobian_weight = 0.5 * dshift * dfreq * dphi * theta_weights(i_theta) * (frequencies_centered.^2);
    jacobian_weight = repmat(jacobian_weight, nphi, 1);
    jacobian_weights = [jacobian_weights; jacobian_weight];
  endfor
  
  jacobian_weights = [jacobian_weights; ones(size_append, 1)*(dfreq^3)];  % append jacobian weights (of appended edges)
  printf("Done.\n");
  fflush(stdout);
  
  % values --------------------------------------------------------
  printf("Computing values of for the Fourier transform at nodes...");
  fflush(stdout);
  
  % constant shift correction matrix of size ntotal x nphi
  ft_shift_correction_matrix = zeros(ntotal, nphi);
  ft_shift_correction_matrix(:, :) = repmat(ft_shift_correction, 1, nphi);
  
  for i_theta = 1 : ntheta
    % fft along dimension 'nshift'
    rt_theta_slice =  permute(rt_matrix(i_theta, :, :), [3 2 1]);                      % ntotal x nphi
    fft_theta_slice = fft(rt_theta_slice, ntotal, 1);                                  % fft along dimension : shift
    fft_theta_slice = fftshift(fft_theta_slice .* ft_shift_correction_matrix , 1);     % fftshift along dimension : ntotal
    
    % reshaping slice to vector
    values_add = reshape(fft_theta_slice, nphi * ntotal, 1);
    values = [values; values_add];
  endfor
  
  values = [values; zeros(size_append, 1)];
  printf("Done.\n");
  fflush(stdout);
  clear ft_shift_correction_matrix;
  % END ------------------------------------------------------------------------
  
  clear rt_matrix, phi_vec, theta_vec, frequencies_centered;
  
  printf("Done.\n");
  fflush(stdout); 
  
endfunction 
