% 1D Fourier transform of the data given by Radon transforms in 3D

function [nodes, values, jacobian_weights] = RtFt_3d(filename, nphi, ntheta, nshift, rsupp, 
  padding_coeff = 4)
% depends on lgwt.m

% This script reads data given by Radon transforms in 3D from file and performs 
% 1D Fourier transform along shift argument. Returns arguments (nodes) of the 
% aforementioned 1D Fourier integral and its respective values (values). 

% Usage of the script
% filename : file where the data is stored in the form [sigma, phi, theta]
% nphi : number of projections in azimuth angle [0, 2*pi)
% ntheta : number of projections in polar angle (0, pi)
% nshift : number of hyperplanes per one direction
% rsupp : radius of the support, where the test function is lying
% padding_coeff (4) : number of frequencies per direction (padding_coeff / dshift)
 
% Geometry assumes using uniform distribution of angles 'phi' and Gaussian polar 
% angles 'theta', i.e. theta_j = arccos(t_j), (t_j, j = 1, ntheta) - Gaussian points on [-1, 1] 
  
  rt = csvread(filename)(: , 4);                                 % values of the weighted Radon transforms
  rt_matrix = reshape(rt, ntheta, nphi, nshift);                 % reshape data into matrix [theta, phi, shift]
  
  phi_vec = vec(0 : nphi-1) * (2*pi / nphi);                     % directions 'phi' on circle [0, 2*pi)
  [gauss_nodes, theta_weights] = lgwt(ntheta, -1.0, 1.0);        % Gaussian nodes on [-1,1] and weights
  theta_vec = acos(gauss_nodes);                                 % directions 'theta' on (0, pi)
  
  dphi = 2*pi / nphi;
  dshift = 2*rsupp / (nshift-1);
  ntotal = padding_coeff * floor(1 / dshift);   % total number of points with zero padding
  dfreq   = 1 / (ntotal * dshift);              % discretization step in frequency domain
  
  
  frequencies = vec(0 : ntotal-1) * dfreq;
  frequencies_centered = vec(-ntotal/2 : ntotal/2) * dfreq;

  nodes = [];               % array of nodes in frequency domain
  values = [];              % array of values of Fourier integral in nodes 
  jacobian_weights = [];    % array of jacobian multipliers for each node
  
  ft_shift_correction = exp(2*pi * 1i * rsupp * frequencies);    % correction of FFT due to displacement [-1, 1] -> [0, 2]
  
  % get array of nodes ---------------------------------------------------------
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
  end
  
  % get array of jacobian weights ----------------------------------------------
  for i_theta = 1 : ntheta 
    jacobian_weight = 0.5 * dfreq * dphi * theta_weights(i_theta) * (frequencies_centered.^2);
    jacobian_weight = repmat(jacobian_weight, nphi, 1);
    jacobian_weights = [jacobian_weights; jacobian_weight];
  end
  
  % get array of values --------------------------------------------------------
  
  % constant shift correction matrix of size ntotal x nphi
  ft_shift_correction_matrix = zeros(ntotal, nphi);
  ft_shift_correction_matrix(:, :) = repmat(ft_shift_correction, 1, nphi);
  
  for i_theta = 1 : ntheta
    
    % perform fft along dimension 'nshift'
    rt_theta_slice =  permute(rt_matrix(i_theta, :, :), [3 2 1]);                      % ntotal x nphi
    fft_theta_slice = dshift * fft(rt_theta_slice, ntotal, 1);                         % fft along dimension : shift
    fft_theta_slice = fftshift(fft_theta_slice .* ft_shift_correction_matrix , 1);     % fftshift along dimension : ntotal
    
    % add another layer in shift
    fft_full_theta_slice = zeros(ntotal + 1, nphi);                                    % zero matrix to append last layer
    fft_full_theta_slice(1 : ntotal, 1 : nphi) = fft_theta_slice;                      % copy the existing one
    fft_full_theta_slice(ntotal + 1, 1 : nphi) = fft_full_theta_slice(1, 1 : nphi);
    
    % reshape slice (theta=const) back to vector
    values_add = reshape(fft_full_theta_slice, nphi * (ntotal + 1), 1);
    values = [values; values_add];
  end
  % ----------------------------------------------------------------------------
end 
