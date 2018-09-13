% 1D Fourier transform of the data given by Radon transforms in 2D

function [nodes, values, jacobian_weights] = RtFt_2d(filename, nphi, nshift, rsupp, 
  padding_coeff = 4)

% This script reads data given by ray transforms in 2D from file and performs 
% 1D Fourier transform along shift argument. Returns arguments (nodes) of the 
% aforementioned 1D Fourier integral and its respective values (values). 

% Usage of the script
% filename          : file where the data is stored in the form [sigma, phi]
% nphi              : number of projections in azimuth angle [0, 2*pi)
% nshift            : number of rays per one direction
% rsupp             : radius of the support, where the test function is lying
% padding_coeff (4) : number of frequencies per direction (padding_coeff / dshift)
 
% Geometry assumes using uniform distribution of angles 'phi' and Gaussian polar 
% angles 'theta', i.e. theta_j = arccos(t_j), t_j - Gaussian points, j = 1, ntheta 
  
  data = csvread(filename); 
  rt = data(: , 4);                             % values of the weighted ray transforms
  data_matrix = reshape(rt, nphi, nshift);      % reshape data into matrix [theta, phi, shift]
  
 
  dphi = 2*pi/ nphi;
  dshift = 2 * rsupp / (nshift-1);
  ntotal = padding_coeff * floor(1 / dshift);   % total number of points with zero padding (should be odd)
  dfreq   = 1 / (ntotal * dshift);              % discretization step in frequency domain
  
  phi_vec = vec(0 : nphi-1) * dphi;    % directions 'phi' on circle [0, 2*pi)
  frequencies = vec(0 : ntotal-1) * dfreq;
  frequencies_centered = vec(-ntotal/2 : ntotal/2) * dfreq;

  nodes = [];            % nodes in frequency domain
  values = [];           % values of Fourier integral in nodes 
  jacobian_weights = []; % jacobian multipliers for each value
  
  ft_shift_correction = exp(2*pi*1i*rsupp * frequencies);    % correction of FFT due to displacement [-1, 1] -> [0, 2]

  for i_phi = 1 : nphi
      
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
      
      printf("Iteration phi %d\n", i_phi);
      fflush(stdout);
  end
end 