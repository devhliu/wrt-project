%Reconstruction of a test function in 3D from its Fouirer transforms

function test_function = nfft_reconstruct_3d(ngrid, nodes, values, jacobian_weights)
% nfft_reconstruct_ftrt3d_tfunc3d.m
% depends NFFT library (see Chemnitz-TU for NFFT)
  
% Script performs reconstruction of a function from its Fourier transforms at 'nodes' 
% with values 'values'. Mathematically it works as a discretized version of (inverse) Fourier integral 
% in 3D. Result is given as a 3D matrix of size : ngrid x ngrid x ngrid, which, in turn, is a 
% grid on [-1,1)x[-1,1)x[-1,1). 

%Usage of the script
% ngrid : number of points (as intervals) in [-1.0, 1.0)
% nodes : matrix of size (number_of_nodes x 3), where '3' stands for [x, y, z], these are the 
%          points in space where Fourier transform is known
% values : vector of size (number_of_nodes) values of Fourier transform of a function in 'nodes'
% jacobian_weights : volumes at 'nodes' in the Riemann summ of discretized Fourier integral

% Jacobian weightening
  summands = jacobian_weights .* values; 

% Inversion using NFFT
  
  %normalization of 'nodes' so that they belong to torus [-1/2, 1/2)^3
  deltax = 2.0 /  ngrid; 
  nodes_normalized = nodes * deltax;
  nodes_normalized = nodes_normalized';

  % Initialisation of plan (for sufficient amount of memory)
  plan=nfft_init_3d(ngrid, ngrid, ngrid, size(nodes_normalized,2)); 
  
  %Initialization plan for the case of lack of memory
  % not clear what is n and what is cutoff
  %n=2^(ceil(log(max([ngrid; ngrid; ngrid]))/log(2))+1);
  %plan=nfft(3, [ngrid; ngrid; ngrid;], size(nodes_normalized,1), n,n,n,15,bitor(PRE_FG_PSI,bitor(FG_PSI,NFFT_OMP_BLOCKWISE_ADJOINT)),FFTW_MEASURE);
  
  % set nodes in plan
  nfft_set_x(plan, nodes_normalized);
  
  % precomputations
  nfft_precompute_psi(plan);
  
  %set Fourier coefficients
  nfft_set_f(plan, summands);
  
  %inversion
  nfft_adjoint(plan);
  
  %return function values  
  test_function = reshape(real(nfft_get_f_hat(plan)), ngrid, ngrid, ngrid);
  
  nfft_finalize(plan);
  
  %complicated turns (consequences of reshape)
  test_function = permute(test_function, [2 3 1]);  
  test_function = flipdim(test_function, 3);
  test_function = flipdim(test_function, 1);
 
end