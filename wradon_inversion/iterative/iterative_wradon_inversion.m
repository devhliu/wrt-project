% Reconstruction of a test function in 3D from weighted Radon transforms 
% using iterative approach of F. Goncharov, R. Novikov (2017)

function test_function = IterativeWradonInv3D(general_params, init_data_params, ...
  weights_data_params, stopping_params)
% iterative_wradon_inversion.m
% depends on read_func_grid3d.m, stopping_rule_check.m, operatorQ.m

% general_params
%     ngrid -- integer, number points on the grid per dimension (on [-1,1], including the ends) 
% init_data_params
%     file_name_init -- file with the data of the test-function on the grid [ngrid x ngrid x ngrid]
% weight_data_params:
%     degree -- integer, highest order of harmonics used
%     file_name_array -- array of files where coefficients of the weight are stored
%                        coefficients are tensors [ngrid x ngrid x ngrid]
% stopping_params
%     flag -- FLAG_NITERATIONS, FLAG_PRECISION, FLAG_COMPOSITE
%     flag_values{1,2} -- integer, double
%         FLAG_NITERATIONS -- algorithm stops after 'stopping_iteration' iterations
%                             stopping_iteration = flag_values{1} 
%         FLAG_PRECISION   -- algorithm stops reaching given precision |u_n - u_(n-1)| / | u_n-1 | < stopping_precision
%                             stopping_precision = flag_values{2}
%         FLAG_COMPOSITE   -- algorithm stops either after 'stopping_iteration' iterations either 
%                             reaching given precision |u_n - u_(n-1)| / | u_n-1 | < stopping_precision
%                             stopping_iteration = flag_values{1} 
%                             stoppping_precision = flag_values{2}
%
%     for a good choice of number of harmonics algo converges geometrically to quasi-optimum
%     set eps as precision to quasi-solution - it implies the exact number of iterations you need to do
%     to reach the quasi-optimum 

% ----------------------------read initial function ----------------------------
  ngrid = general_params{1};                                                    % size of the grid
  F0 = ReadFuncGgrid3D(init_data_params{1});                                    % read initial function from file

  if (degree == 0) 
    test_function = F0;
    return;
  endif
% ---------------------------- read weight coeffs. -----------------------------
  degree = weight_params{1};                                                    % read order of expansion used
  weights = zeros(ngrid, ngrid, ngrid, (degree + 1)*(2*degree + 1));
  weights_files = weight_params{2};                                             % files with weights
  for i = 1 : (degree + 1)*(2*degree + 1)                                       % read weights (corresponding to harmonics)
    weights(:, :, i) = ReadFuncGrid3D(ngrid, weights_files(i));                                                                            
  endfor
 
% -- Precomputations (normalization of weights, spherical harmonics matrices) --

% normalization of weights (w_{k,n} / w_{0,0}), precomputation of all harmonics
  weights_normalized = zeros(ngrid, ngrid, ngrid, (degree + 1)*(2*degree + 1) - 1);
  weights_normalized = weights(:, :, :, 2 : end) ./ weights(:, :, :, 1);
  harmonics = GenSphHarmonics(ngrid, degree);                         % spherical harmonics on [ngrid x ngrid x ngrid] 
                                                                      % of all degrees/orders
                                                                      
  u_new = zeros(ngrid, ngrid, ngrid);  % two layer iterative scheme -- future function layer
  u_old = F0;                          % current function layer

  iteration = 0;  % iteration counter
  
% ------------------------ Cycle computations ----------------------------------
  do 
    u_new = F0 + IntOperatorQ(ngrid, degree, weights_normalized, harmonics, u_old);
    u_old = u_new;
    iteration = iteration + 1;
  
    cycle_params = cell(2);
    cycle_params{1} = iteration; 
    cycle_params{2} = norm(u_new - u_old, 2) / norm(u_old, 2);                  % current precision from the new iteration
  while (stopping_rule(stoppping_params, cycle_params) == false)
# ------------------------ end of cycle computations ---------------------------
  test_function = u_old ./ weights_coeffs(:, :, :, 1);                          % return value of the test_function
endfunction 

