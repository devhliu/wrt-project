function adjoint = AdjointRadonTransform3D(filename, ngrid, method, nphi, ntheta, nshift, rsupp) 
  
  % read (weighted) Radon transforms from file
  rtransforms = csvread(filename)(: , 4);                                 % values of the (weighted) Radon transforms
  rtransforms = reshape(rtransforms, ntheta, nphi, nshift);               % reshape data into matrix [theta, phi, shift]  
  
  % set standard variables
  dphi   = 2 * pi / nphi;
  dshift = 2 * rsupp / (nshift-1);
  
  % set grid in Radon space
  angles_phi = (0 : nphi-1) * (2*pi / nphi);                            % directions 'phi' on circle [0, 2*pi)
  [gauss_nodes, theta_weights] = lgwt(ntheta, -rsupp, rsupp);        % Gaussian nodes on [-1,1] and weights
  angles_theta = acos(gauss_nodes);                                     % directions 'theta' on (0, pi)
  clear gauss_nodes;
  base_shifts = linspace(-rsupp, rsupp, nshift);
  
  adjoint = zeros(ngrid^3, 1);
  
  % create grid matrices
  lin = linspace(-rsupp, rsupp, ngrid);
  [X,Y,Z] = meshgrid(lin, lin, lin);
  XX = reshape(X, ngrid^3, 1); 
  YY = reshape(Y, ngrid^3, 1);
  ZZ = reshape(Z, ngrid^3, 1);
  pts = [XX YY ZZ]; % size ngrid^3 x 3
  
  % integral over sphere
  for i_phi = 1 : nphi
    for i_theta = 1 : ntheta
      
      % get angles
      phi = angles_phi(i_phi);
      theta = angles_theta(i_theta);
      
      % compute scalar products with direction for all points in the grid
      direction = vec([sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)]);
      shifts = pts * direction;
      
      % adjoint Radon t-m at all points of the grid for fixed direciton 
      radon_interpolated = interp1(base_shifts, vec(rtransforms(i_theta, i_phi, :)), shifts, method, 0);      
      
      % add term
      adjoint += radon_interpolated * dphi * theta_weights(i_theta);
      
    endfor
  endfor
  
  adjoint = reshape(test_function_ad, ngrid, ngrid, ngrid);
endfunction