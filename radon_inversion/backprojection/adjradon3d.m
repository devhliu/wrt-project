function adjoint = adjradon3d_test(filename, ngrid, nphi, ntheta, nshift, interp_method=1, rsupp=1.0, expand_factor=2.0) 
  
% adjradon3d.m 

% Function takes on input Radon transforms of a certain function and computes its adjoint Radon transform in 3D.  
% The Radon data are assumed to be given on grid uniform in shifts and uniform-gauss in directions on the sphere.

% INPUT PARAMETERS : 
%     filename : file with Radon transforms
%     ngrid    : number of pixels in one direction in the XYZ-domain
%     ntheta   : 
%     nshift   : 
%     interp_method : interpolation method, because Radon transforms are given on a discrete grid
%                     an interpolated data is needed for computation of the adjoint integral
%                     possible methdods : 1 is "linear" (default), 0 is "nearest"
%
%     rsupp    : length of the XYZ-domain in one direction (by default = 1.0)
%     expand_factor : compute adjoint transform on the grid on the size [expand_factor * ngrid]^3

% RETURN VALUE : 
%     A three-dimensional matrix is returned of size 'expand_factor*ngrid' where elements stand for the adjoint 
%     transform at the corresponding geometrical points.
  
  
  % read (weighted) Radon transforms from file
  rt = csvread(filename);                        % read values of the (weighted) Radon transforms
  rt = rt(:, size(rt)(2));
  rt = reshape(rt, ntheta, nphi, nshift);        % reshape into matrix [theta, phi, shift]  
  
  % set standard variables
  dphi   = 2 * pi / nphi;
  dshift = 2 * rsupp / (nshift-1);
  
  % set grid in Radon space
  % get equatorial angle
  angles_phi = (0 : nphi-1) * (2*pi / nphi);                         % directions 'phi' on circle [0, 2*pi)
  % get latitude angles (Gauss-Quadrature rule)
  [nodes_gauss, weights_theta] = lgwt(ntheta, -rsupp, rsupp);        % Gaussian nodes on [-1,1] and weights
  angles_theta = acos(nodes_gauss);                                  % directions 'theta' on (0, pi)
  clear nodes_gauss;
  shifts = linspace(-rsupp, rsupp, nshift);
 
  % create grid XYZ
  lin = (-(expand_factor/2)*ngrid : (expand_factor/2)*ngrid) * dshift;   % adjoint may be computed on bigger grid to avoid border effects in derivations
  L = length(lin);
  [XX,YY] = meshgrid(lin, lin);
  XX = reshape(XX, L^2, 1);
  YY = reshape(YY, L^2, 1);
  
  adjoint = zeros(L^3, 1);
  for i_z= 1 : L
    
    points = [XX YY ones(L^2,1)*lin(i_z)];                   % size of the grid L x L
    adjoint_slice = zeros(L^2, 1);
    
    printf("i_z=%d\n", i_z);
    fflush(stdout); 
    
  % integral over sphere for each z
    for i_phi = 1 : nphi    
      %printf("i_phi=%d\n", i_phi);
      %fflush(stdout);
      for i_theta = 1 : ntheta
      
        % get angles
        phi = angles_phi(i_phi);
        theta = angles_theta(i_theta);
      
        % compute scalar products with direction for all points in the grid
        direction = vec([sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)]);
        temp_shifts = points * direction;
      
        %leave points which are in the domain of integration
        temp_shifts_index = find((rsupp > abs(temp_shifts)));
        temp_shifts_interpolate = temp_shifts(temp_shifts_index); % shifts where interpolation is needed
      
        % adjoint Radon at all points of the z-slice for direciton (phi, theta)
        add_interpolated = qinterp1(shifts, vec(rt(i_theta, i_phi, :)), temp_shifts_interpolate, interp_method);
      
        % add
        add = zeros(L^2, 1);
        add(temp_shifts_index) = add_interpolated;
        adjoint_slice += add * dphi * weights_theta(i_theta);
        
      endfor
    endfor
    adjoint(1 + (i_z-1)*L^2 : i_z * L^2) = adjoint_slice;
  endfor 
  
  adjoint = reshape(adjoint, L, L, L);
endfunction