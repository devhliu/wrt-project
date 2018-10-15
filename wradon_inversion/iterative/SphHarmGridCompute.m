% function computes values of spherical harmonics on a cubic grid [-1,1]^3 

function SphHarmonicsGrid = GenSphHarmonics(ngrid, degree)
  
  linx = linspace(-1, 1, ngrid);
  liny = linspace(1, -1, ngrid);  % axis Y must be inverted
  linz = linspace(-1, 1, ngrid);
  [X,Y,Z] = meshgrid(linx, liny, linxz);                                        % grid in 3D
                                                                               
  [PHI, THETA, ] = cart2sph(X, Y, Z);                                           % grid in spherical coordinates
  THETA = THETA - pi/2;                                                         % shift latitude angle from (-pi/2, pi/2) ...
                                                                                % to (0,pi)                                                                  
  SphHarmonicsGrid = zeros(ngrid, ngrid, ngrid, (degree + 1)*(2*degree + 1) - 1); 
                                                                                % array of spherical harmonics on the grid
  current_harmonic = 1;  % iterator for spherical harmonics
  for k = 1 : degree
    for m = -k : k
      current_harmonic_func = @(theta, phi) Ykm(k, m, theta, phi); % lambda function for a given harmonic
      SphHarmonicsGrid(:, :, :, current_harmonic) = arrayfun(current_harmonic_func, THETA, PHI);
      current_harmonic = current_harmonic + 1;
    endfor
  endfor
endfunction
  