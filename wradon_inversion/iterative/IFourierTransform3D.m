% Computes inverse Fourier transform of function 'u' 
% given in the unit cube [-1, 1]^3; where ngrid is the number of points
% on the unit dimension of the cube including endpoints

function iFu = IFourierTransform3D(ngrid, u)
  
  dx = 2.0 / (ngrid - 1);    % real space shift
  ldomain = 2.0 + dx;        % size of extended space domain
  dfreq = 1 / (ngrid * dx);  % frequency step to obtain DFT
  
  shift_frequencies_xi1 = zeros(ngrid, ngrid, ngrid);  
  shift_frequencies_xi2 = zeros(ngrid, ngrid, ngrid);
  shift_frequencies_xi3 = zeros(ngrid, ngrid, ngrid);
  for i = 1 : ngrid
    shift_frequencies_xi1(:, i, :) = (i - 1) * dfreq;
    shift_frequencies_xi2(i, :, :) = (i - 1) * dfreq;
    shift_frequencies_xi3(:, :, i) = (i - 1) * dfreq;
  endfor
  
  shift_argument = shift_frequencies_xi1 .+ shift_frequencies_xi2 ...
                 .+ shift_frequencies_xi3;
  shift_operator = exp((-2)*pi * 1i * ldomain * shift_argument);
  
  iFu = (dx^3) * ifftn(u);  % ifft on the unit cube
  iFu = fftshift(iFu .* shift_operator);
  
endfunction