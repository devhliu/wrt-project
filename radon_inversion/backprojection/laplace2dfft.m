function output = laplace2dfft(u, ngrid, period) 
  
  % Laplace transform of a three-dimensional image of size ngrid via FFT
  
  Y = fftn(u); % 3D DFT of the signal 
  frequence_factor = zeros(ngrid, 1);
  for k = 0 : ngrid-1
    if (k <= (ngrid/2))
      frequence_factor(k + 1) = k;
    else 
      frequence_factor(k + 1) = k - ngrid;
    endif
   endfor
   
  freq_factor_dim1 = zeros(ngrid, ngrid);  
  freq_factor_dim2 = zeros(ngrid, ngrid);
  
  [freq_factor_dim1, freq_factor_dim2] = meshgrid(frequence_factor, frequence_factor);
  
  U = -(2*pi / period)^2*(freq_factor_dim1.^2 + freq_factor_dim2.^2);
  
  output = ifftn(U .* Y);
  
endfunction