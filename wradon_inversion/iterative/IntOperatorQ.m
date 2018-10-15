% Application of integral operator Q_{W,D,m} to function 'u', where 
% W is the weight of finite number of harmonics, D is a unit ball in 3D,
% m - highest degree of harmonics expansion; [Goncharov, Novikov, 2016]

function Qu = IntOperatorQ(ngrid, degree, weights_normalized, harmonics, cutoff, u)
  
  % domain cutoff 
   U = cutoff .* u;
   argVect = zeros(ngrid, ngrid, ngrid, (degree + 1)*(2*degree + 1) - 1); % initialization of an argument vector
  
  % 1) mutipliciation by normalized weights
   argVect = weights_normalized .* U; 
  
  % 2) direct Fourier transforms
  for i = 1 : ((degree + 1)*(2*degree + 1) - 1)
    argVect(:, :, i) = FourierTransform3D(ngrid, argVect(:, :, i));
  endfor
  
  % 3) multiplication by harmonics in Fourier space 
   argVect = harmonics .* argVect;
  
  % 4) summation of terms
   argVect = sum(argVect, 4);
   
  % 5) inverse Fourier transform
   argVect = real(IFourierTransform3D(ngrid, argVect));
  
  % 6) domain cutoff  
   argVect = cutoff .* argVect;
   
  % normalization (check if it is required)
  % argVect = (1/2*pi) * argVect;
  
  % output 
   Qu = argVect;
endfunction 