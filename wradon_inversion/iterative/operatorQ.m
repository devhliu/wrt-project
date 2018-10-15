% Application of integral operator Q_{W,D,m} to function 'u', where 
% W is the weight of finite number of harmonics, D is a unit ball in 3D,
% m - highest degree of harmonics expansion; [Goncharov, Novikov, 2016]

function Qu = IntOperatorQ(ngrid, degree, weights_normalized, harmonics, u)
  % domain cutoff 
   U = cutoff .* u;
  % mutipliciation by normalized weights
   argVect = weights_normalized .* U; 
  % direct Fourier transforms
    ...
  % multiplication by harmonics in Fourier space 
   argVect = harmonics .* argVect;
  
  % summation of terms
   argVect = sum(argVect, 4);
  
  % inverse Fourier transform
   ...
  
  % domain cutoff  
   argVect = cutoff .* argVect;
   
  % normalization
  % TODO: check if it is required when phase factor is equal to 2*PI (see FFT)
  
  Qu = argVect;
  
endfunction 