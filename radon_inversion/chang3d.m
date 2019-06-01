function reconstruction = chang3d(ngrid, w0, input_filename, nphi, ntheta, nshift, rsupp, padd)
  
  [nodes, values, jweights] = rtft3D(input+filename, nphi, ntheta, nshift, rsupp, padd);
  reconstruction = nfft_reconstruct_3d(ngrid, nodes, values, jweights) ./ w0;
  
endfunction