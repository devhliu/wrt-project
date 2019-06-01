function reconstruction = chang2d(ngrid, w0, input_filename, nphi, nshift, rsupp, padd)
  
  [nodes, values, jweights] = rtft2D(input_filename, nphi, nshift, rsupp, padd);
  reconstruction = nfft_reconstruct_2d(ngrid, nodes, values, jweights) ./ w0;
  
endfunction