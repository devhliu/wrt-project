% characteristic function of a unit ball in the unit cube in 3D

function cutoff = GenCutoffFunc(ngrid)
  
  lin = linspace(-1, 1, ngrid);
  [X,Y,Z] = meshgrid(lin, lin, lin);
  cuttoff = ((X.^2 + Y.^2 + Z.^2)< 1);
  
endfunction