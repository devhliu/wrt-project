# Weighted ray / Radon transforms in 3D

<p float="center">
  <img src="https://github.com/fedor-goncharov/Weighted-ray-Radon-transforms-in-3D/blob/master/pictures/k_comparison_output.gif" width="400" />
  <img src="https://github.com/fedor-goncharov/Weighted-ray-Radon-transforms-in-3D/blob/master/pictures/shepp_logan_reduction.gif" width="400" />
</p>

This is project is a part of my phd thesis conducted under the supervision of professor [Roman Novikov](http://www.cmap.polytechnique.fr/~novikov/)

Here we aim to develop new inversion methods for weighted (generaliezed) Radon transforms. 
The latter are of particular importance in various applications in the domain of inverse 
problems (e.g., in tomographies, geophysics). In particular, we work on methods which 
could be numerically more stable against the noise in tomographical data compared to existing methods. 

The very precise theoretical explanation of given algorithms is in [[2,3]](http://www.cmap.polytechnique.fr/~fedor.goncharov/publications.html).
Also my personal page is [here](http://www.cmap.polytechnique.fr/~fedor.goncharov/).

## Sructure of the project

Here I briefly explain for what each program is intended for. The details about their input/output/parameters/usage/.../ 
you can find in '(respective_folder)/README.md'.

  * #### nfft_inversion -- Octave/Matlab scripts for inversion of Radon transforms in 2D/3D using projection theorem 
  
  * #### radon_analytic -- computations of Radon transforms in 3D of a test-function with compact support in 
  a three-dimensional unit ball; the main feature is that the expression of the test-function must be given by an analytical expression 
  in a separate C file
  
  * #### radon_grid -- computations of Radon transforms of a test-function with compact support in 
  a three-dimensional unit ball; here the test-function must be given by values on a discrete grid in a three-dimensional cube 
  
  * #### ray_analytic -- computations of ray transforms in 3D in a layer-by-layer sampling scheme of a test-function with compact support in a three-dimensional unit ball (see also the README.md inside project for information about 'layer-by-layer sampling scheme'); the test-function must be given by an analytical expression in a separate C file
  
  * #### ray_grid -- computations of ray transforms in 3D in a layer-by-layer sampling scheme of a test-function with compact support in a three-dimensional unit ball (see the README.md inside project for information about 'layer-by-layer sampling scheme');
  test-function must be given by its values on a discrete grid in a three-dimensional cube
  
  * #### reduction_ray_radon -- reduces the data given by ray transforms in 3D in layer-by-layer sampling scheme to the 
  data given by Radon transforms in 3D
  
## Description of programs 

// TODO 

## Future plans

// TODO

