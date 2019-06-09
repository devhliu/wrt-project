## General information

This program computes Radon transforms in 3D along planes of a function which realized by an analytic expression in a separate file.  File "test_function.c" contains a template of such realization which should be implemented by the user prior to compile.

**IMPORTANT:** This program computes Radon transforms for all planes in 3D which belong to a certain grid. The planes are parametrized by triples (s, φ, θ), where s is the distance from the origin to the plane, (φ, θ) shperical angles for the normal vector of the plane. In particular, the normal vector is given by: (sin(θ)cos(φ), sin(θ)sin(φ), sin(θ)), φ is in [0,2π], θ is 
in (0, π). Variables s, φ are distributed uniformly in [-1,1], (0,2π), respectively. Variable θ corresponds to Gauss-quadrature rule: 
      θ_j = arccos(t_j), where t_j are uniform points on [-1,1].

The properties of grids (s, φ, θ) are set by user in the config file (binary/config.txt).

**OUTPUT FORMAT:** The Radon transforms are returned in a separate csv file in the following format: 

      for islice = (0: nslice-1)
        for iphi = (0: nshift-1)
          for itheta = (0: nhpi-1) 
             foutput(output_file, "%lf, %lf, %lf, %lf\n", z[islice], s[ishift], φ[iphi], Rf(z[islice], s[ishift], φ[iphi]);
          endfor
        endfor
      endfor



## Requirements 

The programs here are designed to work under Unix operating systems.  
To compile the project on your computer you need to have installed:  

GCC compiler, OpenMP libraries, GNU GSL libraries (+2.5)

## Compilation / Installation

  1) Go to 'src' directory
  
  2) Run Makefile
      ```
        make install
      ```
  3) Clean directory from object files (optional):
  
      ```
        make clean 
      ```
  If you want to use other test-function then you have to change the file
  'test_function.c' and repeat steps (1-2), possibly setting a new name for binary in the Makefile.


### Config file 

The config file (config.txt) is stored in /binary. 

The purpose of the config file for parameters (-p, --parameters [file] options) is to provide to the program information 
on the grid of rays. 

1. The second line contains the number of longitude angles φ which are uniformly spaced in [0,2π].
2. The second line contains the number of latitude angles θ which are distributed in [0,π] according to Gauss-quadrature rule. 
3. The third line contains a number shifts 's' for planes with a fixed normal. Shifts are spaced uniformly along [-1,1] (including endpoints -1,1).
4. The fourth line contains number of grid points 'ngrid' on each plane of integration. Each plane has size ngrid x ngrid. 
This parameter affects precision of integration.

The length of each line should not exceed 128 symbols

### Example of a config file

> 256			: number of longitude angles  
> 128			: number of latitude angles  
> 129			: number of steps per fixed direction 
> 129			: number of points on the plane in one dimension (plane will have grid of ngrid x ngrid points)

### Output

The output stored in a specified output file (-o --output) in a CSV format in the following order:  

**[s], [φ], [θ], [Rf(s,φ,θ)]**  

Example:  
> -1.000000, 0.000000, 3.122878, -0.009362  
> -1.000000, 0.000000, 3.098635, -0.009187  
> -1.000000, 0.000000, 3.074249, -0.008878

### Examples of test-functions

One can find in folder 'examples-test-functions' some examples of realizations of some simple test-functions.  
In order to try them, rename any of those files to 'test_function.c' and copy them to the main directory of 'radon_analytic'.     Then proceed with steps in 'Compilation / Installation' in order to obtain a compiled binary for a given
test-function. 


