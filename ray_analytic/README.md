
## General description 

This program computes ray transforms in 3D of a function which realized by an analytic expression in a separate file.  
File "test_function.c" contains a template of such realization.

**IMPORTANT:** This program computes ray transforms not for all rays in 3D, but for the slice-by-slice scheme. 
It means that support of the test-function is sliced into by set of planes parallel to XY and in each 
plane ray transforms are computed. Hence, the data is parametrized by the triple (z, s, $\phi$)

For each such plane ray transforms are computed for rays with uniformly distributed 
directions and shifts (in the plane). The output is stored in a CSV file, where the data is ordered as follows:  
 > [z coordinate of slice], [shift in the plane], [polar angle in the plane], [value of the ray transform]

## Requirements 

The programs here are designed to work under Unix operating systems.  
To compile the project on your computer you need to have installed:  

GCC compiler, OpenMP libraries, GNU GSL libraries (+2.5)

## Compilation / Installation
  
  1. Run Makefile
      ```
        make install
      ```
  2. Clean directory from object files (optional):
  
      ```
        make clean 
      ```
  If you want to generate ray data for other test-function then you have to change the file
  'test_function.c' and repeat steps (1,2) to recompile the binary.
  
## Usage / Examples

(binary) -h --help &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;: Display usage information  
(binary) -p --parameters filename &nbsp;&nbsp;&nbsp;&nbsp;: Read parameters of the ray grid from config file  
(binary) -o --output filename &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;: Write output data to a file  
(binary) -n --nthreads number &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;: Use (number) of OpenMP threads for parallelization  

### Config file requirements

The purpose of the config file for parameters -p (--parameters) is to provide to the program information about the grid  
in ray space in slice-by-slice sampling scheme. The output data is parametrized by plane slice (z coordinate) and 
by the cooridnates of a ray in this plane (polar angle, shift). Also, though the test-function is
given by an analytical expression its integration over planes necessarily requires discretization which is also must be 
given in the config file. More precisely, a test-function is always assumed to be supported in a unit ball in 3D which 
lies inside the unit cube [-1,1]^3. So the parameter to be specified is a number of points per dimension in the uniform rectangular grid inside this unit cube.

1) The first line contains a number of shifts which are positioned uniformly along [-1,1] (for a fixed direction in 
a fixed plane).  
1) The second line contains a number of polar angles which are positioned uniformly along [0,2pi].  
3) The third line contains a number of z-slices which are positioned uniformly along [-1,1]. This value 
is also corresponds to a number of points per dimension on [-1,1] for the grid on the unit cube [-1,1]^3.

Some comments are allowed after each line, however the length of each line should not exceed 128 symbols.

### Example of a config file

> 129			: number of steps per fixed direction  
> 256			: number of directions, i.e., longitude angles  
> 129			: number of z-slices and points on the grid per dimension  


### Output

The output stored in a specified output file (-o --output) in a CSV format in the following order:  
**[z coordinate of slice], [shift in the plane], [polar angle in the plane], [value of the ray transform]**  

Example:  
> -1.000000, -1.000000, 5.571418, 0.000000  
> -1.000000, -1.000000, 5.595962, 0.000000  
> -1.000000, -1.000000, 5.620506, 0.000000  

### Examples of test-functions

One can find in folder 'test_functions' some examples of realizations of test-functions.  
In order to try them, rename any of those files to 'test_function.c' and copy them to the main directory of 'ray_analytic'.  
Then proceed with steps in 'Compilation / Installation' in order to obtain a compiled binary for a given
test-function. 


