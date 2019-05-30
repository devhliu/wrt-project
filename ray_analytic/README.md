
## General description 

This program computes ray transforms in 3D of a function which realized by an analytic expression in a separate file.  
File "test_function.c" contains a template of such realization which should be completed by the user before to compile.

**IMPORTANT:** This program computes ray transforms not for all rays in 3D, but for the slice-by-slice scheme. 
It means that support of the test-function is sliced into by set of planes parallel to XY and in each 
plane ray transforms are computed. Hence, the rays (and the data) are parametrized by the triple (z, s, φ), where z is the coordinate of the slice (z=const), s - distance from the origin to the line, φ - angle for normal to the line. 
(z, s, φ) are uniformly distributed in [-1,1]x[-1,1]x(0,2π), where the number of points is defined by the user in 
the config file (binary/config.txt).

**FORMAT** The ray transforms are returned in a separate file, in the following format: 

      for islice = (0: nslice-1)
        for ishift = (0: nshift-1)
          for iphi = (0: nhpi-1) 
             foutput(output_file, "%lf, %lf, %lf, %lf\n", z[islice], s[ishift], φ[iphi], Pf(z[islice], s[ishift], φ[iphi]);
          endfor
        endfor
      endfor


## Requirements 

The programs here are designed to work under Unix operating systems. To compile the project on your 
computer you need to have installed:  

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


### Config file 

The purpose of the config file for parameters (-p, --parameters [file] options) is to provide to the program information 
on the grid of rays. 

1. The first line contains a number of shifts which are positioned uniformly along [-1,1] (including endpoints -1,1).  
2. The second line contains a number of polar angles which are positioned uniformly along [0,2π].
3. The third line contains a number of z-slices which are positioned uniformly along [-1,1] (including endpoints -1,1).

The length of each line should not exceed 128 symbols.

### Example of a config file

A config (config.txt) file is already placed in binary folder. It has the following format: 

> 129			: number of steps per fixed direction  
> 256			: number of directions, i.e., longitude angles  
> 129			: number of z-slices and points on the grid per dimension  


### Output

The output is stored in an output file (-o, --output [file] options) in a CSV format in the following order: 
 **[z], [s], [φ], [Pf(z, s, φ)]**  

Example of output:  

> -1.000000, -1.000000, 5.571418, 0.000000  
> -1.000000, -1.000000, 5.595962, 0.000000  
> -1.000000, -1.000000, 5.620506, 0.000000  

### Examples of test-functions

One can find in folder 'examples-test-functions' some realizations of test-functions.  
In order to try them, rename any of those files to 'test_function.c' and copy them to the main directory.  
Then proceed with steps in 'Compilation / Installation' in order to obtain a compiled binary for a given
test-function. 



