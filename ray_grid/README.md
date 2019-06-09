## General description 

This program computes ray transforms in 3D of a test-function which is given by its values in a separate file.  

**IMPORTANT:** This program computes ray transforms not for all rays in 3D, but for the slice-by-slice scheme. 
It means that support of the test-function is sliced into by set of planes parallel to XY and in each 
plane ray transforms are computed. Hence, the rays (and the data) are parametrized by the triple (z, s, φ), where z 
is the coordinate of the slice (z=const), s - distance from the origin to the line, φ - angle for normal to the line. 
(z, s, φ) are uniformly distributed in [-1,1]x[-1,1]x(0,2π), where the number of points is defined by the user in 
the config file (binary/config.txt).


**TEST-FUNCTION INPUT:** 
The test-function is given by its values on the uniform discrete grid in [-1,1]^3. The size of the grid is set 
in the config file (binary/config.txt) in variable 'ngrid'. The values of the test-function must be given in a 
separate file in the following order:

	for ix = (0: ngrid-1)
	  for iy = (0: ngrid-1) 
	    for iz = (0: ngrid-1)
	        fprintf(test-function_file, "%lf\n", test_function(ix, iy, iz));
             endfor
	  endfor
	endfor


**OUTPUT:**
The ray transforms are returned in a separate file in the following format: 

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

1. The first line contains the number of points in the uniform grid in [-1,1]^3 in one dimension (endpoints are included).
2. The second line contains the number of z-slices which are positioned uniformly along [-1,1] (including endpoints -1,1).
2. The third line contains the number of shifts which are positioned uniformly along [-1,1] (including endpoints -1,1).  
3. The fourth line contains the number of polar angles which are positioned uniformly along [0,2π].

The length of each line should not exceed 128 symbols.

### Example of a config file

The config file (config.txt) is already placed in binary folder. It has the following format: 

> 129			: number of steps in the uniform grid in [-1,1]^3 in one dimension (endpoints are included)
> 129			: number of steps per fixed direction  
> 256			: number of directions, i.e., longitude angles  
> 129			: number of z-slices and points on the grid per dimension  


### Output

The output is stored in an output file (-o, --output [file] options) in a CSV format in the following format: 

 **[z], [s], [φ], [Pf(z, s, φ)]**  

Example of output:  

> -1.000000, -1.000000, 5.571418, 0.000000  
> -1.000000, -1.000000, 5.595962, 0.000000  
> -1.000000, -1.000000, 5.620506, 0.000000  
