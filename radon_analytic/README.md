## General information

This program performs analytic evaluations of classical Radon transforms in 3D (along 2D planes) of analytic functions whose prototype should be realized in "test_function.c".  

File "test_function.c" contains a template of such realization. Note that ony the function with name "test_function" will be used for computations. After realizaiton of your test function you have to compile the code so it can be used.  

## Requirements 

The programs here are designed to work under Unix operating systems.  
To compile the project on your computer you need to have installed:  

GCC compiler, OpenMP libraries, GNU GSL libraries

## Compilation / Installation
  1) Go to 'src' directory:  
        ```
          cd src
        ```
  2) Open Makefile and set the name of the output file:
        open Makefile in any text editor and set
        ```
          ONAME=(output name of your binary)
        ```
  
  3) Run Makefile
      ```
        make install
      ```
  4) Clean directory from object files (optional):
  
      ```
        make clean 
      ```
  If you want to generate data for other test-function then you have to change the file
  'test_function.c' and repeat steps (1-4), possibly setting a new name in 'ONAME'
  
## Usage / Examples


