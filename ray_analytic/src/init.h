#ifndef _H_RAY_INIT
#define _H_RAY_INIT

/* parameters utilites */

int init_parameters_alloc(int** parameters);
int init_parameters(char* filename, int* parameters);
int clean_parameters(int* parameters);

/* ray utilites */

int init_ray_grid_alloc(int nshift, int nphi, double** shift, double** phi);
  /*
   * shift -- uniform points on [-1, 1]
   * phi   -- uniform grid on S^1 = [0,2pi]
   */
int init_ray_grid(int nshift, int nphi, double* shift, double* phi);
int clean_ray_grid(double* shift, double* phi);


/* Threads filenames utilites */

int init_chunks_filenames_alloc(char* filename, int nchunks, char*** nchunks_filenames, int size);
int init_chunks_filenames(char* filename, int nchunks, char** nchunks_filenames);
int clean_chunks_files(int nthreads, char** chunks_filenames);
int clean_chunks_filenames(int nchunks, char** nchunks_filenames);

#endif