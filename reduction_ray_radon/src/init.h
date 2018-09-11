#ifndef _H_REDUCTION_INIT
#define _H_REDUCTION_INIT

/* parameters utilites */

int init_parameters_alloc(int** parameters);
int init_parameters(char* filename, int* parameters);
int clean_parameters(int* parameters);

/* ray transforms utilities */

int init_ray_grid_alloc(int nshift, int nphi, int ngrid, double**** values);		// allocates memory for ray transforms
int init_ray_values(char* filename, int nshift, int nphi, int ngrid, double*** values); // reads ray transforms from file into array values[ngrid][shift][phi]
int clean_ray_values(int nshift, int nphi, int ngrid, double*** values);                // cleans allocated memory


int init_radon_grid_alloc(int nshift, int nphi, int ntheta, double** shift, double **phi, double **theta);
  /*
   * shift -- uniform grid   on [-1, 1]
   * theta -- Gauss   grid   on [0, PI]
   * phi   -- uniform grid   on S^1 = [0,2*PI]
   * 
   * depends on GSL library
   */
int init_radon_grid(int nshift, int nphi, int ntheta, double* shift, double* phi, double* theta);
int clean_radon_grid(double* shift, double* phi, double* theta);


/* Threads filenames utilites */

int init_chunks_filenames_alloc(char* filename, int nchunks, char*** nchunks_filenames, int size);
int init_chunks_filenames(char* filename, int nchunks, char** nchunks_filenames);
int clean_chunks_files(int nthreads, char** chunks_filenames);
int clean_chunks_filenames(int nchunks, char** nchunks_filenames);


#endif 