#ifndef _H_REDUCTION_INIT
#define _H_REDUCTION_INIT

/* parameters utilites */
// ---------------------------------------------------------------------------------------//
int parameters_alloc(int** parameters);
int parameters_init(char* filename, int* parameters);
int parameters_clean(int* parameters);

/* ray transforms utilities */
// ---------------------------------------------------------------------------------------//
int ray_values_alloc(int nshift, int nphi, int ngrid, double**** ray_values);		            // allocate memory for ray transforms
int ray_values_init(char* filename, int nshift, int nphi, int ngrid, double*** ray_values, int small_format);     // read ray transforms from file into array ray_values[ngrid][shift][phi]
int ray_values_clean(int nshift, int nphi, int ngrid, double*** ray_values);                    // clean allocated memory


int radon_grid_alloc(int nshift, int nphi, int ntheta, double** shift, double **phi, double **theta);
int radon_grid_init(int nshift, int nphi, int ntheta, double* shift, double* phi, double* theta);
int radon_grid_clean(double* shift, double* phi, double* theta);


/* Threads filenames utilites */
// ---------------------------------------------------------------------------------------//
int  chunks_filenames_alloc(char* filename, int nchunks, char*** nchunks_filenames, int size);
int  chunks_filenames_init(char* filename, int nchunks, char** chunks_filenames);
int  chunks_files_clean(int nthreads, char** chunks_filenames);
int  chunks_filenames_clean(int nchunks, char** chunks_filenames);
void chunks_aggregate(char* output_filename, char ** chunks_filenames, int nchunks);

// ------------------------------ info utilities -----------------------------------------//

void print_usage(FILE* stream, char* program_name, int exit_code);
void generate_readme(char* filename, int nphi, int ntheta, int nshift);

// ---------------------------------------------------------------------------------------//

#endif 
