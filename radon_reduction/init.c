#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <gsl/gsl_integration.h>

#define MAX_COMMAND_SIZE 128

/* parameters utilities */
// ---------------------------------------------------------------------------------------//

int parameters_alloc(int** parameters) {
    *parameters = (int*)malloc(4 * sizeof(int));
    
    return 0;
}

int parameters_init(char* filename, int* parameters) {
    
    FILE *f;
    char buffer[128];
    
    f = fopen(filename, "r");
    if (f == NULL) {
        fprintf(stderr, "Recieved call %s while opening %s. \n", strerror(errno), filename);
        return errno;
    }
    
    int count;
    count = fscanf(f, "%d %[^\n]\n", &(parameters[0]), buffer);
    if (count == EOF) {
        if (ferror(f)) {
            perror("Init parameters : fscanf parameters 0 : EOF failure\n");
        }  else {
            fprintf(stderr, "Error: fscanf parameters 0 : reached end of file before it was expected.\n");
        }
	fclose(f);
        return -1;
    } else if (count != 2) {
        fprintf(stderr, "Error: fscanf parameters 0 : matching failure, fscanf expected 1 integer value.\n");
	fclose(f);
        return -1;
    }
    count = fscanf(f, "%d %[^\n]\n", &(parameters[1]), buffer);
    if (count == EOF) {
        if (ferror(f)) {
            perror("Init parameters : fscanf parameters 1 : EOF failure\n");
        }  else {
            fprintf(stderr, "Error: fscanf parameters 1 : reached end of file before it was expected.\n");
        }
	fclose(f);
        return -1;
    } else if (count != 2) {
        fprintf(stderr, "Error: fscanf parameters 1 : matching failure, fscanf expected 1 integer value.\n");
	fclose(f);
        return -1;
    }
    count = fscanf(f, "%d %[^\n]\n", &(parameters[2]), buffer);
    if (count == EOF) {
        if (ferror(f)) {
            perror("Init parameters : fscanf parameters 2 : EOF failure\n");
        }  else {
            fprintf(stderr, "Error: fscanf parameters 2 : reached end of file before it was expected.\n");
        }
	fclose(f);
        return -1;
    } else if (count != 2) {
        fprintf(stderr, "Error: fscanf parameters 2 : matching failure, fscanf expected 1 integer value.\n");
	fclose(f);
        return -1;
    }
    count = fscanf(f, "%d %[^\n]\n", &(parameters[3]), buffer);
    if (count == EOF) {
        if (ferror(f)) {
            perror("Init parameters : fscanf parameters 3 : EOF failure\n");
        }  else {
            fprintf(stderr, "Error: fscanf parameters 3 : reached end of file before it was expected.\n");
        }
	fclose(f);
        return -1;
    } else if (count != 2) {
        fprintf(stderr, "Error: fscanf parameters 3 : matching failure, fscanf expected 1 integer value.\n");
	fclose(f);
        return -1;
    }
    
    fclose(f);
    return 0; //exit sucess
    
}
int parameters_clean(int* parameters) {
    free(parameters);
    
    return 0;
}

/* read ray transforms */ 
// ---------------------------------------------------------------------------------------//
int ray_values_alloc(int nshift, int nphi, int ngrid, double**** ray_values) {
    int i, j;
    
    *ray_values = (double***)malloc(sizeof(double**) * ngrid);
    for (i = 0; i < ngrid; ++i) {
        (*ray_values)[i] = (double**)malloc(sizeof(double*) * nphi);
    }
    
    for (i = 0; i < ngrid; ++i) {
        for (j = 0; j < nphi; ++j) {
            (*ray_values)[i][j] = (double*)malloc(sizeof(double) * nshift);
        }
    }
    
    return 0;
}
int ray_values_init(char* filename, int nshift, int nphi, int ngrid, double*** ray_values, int small_format) {
    
    FILE *f;
    f = fopen(filename, "r");
    if (f == NULL) {
        fprintf(stderr, "Recieved call %s while opening %s. \n", strerror(errno), filename);
        return errno;
    }
    int i_shift, i_phi, i_ngrid;
    int count; 
    
    int count_format;
    if (small_format) {
        count_format = 1;
    } else {
        count_format = 4;
    }

    for (i_ngrid = 0; i_ngrid < ngrid; ++i_ngrid) {
        for (i_shift = 0; i_shift < nshift; ++i_shift) {
            for (i_phi = 0; i_phi < nphi; ++i_phi) {
               
                //read file with ray data line-by-line 
                if (small_format) {
                    count = fscanf(f, "%lf\n", &(ray_values[i_ngrid][i_phi][i_shift]));
                } else {
                    double tmp_zslice, tmp_shift, tmp_phi;
                    count = fscanf(f, "%lf, %lf, %lf, %lf\n", &tmp_zslice, &tmp_shift, &tmp_phi, &(ray_values[i_ngrid][i_phi][i_shift]));
                }
                if (count == EOF) {
                    if (ferror(f)) {
                        perror("Init test-function values : fscanf read value : EOF failure\n");
                    } else {
                        fprintf(stderr, "Error : init test-function : fscanf : reached the end of file before it was expected.\n");
                    }
		    fclose(f);
                    return -1;
                } else if (count != count_format) {
                    fprintf(stderr, "Error: init test-function : fscanf : matching failure, fscanf expected %d floating-point values.\n", count_format);
		    fclose(f);
                    return -1;
                }    
            }
        }
    }
    
    fclose(f);  
    return 0;
}
int ray_values_clean(int nshift, int nphi, int ngrid, double*** ray_values) {
    int i,j;
    for (i = 0; i < ngrid; ++i) {      // 
        for (j = 0; j < nphi; ++j) {   //  
            free(ray_values[i][j]);        // clean directions phi 
        }
    }
    for (i = 0; i < ngrid; ++i) {    //
        free(ray_values[i]);         // clean shifts
    }
    free(ray_values);                    // clean z -- cross sections
    
    return 0;
}


/* initialize/clean Radon grid */
// ---------------------------------------------------------------------------------------//
int radon_grid_alloc(int nshift,  int nphi, int ntheta, double** shift, double** phi, double** theta) {
    *shift = (double*)malloc(nshift * sizeof(double));
    *phi =   (double*)malloc(nphi   * sizeof(double));
    *theta = (double*)malloc(ntheta * sizeof(double));
    
    return 0;
}
int radon_grid_init(int nshift, int nphi, int ntheta,  double* shift, double* phi, double* theta) {
    /*
     * phi --   uniform grid on S^1 = [0,2pi]
     * theta -- Gauss's angles on [0, pi]
     * shift -- uniform points on [-1, 1]
     */
    
    //init shift
    const double dshift = 2.0 / (nshift - 1);
    int i_shift; 
    for (i_shift = 0; i_shift < nshift; ++i_shift) {
        shift[i_shift] = (-1.0) + i_shift * dshift;
    }
    
    //init phi 
    const double dphi = (2 * M_PI) / nphi;
    int i_phi;
    for (i_phi = 0; i_phi < nphi; ++i_phi) {
        phi[i_phi] = dphi * i_phi;
    }
    
    //init theta (using GSL library to compute Gaussian points)
    int i_theta;
    gsl_integration_glfixed_table *gsl_quad_data = gsl_integration_glfixed_table_alloc(ntheta);
    for (i_theta = 0; i_theta < ntheta; ++i_theta) {
        double gauss_point, gauss_weight;
        gsl_integration_glfixed_point(-1.0, 1.0, i_theta, &gauss_point, &gauss_weight, gsl_quad_data);
        theta[i_theta] = acos(gauss_point);
    }
    
    gsl_integration_glfixed_table_free(gsl_quad_data);
    
    return 0;
}
int radon_grid_clean(double* shift, double* phi, double* theta) {
    if (shift != NULL)
        free(shift);
    if (phi != NULL)
        free(phi);
    if (theta != NULL)
        free(theta);
    
    return 0;
}

/* filenames utilities */
// ---------------------------------------------------------------------------------------//
int chunks_filenames_alloc(char* filename, int nchunks, char*** chunks_filenames, int size) {
    (*chunks_filenames) = (char**)malloc(sizeof(char*) * nchunks);
    int i;
    for (i = 0; i < nchunks; ++i) {
        (*chunks_filenames)[i] = (char*)malloc(sizeof(char) * size);
    }
    
    return 0;
}
int chunks_filenames_init(char* filename, int nchunks, char** chunks_filenames) {
    int i;
    for (i = 0; i < nchunks; ++i) {
        sprintf(chunks_filenames[i], "%d_%s", i, filename);
    }
    return 0;
}
int chunks_files_clean(int nchunks, char** chunks_filenames) {
    int i_file;
    for (i_file = 0; i_file < nchunks; ++i_file) {
        char buffer[MAX_COMMAND_SIZE];
        sprintf(buffer, "rm -f %s", chunks_filenames[i_file]);
        int sys_call = system(buffer);
        if (sys_call == -1) {
            fprintf(stderr, "Warning : clean chunks files : system call (rm -f ...) : failed to clean file %s.\n", chunks_filenames[i_file]); 
        }
    }
    return 0;
}
int chunks_filenames_clean(int nchunks, char** chunks_filenames) {
    int i;
    for (i = 0; i < nchunks; ++i) {
        free(chunks_filenames[i]); 
    }
    free(chunks_filenames);
    return 0;
}


void chunks_aggregate(char* output_filename, char ** chunks_filenames, int nchunks) {
  
  FILE *out; //output file
  out = fopen(output_filename, "w");
  
  int i_file;
  for (i_file = 0; i_file < nchunks; ++i_file) {
    //write data from local files (*fp) to output file (*out)
    FILE *fp;
    fp = fopen(chunks_filenames[i_file], "r");
    
    //transmit by chunks of 1 mb
    char buffer[1024 * 1024];
    int size;
    do {
      size = fread(buffer, 1, sizeof(buffer), fp);
      if (size <= 0) break;
      fwrite(buffer, 1, size, out);
    } while (size == sizeof(buffer));
    //reached EOF, close local file
    fclose(fp);
 }

 fclose(out);
}

// --------------------------------- info utilities -----------------------------------------------//

void print_usage(FILE* stream, char* program_name, int exit_code) {
  fprintf(stream, "Usage: %s -p (filename) -o (filename) -n (number)\n", program_name);
  fprintf(stream, 
      "   -h --help                 Display this usage information.\n"
      "   -p --parameters filename  Read parameters of the input/output sampling grid for ray/Radon transforms.\n"
      "   -i --input filename       Read data of ray transforms from file.\n"
      "   -s --format               Choose small format of data in input file.\n"
      "   -o --output filename      Write output data to the file.\n"
      "   -n --nthreads number      Number of OpenMP threads for parallelization.\n");
  exit(exit_code);
}


void generate_readme(char* filename, int nphi, int ntheta, int nshift) {
  FILE *doc_file;
  char doc_filename [ FILENAME_MAX ];
  
  sprintf(doc_filename, "readme_%s", filename);

  doc_file = fopen(doc_filename, "w");
  fprintf(doc_file, "  Reduction of (weighted) ray transforms (%s) to (weighted) Radon transforms.\n", filename);
  fprintf(doc_file, "  Format of the data: [s], [phi], [theta], [Rf]\n");
  fprintf(doc_file, "  Parameters of the data:\n");
  fprintf(doc_file, "      number of shifts [s]    : %d.\n", nshift);
  fprintf(doc_file, "      number of angles [phi]  : %d,\n", nphi);
  fprintf(doc_file, "      number of angles [theta]: %d,\n", ntheta);
    
  fclose(doc_file);
}
