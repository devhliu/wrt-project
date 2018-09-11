#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <gsl/gsl_integration.h>

#define MAX_COMMAND_SIZE 64

/* parameters utilities */ 

int init_parameters_alloc(int** parameters) {
  *parameters = (int*)malloc(4 * sizeof(int));
  
  return 0;
}
int init_parameters(char* filename, int* parameters) {
  
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
     return -1;
  } else if (count != 2) {
         fprintf(stderr, "Error: fscanf parameters 0 : matching failure, fscanf expected 1 integer value.\n");
    return -1;
  }
  count = fscanf(f, "%d %[^\n]\n", &(parameters[1]), buffer);
  if (count == EOF) {
     if (ferror(f)) {
        perror("Init parameters : fscanf parameters 1 : EOF failure\n");
     }  else {
        fprintf(stderr, "Error: fscanf parameters 1 : reached end of file before it was expected.\n");
     }
     return -1;
  } else if (count != 2) {
         fprintf(stderr, "Error: fscanf parameters 1 : matching failure, fscanf expected 1 integer value.\n");
    return -1;
  }
  count = fscanf(f, "%d %[^\n]\n", &(parameters[2]), buffer);
  if (count == EOF) {
     if (ferror(f)) {
        perror("Init parameters : fscanf parameters 2 : EOF failure\n");
     }  else {
        fprintf(stderr, "Error: fscanf parameters 2 : reached end of file before it was expected.\n");
     }
     return -1;
  } else if (count != 2) {
         fprintf(stderr, "Error: fscanf parameters 2 : matching failure, fscanf expected 1 integer value.\n");
    return -1;
  }
  count = fscanf(f, "%d %[^\n]\n", &(parameters[3]), buffer);
  if (count == EOF) {
     if (ferror(f)) {
        perror("Init parameters : fscanf parameters 3 : EOF failure\n");
     }  else {
        fprintf(stderr, "Error: fscanf parameters 3 : reached end of file before it was expected.\n");
     }
     return -1;
  } else if (count != 2) {
         fprintf(stderr, "Error: fscanf parameters 3 : matching failure, fscanf expected 1 integer value.\n");
    return -1;
  }
  
  fclose(f);
  return 0; //exit sucess
  
}
int clean_parameters(int* parameters) {
  free(parameters);
  
  return 0;
}

/* read ray transforms */ 

int init_ray_grid_alloc(int nshift, int nphi, int ngrid, double**** values) {
  int i, j;
  
  *values = (double***)malloc(sizeof(double**) * ngrid);
  for (i = 0; i < ngrid; ++i) {
      (*values)[i] = (double**)malloc(sizeof(double*) * nphi);
  }
  
  for (i = 0; i < ngrid; ++i) {
      for (j = 0; j < nphi; ++j) {
          (*values)[i][j] = (double*)malloc(sizeof(double) * nshift);
      }
  }
  
  return 0;
}
int init_ray_values(char* filename, int nshift, int nphi, int ngrid, double*** values) {
  
  FILE *f;
  f = fopen(filename, "r");
  
  if (f == NULL) {
    fprintf(stderr, "Recieved call %s while opening %s. \n", strerror(errno), filename);
    return errno;
  }
  
  int i_shift, i_phi, i_ngrid;
  int count; 
  
  for (i_ngrid = 0; i_ngrid < ngrid; ++i_ngrid) {
      for (i_shift = 0; i_shift < nshift; ++i_shift) {
          for (i_phi = 0; i_phi < nphi; ++i_phi) {
	      //read file line-by-line 
	      double tmp_zslice, tmp_shift, tmp_phi;
	      count = fscanf(f, "%lf, %lf, %lf, %lf\n", &tmp_zslice, &tmp_shift, &tmp_phi, &(values[i_ngrid][i_phi][i_shift]));
	      if (count == EOF) {
	         if (ferror(f)) {
	         perror("Init test-function values : fscanf read value : EOF failure\n");
	      } else {
	         fprintf(stderr, "Error : init test-function : fscanf : reached the end of file before it was expected.\n");
	      }
	      return -1;
	      } else if (count != 4) {
	             fprintf(stderr, "Error: init test-function : fscanf : matching failure, fscanf expected 4 floating-point values.\n");
	             return -1;
	      }
	      
          }
      }
  }
      
  fclose(f);  
  return 0;
}
int clean_ray_values(int nshift, int nphi, int ngrid, double*** values) {
  int i,j;
  for (i = 0; i < ngrid; ++i) {      // 
      for (j = 0; j < nphi; ++j) {   //  
          free(values[i][j]);        // clean directions phi 
      }
  }
  for (i = 0; i < ngrid; ++i) {    //
      free(values[i]);             // clean shifts
  }
  free(values);                    // clean z -- cross sections
  
  return 0;
}


/* initialize/clean Radon grid */

int init_radon_grid_alloc(int nshift,  int nphi, int ntheta, double** shift, double** phi, double** theta) {
    *shift = (double*)malloc(nshift * sizeof(double));
    *phi =   (double*)malloc(nphi   * sizeof(double));
    *theta = (double*)malloc(ntheta * sizeof(double));

    return 0;
}
int init_radon_grid(int nshift, int nphi, int ntheta,  double* shift, double* phi, double* theta) {
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
int clean_radon_grid(double* shift, double* phi, double* theta) {
   if (shift != NULL)
     free(shift);
   if (phi != NULL)
     free(phi);
   if (theta != NULL)
     free(theta);
   
   return 0;
}

/* filenames utilities */

int init_chunks_filenames_alloc(char* filename, int nchunks, char*** chunks_filenames, int size) {
  (*chunks_filenames) = (char**)malloc(sizeof(char*) * nchunks);
  int i;
  for (i = 0; i < nchunks; ++i) {
      (*chunks_filenames)[i] = (char*)malloc(sizeof(char) * size);
  }
  
  return 0;
}
int init_chunks_filenames(char* filename, int nchunks, char** chunks_filenames) {
  int i;
  for (i = 0; i < nchunks; ++i) {
    sprintf(chunks_filenames[i], "%d_%s", i, filename);
  }
  return 0;
}
int clean_chunks_files(int nchunks, char** chunks_filenames) {
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
int clean_chunks_filenames(int nchunks, char** chunks_filenames) {
  int i;
  for (i = 0; i < nchunks; ++i) {
    free(chunks_filenames[i]); 
  }
  free(chunks_filenames);
  return 0;
}
