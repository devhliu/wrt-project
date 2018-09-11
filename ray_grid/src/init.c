#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>

#define MAX_COMMAND_SIZE 64

/* Grid parameters */ 

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
        fprintf(stderr, "Error: fscanf parameters 0 : reached end of file before it was expected. \n");
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
        fprintf(stderr, "Error: fscanf parameters 1 : reached end of file before it was expected. \n");
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
        fprintf(stderr, "Error: fscanf parameters 2 : reached end of file before it was expected. \n");
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
        fprintf(stderr, "Error: fscanf parameters 3 : reached end of file before it was expected. \n");
     }
     return -1;
  } else if (count != 2) {
         fprintf(stderr, "Error: fscanf parameters 3 : matching failure, fscanf expected 1 integer value.\n");
         return -1;
  }
  
  fclose(f);
  return 0;	//exit sucess
  
}

int clean_parameters(int* parameters) {
  free(parameters);
  return 0;
}

/* Initial values of the test-function */

int init_values_alloc(int ngrid_tfunc, double**** values) {
  int i,j;
  
  *values = (double***)malloc(sizeof(double**) * ngrid_tfunc);
  for (i = 0; i < ngrid_tfunc; ++i) {
      (*values)[i] = (double**)malloc(sizeof(double*) * ngrid_tfunc);
  }
  
  for (i = 0; i < ngrid_tfunc; ++i) {
      for (j = 0; j < ngrid_tfunc; ++j) { 
          (*values)[i][j] = (double*)malloc(sizeof(double) * ngrid_tfunc);
      }
  }
  return 0;  
}

int init_values(char* filename, int ngrid_tfunc, double*** values) {
  FILE *f;
  
  f = fopen(filename, "r");
  
  if (f == NULL) {
     fprintf(stderr, "Recieved call %s while opening %s. \n", strerror(errno), filename);
     return errno;
  }
  
  int xi, yi, zi;
  int count;
  for (xi = 0; xi < ngrid_tfunc; ++xi)
    for (yi = 0; yi < ngrid_tfunc; ++yi)
      for (zi = 0; zi < ngrid_tfunc; ++zi) {
	
	//read file line-by-line 
	count = fscanf(f, "%lf\n", &(values[xi][yi][zi]));
	if (count == EOF) {
	   if (ferror(f)) {
	      perror("Init test-function values : fscanf read value : EOF failure\n");
	   } else {
	      fprintf(stderr, "Error : init test-function : fscanf : reached the end of file before it was expected.\n");
	   }
	   return -1;
	} else if (count != 1) {
	      fprintf(stderr, "Error: init test-function : fscanf : matching failure, fscanf expected 1 floating-point value.\n");
	      return -1;
	}
     }
      
  fclose(f);  
  return 0;	//exist sucess
}
int clean_values(int ngrid_tfunc, double*** values) {
  int i,j;
  for (i = 0; i < ngrid_tfunc; ++i) {
      for (j = 0; j < ngrid_tfunc; ++j) {
          free(values[i][j]);
      }
  }
  for (i = 0; i < ngrid_tfunc; ++i) {
      free(values[i]);
  }
  free(values);
  
  return 0;
}


/* Radon grid */

int init_ray_grid_alloc(int nshift, int nphi, double** shift, double** phi) {
  *shift = (double*)malloc(nshift * sizeof(double));  
  *phi = (double*)malloc(nphi * sizeof(double));    
  return 0;
}
int init_ray_grid(int nshift, int nphi, double* shift, double* phi) {
  
  /*
   * phi --   uniform grid on S^1 = [0,2pi]
   * theta -- Gauss's angles on [0, pi]
   * shift -- uniform points on [-1, 1]
   */
  
  //init shift
  const double dshift = 2.0 / (nshift - 1);
  int i_shift; 
  for (i_shift = 0; i_shift < nshift; ++i_shift)
    shift[i_shift] = (-1.0) + i_shift * dshift;
  
  //init phi 
  const double dphi = (2 * M_PI) / nphi;
  int i_phi;
  for (i_phi = 0; i_phi < nphi; ++i_phi)
    phi[i_phi] = dphi * i_phi;
  
  return 0; // code sucess
}
int clean_ray_grid(double* shift, double* phi) {
   
   if (shift != NULL)
     free(shift);
   
   if (phi != NULL)
     free(phi);
   
   return 0;
}

/* Threads filenames */

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
