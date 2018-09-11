#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <assert.h>
#include <omp.h>
#include <sys/time.h> //measure physical time of evaluation
#include "init.h"
#include "integral.h"

#define FILENAME_MAXSIZE 64
/*
 * The program reads the grid data for the test-function from file (-i filename)
 * and computes its ray transforms in slices of planes (parallel to XY). The set of rays 
 * in each slicing plane corresponds to 'uniform geometry' (uniform shifts x uniform angles on circle)
 * The parameters of the input test function an parameters of ray transforms are given in a 
 * separate configuration file (-p filename).
 * 
 * IMPORTANT : It is always assumed that the test-function is compactly supported with support in the
 * centered unit ball. Otherwise the generated data will be incomplete and this will produce artifacts in reconstructions.
 *
 * The results of computations are stored in a separate file (-o filename). 
 * 
 * Output is stored either in CSV.
*/

/* ***************************************************************************
 * ****************************** HANDLING FILES *****************************
 * ***************************************************************************
 */
void data_generate_readme(char* filename, int nz, int nshift, int nphi) {
  FILE *file;
  char readme_filename [ FILENAME_MAXSIZE ];
  
  sprintf(readme_filename, "readme_%s", filename);

  file = fopen(readme_filename, "w");
  fprintf(file, "  Ray transforms in 3D for test-function given in %s\n", filename);
  fprintf(file, "  Format of the data: [z], [phi], [shift], [Pf]\n");
  fprintf(file, "  Parameters of the data:\n");
  fprintf(file, "    number of layers OZ direction  [nz]   : %d,\n", nz);
  fprintf(file, "    number of angles on circle     [phi]  : %d,\n", nphi);
  fprintf(file, "    number of shifts per direciton [s]    : %d.\n", nshift);
    
  fclose(file);
}

void data_aggregate_chunks(char* output_filename, char ** chunks_filenames, int nchunks) {
  
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
    //printf("  Data from %s has been sucessfully transmitted to %s\n", chunks_filenames[i_file], output_filename);
 }
}

/* ***************************************************************************
 * ****************************** HANDLING INPUT *****************************
 * ***************************************************************************
 */


const char* program_name;

/* Prints usage information to STREAM (stdout or stderr), and 
 * exits the program with EXIT_CODE. Does not return. */

void print_usage(FILE* stream, int exit_code) {
  fprintf(stream, "Usage: %s -p (filename) -i (filename) -o (filename) -n (number)\n", program_name);
  fprintf(stream, 
	  "   -h --help                 Display usage information.\n"
	  "   -p --parameters filename  Read parameters of the test-function and 3D ray transforms from file.\n"
	  "   -i --input filename       Read test-functiom data from file (according to parameters in '-p').\n"
	  "   -o --output filename      Place output data to a file.\n"
	  "   -n --nthreads number      Number of OpenMP threads for parallelization of calculations.\n");
  exit(exit_code);
}

/* ***************************************************************************
 * ****************************** ENTRY POINT ********************************
 * ***************************************************************************
 */

int main(int argc, char * argv[]) {
  
  if (argc == 1) {
     fprintf(stderr, "%s: arguments required. Try '-h' or '--help'.\n", argv[0]);
     exit(EXIT_FAILURE);
  }
  
  int next_option;
  const char* short_options = "hp:i:o:n:";
  const struct option long_options[] = {
    {"help", no_argument, 0, 0}, 
    {"parameters", required_argument, 0, 0},
    {"input", required_argument, 0, 0}, 
    {"output", required_argument, 0, 0},
    {"nthreads", required_argument, 0, 0},
    {NULL, 0, NULL, 0}				/* Required at the end of array */
  };

  char* parameters_filename;
  char* input_filename; 
  char* output_filename;
  int nthreads = 1;
  
  program_name = argv[0];

  /* Reading long options */ 
  do {
    next_option = getopt_long(argc, argv, short_options, long_options, NULL);
    
    switch (next_option) {
      case 'h':
	print_usage(stdout, 0);
      case 'p':
	parameters_filename = optarg;
	assert(optarg[0] != '-');
	printf("Parameters file: %s\n", optarg);
	break;
      case 'i':
	input_filename = optarg;
	assert(optarg[0] != '-');
	printf("Input file with test-function data: %s\n", optarg);
	break;
      case 'o':
	output_filename = optarg;
	assert(optarg[0] != '-');
	printf("Output file with ray data: %s\n", optarg);
	break;
      case 'n':
	nthreads = atoi(optarg);
	assert(nthreads > 0);
	printf("Number of threads: %d\n", nthreads);
	break;
      case ':': /* missing argument */ 
	fprintf(stderr, "%s: option '-%c' requires an argument\n",
		argv[0], optopt);
	print_usage(stderr, 1);
      case '?':
	print_usage(stderr, 1);
      case -1: /* Done with the options. */
	break;
      default:
	abort();
    }
  } while (next_option != -1);
  
  if (optind < argc) {
    printf("Warning: non-option ARGV-arguments are omitted. See '-h' or '--help'.\n"); 
  }
  
  /**************************************** START INITIALIZATIONS *****************************************/
  /* Init configuration parameters of the grid and of Radon data from file */
  printf("  Reading parameters from %s...", parameters_filename);
  int *parameters;
  if (init_parameters_alloc(&parameters) != 0) {
     perror("Aborting : init_parameters_alloc.\n");
     exit(EXIT_FAILURE);
  }
  
  if (init_parameters(parameters_filename, parameters) != 0) {
    perror("Aborting : init_parameters.\n");
    exit(EXIT_FAILURE);
  }
  printf("Done.\n");
  
  int ngrid_tfunc = parameters[0];	// number of points per dimension (in [-1,1]) in the grid for the test-function
  int nshift      = parameters[1];
  int nphi        = parameters[2];
  int ngrid       = parameters[3];	// number of z-slices; (affects density of points for and integral)
  
  
  /* Init values of the test function on the grid */
  printf("  Initializing test-function data...");
  double*** values;
  if (init_values_alloc(ngrid_tfunc, &values) != 0) {
     perror("Aborting : init_values_alloc.\n");
     exit(EXIT_FAILURE);
  }
  if (init_values(input_filename, ngrid_tfunc, values) != 0) {
     perror("Aborting : init_values.\n");
     exit(EXIT_FAILURE);
  }
  printf("Done.\n");
  
  
  /* Init grid in Radon space */
  double *phi, *shift;
  printf("  Initializing the grid in ray transforms space...");
  if (init_ray_grid_alloc(nshift, nphi, &shift, &phi) != 0) {
     perror("Aborting : init_ray_grid_alloc.\n");
     exit(EXIT_FAILURE);
  }
  if (init_ray_grid(nshift, nphi, shift, phi) != 0 ) {
     perror("Aborting : init_ray_grid.\n");
     exit(EXIT_FAILURE);
  }
  printf("Done.\n\n");
  
  /************************************** END OF INITIALIZATIONS ******************************************/
  
  
  /************************* START RADON TRANSFORM WITH PARALLELIZATION VIA OPENMP ************************/
  //start timer
  struct timeval start, end;
  gettimeofday(&start, NULL);  

  
  /* Init chunks filenames (nchunks = nthreads). Chunk filename %d_(output_filename) */
  char** chunks_filenames;
  init_chunks_filenames_alloc(output_filename, nthreads, &chunks_filenames, FILENAME_MAXSIZE);
  init_chunks_filenames(output_filename, nthreads, chunks_filenames);
  
  omp_set_num_threads(nthreads);
  #pragma omp parallel
  {
        int thread_ID = omp_get_thread_num();

        //set equally parts of shifts for each thread
	int block_size = ngrid_tfunc / nthreads + (thread_ID < (ngrid_tfunc % nthreads) ? 1 : 0);
	
        int i_slice_min = (thread_ID < (ngrid_tfunc % nthreads) ? (ngrid_tfunc / nthreads + 1) * thread_ID : (ngrid_tfunc / nthreads)*thread_ID + (ngrid_tfunc % nthreads)),
	    i_slice_max = i_slice_min + block_size;
	    
        printf("    Thread %d. Domain of work (z slices): start %d, end %d, size %d\n", thread_ID, i_slice_min, i_slice_max-1, block_size);
        
	
        //open chunk ouptut file
        FILE *foutput;
        foutput = fopen(chunks_filenames[thread_ID], "w");
	int i_zslice, i_phi, i_shift;
	
	double zslice;
        const double dzslice = 2.0 / (ngrid_tfunc - 1);
        //iterate ray integrals over z-slices
        for (i_zslice = i_slice_min; i_zslice < i_slice_max; ++i_zslice) {      
	    zslice = (-1.0) + i_zslice * dzslice;
	    
	    for (i_shift = 0; i_shift < nshift; ++i_shift) {
	     
                //iterate ray integrals over longitude angle
                for (i_phi = 0; i_phi < nphi; ++i_phi) {
	
	            //ray integral
		    double rt = ray_transform(values, ngrid_tfunc, 
						 zslice, shift[i_shift], phi[i_phi], ngrid);
		    fprintf(foutput, "%lf, %lf, %lf, %lf\n", zslice, shift[i_shift], phi[i_phi], rt);
                }
            }
        }
        
	fclose(foutput);
	printf("    Thread %d. Job is done.\n", thread_ID);
  } //threads finish
  
  //stop timer
  gettimeofday(&end, NULL);
  double delta_t = ((end.tv_sec  - start.tv_sec) * 1000000u + 
         end.tv_usec - start.tv_usec) / 1.e6;
  printf("\n");
  printf("  Number of threads: %d, elapsed time: %3.1f sec\n", nthreads, delta_t);
  
  /************************* END RADON TRANSFORM WITH PARALLELIZATION VIA OPENMP **************************/
  
  
  /*********************************** START AGGREGATING DATA TO ONE FILE *********************************/
  
  printf("  Generating readme...");
  data_generate_readme(output_filename, ngrid, nshift, nphi);
  printf("  Done.\n");
  printf("  Aggregating chunks to one output file...");
  data_aggregate_chunks(output_filename, chunks_filenames, nthreads);
  printf("Done.\n");
  
  /*********************************** END AGGREGATING DATA TO ONE FILE ***********************************/
  

  /******************************************* START CLEANING MEMORY **************************************/
  
  printf("  Cleaning allocated memory...");
  
  // clean files
  if (clean_chunks_files(nthreads, chunks_filenames) != 0) {
     fprintf(stderr, "Warning : clean_chunks_files.\n");
  }
  if (clean_chunks_filenames(nthreads, chunks_filenames) != 0) {
     fprintf(stderr, "Warning : clean_chunks_filenames.\n");
  }
  
  // clean memory
  if (clean_parameters(parameters) != 0) {
     fprintf(stderr, "Warning : clean_parameters.\n");
  }
  if (clean_values(ngrid_tfunc, values) != 0) {		//clean data for the test-function
     fprintf(stderr, "Warning : clean_values.\n");
  }
  if (clean_ray_grid(shift, phi) != 0) {
     fprintf(stderr, "Warning : clean_ray_grid.\n");
  }
  
  printf("Done.\n");
  
  /******************************************* END CLEANING MEMORY ****************************************/
  
  exit(EXIT_SUCCESS);
}
