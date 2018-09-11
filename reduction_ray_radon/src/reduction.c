#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>   // input options handling
#include <assert.h>
#include <omp.h>      // parallelization in the set of 2d-planes
#include <sys/time.h> // measure of parallel time evaluation

#include "init.h"
#include "slicing.h"


#define FILENAME_MAXSIZE 64

/*
  * Program performes reduction of data given by (weighted) ray transforms in 3D 
  * to the data given by (weighted) Radon transforms along 2d-planes in 3D. Such 
  * reduction is a consequence of fibration of plane into parallel lines (even for
  * weighted ray-Radon transforms).
  * Parameters of ray transforms on the input and of Radon transforms on the output are 
  * stored in a configuration file (-p filename). 
  * 
  * IMPORTANT : It is always assumed that the test-function is compactly supported with support in the
  * centered unit ball. Otherwise the generated data will be incomplete and this will produce artifacts in 
  * reconstructions.
  *
  * Output is stored in CSV format.
  * 
  */
  
const char* program_name;

/* Prints usage information to STREAM (stdout or stderr), and 
 * exits the program with EXIT_CODE. */

void print_usage(FILE* stream, int exit_code) {
  fprintf(stream, "Usage: %s -p (filename) -o (filename) -n (number)\n", program_name);
  fprintf(stream, 
	  "   -h --help                 Display this usage information.\n"
	  "   -p --parameters filename  Read parameters of the input/output sampling grid for ray/Radon transforms.\n"
	  "   -i --input filename       Read data of ray transforms from file.\n"
	  "   -o --output filename      Write output data to the file.\n"
	  "   -n --nthreads number      Number of OpenMP threads for parallelization.\n");
  exit(exit_code);
}


void generate_data_readme(char* filename, int nphi, int ntheta, int nshift) {
  FILE *doc_file;
  char doc_filename [ FILENAME_MAXSIZE ];
  
  sprintf(doc_filename, "readme_%s", filename);

  doc_file = fopen(doc_filename, "w");
  fprintf(doc_file, "  Radon transforms from ray transforms given in %s\n", filename);
  fprintf(doc_file, "  Format of the data: [s], [phi], [theta], [Rf]\n");
  fprintf(doc_file, "  Parameters of the data:\n");
  fprintf(doc_file, "    number of shifts [s]    : %d.\n", nshift);
  fprintf(doc_file, "    number of angles [phi]  : %d,\n", nphi);
  fprintf(doc_file, "    number of angles [theta]: %d,\n", ntheta);
    
  fclose(doc_file);
}

void aggregate_data_onefile(char* output_filename, char ** chunks_filenames, int nchunks) {
  
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
}



int main(int argc, char ** argv) {
  
  //------------------------------------ARGUMENTS HANDLING------------------------------------------------------------// 
  if (argc == 1) {
     fprintf(stderr, "%s: arguments required. Try '-h' or '--help'.\n", argv[0]);
     exit(EXIT_FAILURE);
  }
  
  program_name = argv[0];
  
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
  
  //------------------------------ Reading options-------------------------------// 
  do {
    next_option = getopt_long(argc, argv, short_options, long_options, NULL);
    
    switch (next_option) {
      case 'h':
	print_usage(stdout, 0);
      case 'p':
	parameters_filename = optarg;
	assert(optarg[0] != '-');
	printf("Parameters file: %s\n", parameters_filename);
	break;
      case 'i':
	input_filename = optarg;
	assert(optarg[0] != '-');
	printf("Input file with ray data: %s\n", input_filename);
	break;
      case 'o':
	output_filename = optarg;
	assert(optarg[0] != '-');
	printf("Output file with Radon data: %s\n", output_filename);
	break;
      case 'n':
	nthreads = atoi(optarg);
	assert(nthreads > 0);
	printf("Number of threads: %d\n", nthreads);
	break;
      case ':': /* missing argument */ 
	fprintf(stderr, "%s: option '-%c' requires an argument\n",
		argv[0], optopt);
	print_usage(stderr, EXIT_FAILURE);
      case '?':
	print_usage(stderr, EXIT_FAILURE);
      case -1: /* Done with the options. */
	break;
      default:
	abort();
    }
  } while (next_option != -1);
  
  if (optind < argc) {
    printf("Warning: non-option ARGV-arguments are omitted. See '-h' or '--help'.\n"); 
  }
  
  //-------------------------------------- INITIALIZATIONS ------------------------------------------------------------------------------/
  //---------------------- Read parameters of the ray/Radon grids from parameters file --------------------------------------------------/
  
  printf("  Reading parameters from %s...", parameters_filename);
  int *parameters;
  if (init_parameters_alloc(&parameters) != 0) { //allocate memory for parameters
     perror("Aborting : init_parameters_alloc.\n");
     exit(EXIT_FAILURE);
  }
  
  if (init_parameters(parameters_filename, parameters) != 0) {
     perror("Aborting : init_parameters.\n");
     exit(EXIT_FAILURE);
  }
  printf("Done.\n");
  
  int nshift = parameters[0];
  int nphi   = parameters[1];
  int ntheta = parameters[2];
  int ngrid  = parameters[3];		// number of z-slices and points on grid
  
  
  // --------------------------------------------Read ray transform data from file --------------------------------------------------------/
  
  
  printf("  Initializing data for ray transforms...");
  double*** values;
  if (init_ray_grid_alloc(nshift, nphi, ngrid, &values) != 0) {
     perror("Aborting : init_ray_grid_alloc.\n");
     exit(EXIT_FAILURE);
  }
  if (init_ray_values(input_filename, nshift, nphi, ngrid, values) != 0) {
     perror("Aborting : init_ray_values.\n");
     exit(EXIT_FAILURE);
  }
  printf("Done.\n");
  
  //---------------------------------------------------------------------------------------------------------------------------------------/
  
  //------------------------------------------------------ init grid in Radon space -------------------------------------------------------/
  double *shift, *phi, *theta;
  printf("  Initializing the grid in Radon space...");
  if (init_radon_grid_alloc(nshift, nphi, ntheta, &shift, &phi, &theta) != 0) {		// allocate memory for arrays : shift, theta, phi
     perror("Aborting : init_radon_grid_alloc.\n");
     exit(EXIT_FAILURE);
  }
  if (init_radon_grid(nshift, nphi, ntheta, shift, phi, theta) != 0) {		        // set values of : shift, theta, phi
     perror("Aborting : init_radon_grid.\n");
     exit(EXIT_FAILURE);
  }
  printf("Done.\n\n");
  
  //------------------------------------ END OF INITIALIZATIONS ------------------------------------------------------------------------/
  
  //------------------------ START REDUCTION RAY TRANSFORMS -> RADON TRANSFORMS (PARALLELIZATION VIA OPENMP) ---------------------------/
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

        //set equally shifts of planes for each thread
	int block_size = nshift / nthreads + (thread_ID < (nshift % nthreads) ? 1 : 0);
	
        int shift_min = (thread_ID < (nshift % nthreads) ? (nshift / nthreads + 1)*thread_ID : (nshift / nthreads)*thread_ID + (nshift % nthreads)),
	    shift_max = shift_min + block_size;
	    
        printf("    Thread %d. Domain of work (shifts): start %d, end %d, size %d\n", thread_ID, shift_min, shift_max-1, block_size);
        
	
        //open ouptut file
        FILE *foutput;
	foutput = fopen(chunks_filenames[thread_ID], "w");
	
	int i_shift, i_phi, i_theta;
    
        //ray integrals in the given OZ slice
        for (i_shift = shift_min; i_shift < shift_max; ++i_shift) {
      
	    //iterate ray integrals over shifts in slice plane 
	    for (i_phi = 0; i_phi < nphi; ++i_phi) {
	      
		//iterate ray integrals over polar angle in slice plane
		for (i_theta = 0; i_theta < ntheta; ++i_theta) {	
	            
		    //reduction of ray integrals to Radon transform over the plane (shift[i_shift], phi[i_phi], theta[i_theta])
		    double Rt = Rt_slicing(i_shift, i_phi, i_theta,
				  shift, phi, theta,
				  nshift, nphi, ntheta, 
				  values, 
				  ngrid);
		    //write result to file
		    fprintf(foutput, "%lf, %lf, %lf, %lf\n", shift[i_shift], phi[i_phi], theta[i_theta], Rt);
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
  
  //------------------------ END REDUCTION RAY TRANSFORMS -> RADON TRANSFORMS (PARALLELIZATION VIA OPENMP) ---------------------------/
  
  
  //------------------------------------------ START AGGREGATING DATA TO ONE FILE ----------------------------------------------------/
  
  generate_data_readme(output_filename, nshift, nphi, ntheta);
  printf("  Aggregating chunks to one output file...");
  aggregate_data_onefile(output_filename, chunks_filenames, nthreads);
  printf("Done.\n");
  
  //---------------------------------------------- END AGGREGATING DATA TO ONE FILE --------------------------------------------------/
  

  //------------------------------------------------------- START CLEANING MEMORY ----------------------------------------------------/
  
  printf("  Cleaning allocated memory...");
  if (clean_chunks_files(nthreads, chunks_filenames) != 0) {
     fprintf(stderr, "Warning : clean_chunks_files.\n");
  }
  if (clean_chunks_filenames(nthreads, chunks_filenames) != 0) {
     fprintf(stderr, "Warning : clean_chunks_filenames.\n");
  }
 
  if (clean_parameters(parameters) != 0) {
     fprintf(stderr, "Warning : clean_parameters.\n");
  }
  if (clean_ray_values(nshift, nphi, ngrid, values) != 0) {
     fprintf(stderr, "Warning : clean_ray_values.\n");
  }
  if (clean_radon_grid(shift, phi, theta) != 0) {
     fprintf(stderr, "Warning : clean_radon_grid.\n");
  }
  printf("Done.\n");
  
  //----------------------------------------------------------- END CLEANING MEMORY --------------------------------------------------/
  
  exit(EXIT_SUCCESS);
}