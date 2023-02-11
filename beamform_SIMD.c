// 3D Ultrasound beamforming baseline code for EECS 570 
// Created by: Richard Sampson, Amlan Nayak, Thomas F. Wenisch
// Revision 1.0 - 11/15/16

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include <pthread.h>
#include <immintrin.h>
#pragma intel optimization_parameter target_arch=AVX512
#pragma intel optimization_parameter target_arch=KNC

#define ALIGN 64
#define N 16

#define PARTITION 16
#define STEP 4


// ===========================================================
  
  int size;
/* Variables for transducer geometry */
	int trans_x = 32; // Transducers in x dim
	int trans_y = 32; // Transducers in y dim
	
	float *rx_x; // Receive transducer x position
	float *rx_y; // Receive transducer y position
	float rx_z = 0; // Receive transducer z position

	int data_len = 12308; // Number for pre-processed data values per channel
	int offset = 0; // Offset into rx_data
	float *rx_data; // Pointer to pre-processed receive channel data

	float tx_x = 0; // Transmit transducer x position
	float tx_y = 0; // Transmit transducer y position
	float tx_z = -0.001; // Transmit transducer z position

	/* Variables for image space points */
	int point; // Index into image space

	float *point_x; // Point x position
	float *point_y; // Point y position
	float *point_z; // Point z position

	int pts_r = 1560; // Radial points along scanline

	float *image_pos; // Pointer to current position in image
	float *image;  // Pointer to full image (accumulated so far)
  
  float *dist_tx; // Transmit distance (ie first leg only)
  
  const float idx_const = 0.000009625; // Speed of sound and sampling rate, converts dist to index
	const int filter_delay = 140; // Constant added to index to account filter delay (off by 1 from MATLAB)
 
  int rx_range; // = trans_x * trans_y;
  int range;

  // -------------------------Something new!!!-------------------------
  float rcp_idx_const;
  const float new_filter_delay = 140.5;
  __m512 ridx_vec; // = _mm512_set1_ps (rcp_idx_const);
  __m512 filter_vec; // = _mm512_set1_ps (new_filter_delay);
  
  
struct thread_args {
  int start;
  
  __m512 xx;// = _mm512_set1_ps (0.0);  // x intermediate value
  __m512 yy;// = _mm512_set1_ps (0.0);  // y intermediate value
  
  __m512i ii; // INDEX vector
  __m512i offset_vec; // offset vector
  
  __m512 rx_x_vec;
  __m512 rx_y_vec;
  __m512 rx_z_vec;
};

// ------------------------------------------------------------

void *calcDistance (void * _arg){//void *_arg) {
  
  struct timeval tv;
  uint64_t start_time;

  gettimeofday(&tv,NULL);
    	start_time = tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
    	
  //struct thread_args * arg = (struct thread_args *) _arg;
  
  struct thread_args * arg = (struct thread_args *) _arg;
  
  //int start = *((int *) _start);
  int start = arg->start;
  int end = arg->start + range;
  
  
  

  //float x_comp, y_comp, z_comp;
  //float dist;
  int it_rx;
  int point, index;
  int offset = 0;
  float *image_pos;// = image; // Pointer to current position in image
  float *image_start = &(image[start]);
  
  // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
 
  
  
  int i;
     
  // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  
  
	for (it_rx = 0; it_rx < rx_range; it_rx++) {

		image_pos = image_start; // Reset image pointer back to beginning
   
    arg->rx_x_vec = _mm512_set1_ps (rx_x[it_rx]);
    arg->rx_y_vec = _mm512_set1_ps (rx_y[it_rx]);
    
    arg->offset_vec = _mm512_set1_epi32 (offset);

    /*
		// Iterate over entire image space
    for (point = start; point < end; point++) {
    
					x_comp = rx_x[it_rx] - point_x[point];
					x_comp = x_comp * x_comp;
					y_comp = rx_y[it_rx] - point_y[point];
					y_comp = y_comp * y_comp;
					z_comp = rx_z - point_z[point];
					z_comp = z_comp * z_comp;

					dist = dist_tx[point] + (float)sqrt(x_comp + y_comp + z_comp);
					index = (int)(dist*rcp_idx_const + new_filter_delay);
                
          *image_pos++ += rx_data[index+offset];
                                             
		}*/
    //px_ptr = &point_x[start];
    //py_ptr = &point_y[start];
    //pz_ptr = &point_z[start];
    
     
    for (point = start; point < end; point+=N) {
      
      // SUBTRACTION: x = rx_x[it_rx] - point_x[point];
      arg->xx = _mm512_sub_ps (arg->rx_x_vec, *(__m512*)(point_x+point));
      
      // MULTIPLICATION: x = x * x;
      arg->xx =  _mm512_mul_ps (arg->xx, arg->xx);
      
      // SUBTRACTION: y = rx_y[it_rx] - point_y[point];
      arg->yy = _mm512_sub_ps (arg->rx_y_vec, *(__m512*)(point_y+point));
      
      // x = x + y*y
      arg->xx = _mm512_fmadd_ps (arg->yy, arg->yy, arg->xx);
      
      // SUBTRACTION: z = rx_z - point_z[point];
      arg->yy = _mm512_sub_ps (arg->rx_z_vec, *(__m512*)(point_z+point));
      
      // x = x + z*z
      arg->xx = _mm512_fmadd_ps (arg->yy, arg->yy, arg->xx);
      
      // SQRT of sum: x = sqrt(x)
      arg->xx = _mm512_sqrt_ps (arg->xx);
      
      // TOTAL DISTANCE: dist = dist_tx[point] + x
      arg->xx = _mm512_add_ps (arg->xx, *(__m512*)(dist_tx+point));
      
      // CALCULATE INDEX (float)
      arg->xx = _mm512_fmadd_ps (arg->xx, ridx_vec, filter_vec);
      
      // Convert to INT
      arg->ii = _mm512_cvtfxpnt_round_adjustps_epi32 (arg->xx, (_MM_FROUND_TO_NEG_INF |_MM_FROUND_NO_EXC), _MM_EXPADJ_NONE );
      //ii = _mm512_cvtps_epi32 (xx);
      
      arg->ii = _mm512_add_epi32 (arg->ii, arg->offset_vec);  // ADD offset
      
      arg->xx = _mm512_i32gather_ps (arg->ii, rx_data, 4);  // GATHER memory
      //xx = _mm512_i32gather_ps (rx_data, ii, 4);  // GATHER memory
      
      //oo = _mm512_load_ps (image_pos);// Load the image value and add new value to it
      
      //xx = _mm512_add_ps (xx, oo);  // ADD with gathered values
      arg->xx = _mm512_add_ps (arg->xx, *(__m512*)(image_pos));
      
      //_mm512_stream_ps (image_pos, xx);  // STORE back to image using non-temporal memory hint
      _mm512_store_ps (image_pos, arg->xx);  // STORE back to image
      
      image_pos += N;
      
    }
    
		offset += data_len;
	} 
  uint64_t end_time;
  uint64_t elapsed;
  
  gettimeofday(&tv,NULL);
    	end_time = tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
    	elapsed = end_time - start_time;

	printf("THREAD loop time (usec): %lld\n", elapsed);
  return 0;
}

// ===========================================================

int main (int argc, char **argv) {

	size = atoi(argv[1]);

	int sls_t = size; // Number of scanlines in theta
	int sls_p = size; // Number of scanlines in phi

	/* Iterators */
	int it_rx; // Iterator for recieve transducer
	int it_r; // Iterator for r
	int it_t; // Iterator for theta
	int it_p; // Iterator for phi

	/* Variables for distance calculation and index conversion */
	float x_comp; // Itermediate value for dist calc
	float y_comp; // Itermediate value for dist calc
	float z_comp; // Itermediate value for dist calc
  
  rx_range = trans_x * trans_y;
  
  rcp_idx_const = 1/idx_const;
  
  ridx_vec = _mm512_set1_ps (rcp_idx_const);
  filter_vec = _mm512_set1_ps (new_filter_delay);
	
	float dist; // Full distance
	
	int index; // Index into transducer data

        FILE* input;
        FILE* output;

	/* Allocate space for data */
	rx_x = (float*) malloc(trans_x * trans_y * sizeof(float));
	if (rx_x == NULL) fprintf(stderr, "Bad malloc on rx_x\n");
	rx_y = (float*) malloc(trans_x * trans_y * sizeof(float));
	if (rx_y == NULL) fprintf(stderr, "Bad malloc on rx_y\n");
	rx_data = (float*) malloc(data_len * trans_x * trans_y * sizeof(float));
	if (rx_data == NULL) fprintf(stderr, "Bad malloc on rx_data\n");
  
  /*
	point_x = (float *) malloc(pts_r * sls_t * sls_p * sizeof(float));
	if (point_x == NULL) fprintf(stderr, "Bad malloc on point_x\n");
	point_y = (float *) malloc(pts_r * sls_t * sls_p * sizeof(float));
	if (point_y == NULL) fprintf(stderr, "Bad malloc on point_y\n");
	point_z = (float *) malloc(pts_r * sls_t * sls_p * sizeof(float));
	if (point_z == NULL) fprintf(stderr, "Bad malloc on point_z\n");
  */
  
  posix_memalign((void*)&point_x, ALIGN, pts_r * sls_t * sls_p * sizeof(float));
  posix_memalign((void*)&point_y, ALIGN, pts_r * sls_t * sls_p * sizeof(float));
  posix_memalign((void*)&point_z, ALIGN, pts_r * sls_t * sls_p * sizeof(float));
  
  posix_memalign((void*)&dist_tx, ALIGN, pts_r * sls_t * sls_p * sizeof(float));
	//dist_tx = (float*) malloc(pts_r * sls_t * sls_p * sizeof(float));
	//if (dist_tx == NULL) fprintf(stderr, "Bad malloc on dist_tx\n");
  
  posix_memalign((void*)&image, ALIGN, pts_r * sls_t * sls_p * sizeof(float));
	//image = (float *) malloc(pts_r * sls_t * sls_p * sizeof(float));
	//if (image == NULL) fprintf(stderr, "Bad malloc on image\n");
  
 
	memset(image, 0, pts_r * sls_t * sls_p * sizeof(float));

	/* validate command line parameter */
	if (argc < 1 || !(strcmp(argv[1],"16") || strcmp(argv[1],"32") || strcmp(argv[1],"64"))) {
	  printf("Usage: %s {16|32|64}\n",argv[0]);
	  fflush(stdout);
	  exit(-1);
	}

	char buff[128];
        #ifdef __MIC__
	  sprintf(buff, "/beamforming_input_%s.bin", argv[1]);
        #else // !__MIC__
	  sprintf(buff, "/n/typhon/data1/home/eecs570/beamforming_input_%s.bin", argv[1]);
        #endif

        input = fopen(buff,"rb");
	if (!input) {
	  printf("Unable to open input file %s.\n", buff);
	  fflush(stdout);
	  exit(-1);
	}	

	/* Load data from binary */
	fread(rx_x, sizeof(float), trans_x * trans_y, input); 
	fread(rx_y, sizeof(float), trans_x * trans_y, input); 

	fread(point_x, sizeof(float), pts_r * sls_t * sls_p, input); 
	fread(point_y, sizeof(float), pts_r * sls_t * sls_p, input); 
	fread(point_z, sizeof(float), pts_r * sls_t * sls_p, input); 

	fread(rx_data, sizeof(float), data_len * trans_x * trans_y, input); 
        fclose(input);

	printf("Beginning computation\n");
	fflush(stdout);
 
  int NUM_THREADS = PARTITION * PARTITION * STEP;  // =1024
  
  int partition_size = size / PARTITION;     // =1 when size=16, =2 when size=32, =4 when size=64
  int step_size = pts_r / STEP;              // =390 if STEP=4
  
  range = partition_size * partition_size * step_size;
  
  //printf("range = %d\n", range);

  pthread_t thread_array[NUM_THREADS]; // Thread array
  
  struct thread_args arg_array[NUM_THREADS]; // Argument array
  
  int start_array[NUM_THREADS];    // Argument array
  int i;
  for (i = 0; i < NUM_THREADS; i++) {
    start_array[i] = i * range;
    arg_array[i].start = i * range;
    arg_array[i].rx_z_vec = _mm512_set1_ps (rx_z);
  }
  
  struct thread_args * arg_ptr = arg_array;

	/* get start timestamp */
 	struct timeval tv;
    	gettimeofday(&tv,NULL);
    	uint64_t start = tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
 
	/* --------------------------- COMPUTATION ------------------------------ */
	/* First compute transmit distance */
  point = 0;
	for (it_t = 0; it_t < sls_t; it_t++) {

		for (it_p = 0; it_p < sls_p; it_p++) {
			for (it_r = 0; it_r < pts_r; it_r++) {

				x_comp = tx_x - point_x[point];
				x_comp = x_comp * x_comp;
				y_comp = tx_y - point_y[point];
				y_comp = y_comp * y_comp;
				z_comp = tx_z - point_z[point];
				z_comp = z_comp * z_comp;

				dist_tx[point++] = (float)sqrt(x_comp + y_comp + z_comp);
      }
    }
	}

  /* Now compute reflected distance, find index values, add to image */
  for(i = 0; i < NUM_THREADS; i++) {
    pthread_create(&thread_array[i], NULL, calcDistance, arg_ptr++);
  }
  
  //gettimeofday(&tv,NULL);
  //  	uint64_t end = tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
  //  	uint64_t elapsed = end - start;
  
  for(i = 0; i < NUM_THREADS; i++) {
    pthread_join(thread_array[i], NULL);
  }


	/* --------------------------------------------------------------------- */

	/* get elapsed time */
    	gettimeofday(&tv,NULL);
    	uint64_t end = tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
    	uint64_t elapsed = end - start;

	printf("@@@ Elapsed time (usec): %lld\n", elapsed);
	printf("Processing complete.  Preparing output.\n");
	fflush(stdout);

	/* Write result to file */
	char* out_filename;
        #ifdef __MIC__
	  out_filename = "/home/micuser/beamforming_output.bin";
        #else // !__MIC__
	  out_filename = "beamforming_output.bin";
        #endif
        output = fopen(out_filename,"wb");
	fwrite(image, sizeof(float), pts_r * sls_t * sls_p, output); 
	fclose(output);

	printf("Output complete.\n");
	fflush(stdout);

	/* Cleanup */
	free(rx_x);
	free(rx_y);
	free(rx_data);
	free(point_x);
	free(point_y);
	free(point_z);
	free(dist_tx);
	free(image);

	return 0;
}
