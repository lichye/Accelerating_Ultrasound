#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <stdint.h>
#include <pthread.h>

#define numThreads 1024
#define trans_x  32 // Transducers in x dim
#define trans_y  32 // Transducers in y dim	
#define data_len 12308 // Number for pre-processed data values per channel
#define tx_x  0 // Transmit transducer x position
#define tx_y  0 // Transmit transducer y position
#define tx_z -0.001 // Transmit transducer z position
const float idx_const = 0.000009625; // Speed of sound and sampling rate, converts dist to index
const int filter_delay = 140; // Constant added to index to account filter delay (off by 1 from MATLAB)
float rx_z = 0; // Receive transducer z position
#define ALIGN 64

/*DATA*/
float *point_x; // Point x position
float *point_y; // Point y position
float *point_z; // Point z position
float *image_pos; // Pointer to current position in image
float *image;  // Pointer to full image (accumulated so far)
float *dist_tx; // Transmit distance (ie first leg only)
FILE* input;
FILE* output;
float *rx_x; // Receive transducer x position
float *rx_y; // Receive transducer y position
float *rx_data; // Pointer to pre-processed receive channel data


struct threadArgs{
	int start;
	int end;
};

void *reflec_dis_compute(void *args) {
	int offset = 0;
	int it_rx;
	float x_comp; 
	float y_comp; 
	float z_comp; 
	float dist; // Full distance
	int index_t;
	int it_r;
	int point;

	struct threadArgs * range = (struct threadArgs *) args;

	for (point = range->start; point < range->end; point++) {
		x_comp = tx_x - point_x[point];
		x_comp = x_comp * x_comp;
		y_comp = tx_y - point_y[point];
		y_comp = y_comp * y_comp;
		z_comp = tx_z - point_z[point];
		z_comp = z_comp * z_comp;
		dist_tx[point] = (float)sqrt(x_comp + y_comp + z_comp);
	}

	for (it_rx = 0; it_rx < trans_x * trans_y; it_rx++) {
		image_pos = image;
		for (it_r = range->start; it_r < range->end; it_r++) {
			x_comp = rx_x[it_rx] - point_x[it_r];
			x_comp = x_comp * x_comp;
			y_comp = rx_y[it_rx] - point_y[it_r];
			y_comp = y_comp * y_comp;
			z_comp = rx_z - point_z[it_r];
			z_comp = z_comp * z_comp;
			dist = dist_tx[it_r] + (float)sqrt(x_comp + y_comp + z_comp);
			index_t = (int)(dist/idx_const + filter_delay + 0.5);
			image_pos[it_r] += rx_data[index_t+offset];	
		}
		offset += data_len;
	}
}

int main (int argc, char **argv) {
	int i;
	int size = atoi(argv[1]);

	int pts_r = 1560; // Radial points along scanline
	int sls_t = size; // Number of scanlines in theta
	int sls_p = size; // Number of scanlines in phi
	/* Allocate space for data */
	rx_x = (float*) malloc(trans_x * trans_y * sizeof(float));
	if (rx_x == NULL) fprintf(stderr, "Bad malloc on rx_x\n");
	rx_y = (float*) malloc(trans_x * trans_y * sizeof(float));
	if (rx_y == NULL) fprintf(stderr, "Bad malloc on rx_y\n");
	
	//posix_memalign((void*)&rx_data, ALIGN, data_len * trans_x * trans_y * sizeof(float));
	rx_data = (float*) malloc(data_len * trans_x * trans_y * sizeof(float));
	if (rx_data == NULL) fprintf(stderr, "Bad malloc on rx_data\n");

	posix_memalign((void*)&point_x, ALIGN, pts_r * sls_t * sls_p * sizeof(float));
	if (point_x == NULL) fprintf(stderr, "Bad malloc on point_x\n");
	
	posix_memalign((void*)&point_y, ALIGN, pts_r * sls_t * sls_p * sizeof(float));
	if (point_y == NULL) fprintf(stderr, "Bad malloc on point_y\n");
	
	posix_memalign((void*)&point_z, ALIGN, pts_r * sls_t * sls_p * sizeof(float));
	if (point_z == NULL) fprintf(stderr, "Bad malloc on point_z\n");

	posix_memalign((void*)&dist_tx, ALIGN, pts_r * sls_t * sls_p * sizeof(float));
	if (dist_tx == NULL) fprintf(stderr, "Bad malloc on dist_tx\n");

	posix_memalign((void*)&image, ALIGN, pts_r * sls_t * sls_p * sizeof(float));
	if (image == NULL) fprintf(stderr, "Bad malloc on image\n");

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






	/* get start timestamp */
 	struct timeval tv;
    	gettimeofday(&tv,NULL);
    	uint64_t start = tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
 
	/* --------------------------- COMPUTATION ------------------------------ */
	/* First compute transmit distance */
	float x_comp; 
	float y_comp; 
	float z_comp; 

	/* Variables for image space points */
	int point = 0;
	int workSpace = sls_t * sls_p * pts_r;

	struct threadArgs work_ranges[numThreads];
	pthread_t child_threads[numThreads];

	int workLength = workSpace/numThreads;
	if(workSpace%workLength!=0){
		printf("There is left data!\n");
	}
	for (i = 0; i < numThreads; i++) {
		work_ranges[i].start = i*workLength;
		work_ranges[i].end = (i+1)*workLength;
	}
	for( i = 0; i < numThreads; i++) {
        	pthread_create(&child_threads[i], NULL, reflec_dis_compute, &work_ranges[i] );
    	}
    	for( i = 0; i < numThreads; i++) {
        	pthread_join(child_threads[i], NULL);
    	}	


	/* --------------------------------------------------------------------- */

	/* get elapsed time */
    	gettimeofday(&tv,NULL);
    	uint64_t end = tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
    	uint64_t elapsed = end - start;
	printf("Processing complete. Preparing output.\nThe elapsed time is: %lu\n", elapsed);
	fflush(stdout);

	/* Write result to file */
	char* out_filename;
        #ifdef __MIC__
	  out_filename = "/home/micuser/beamforming_output.bin";
	  //out_filename = "beamforming_output_mic.bin";

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
